nextflow.preview.dsl=2

/* Process UK Biobank data release into a DuckDB database */

params.enc = "ukb12345.enc"
params.key = "k1234r12345.key"
params.showcase = "https://biobank.ndph.ox.ac.uk/~bbdatan/Data_Dictionary_Showcase.tsv"
params.encoding = "https://biobank.ndph.ox.ac.uk/showcase/ukb/utilx/encoding.dat"

workflow {
    // packed/encrypted UKB release file
    ENC_CH = Channel
        .fromPath(params.enc)
    
    // decryption key
    KEY_CH = Channel
        .fromPath(params.key)
        
    // stage the files
    ENC_KEY_CH = STAGE(ENC_CH, KEY_CH)
    
    // unpack the file
    UKB_ENC_CH = UNPACK(ENC_KEY_CH)

    // data dictionary
    UKB_DOCS_CH = DOCS(UKB_ENC_CH)

     /* Download and process the showcase data dictionary */
    SHOWCASE_CH = Channel
        .fromPath(params.showcase)

    SHOWCASE_FIELDS_CH = DICTIONARY(SHOWCASE_CH)

    // Parse the list of ALL fields and tables (categories)
    TABLES_FIELDS_CH = SHOWCASE_FIELDS_CH
        .flatten()
        .first()
        .splitCsv(header: true)
        .map {['FieldID': it.FieldID, 'Table': it.Table]}

    DICTIONARY_CH = SHOWCASE_FIELDS_CH.flatten().last()

    /* process the fields contained in the data release */
    FIELDS_CH = UKB_DOCS_CH
        .flatten()
        .last()
        .splitCsv(header: ['FieldID'])

    // filter the tables/fields list to fields that are in the data release
    TABLES_FIELDS_INCLUDE_CH = TABLES_FIELDS_CH
        .cross(FIELDS_CH)
        .map {[it[0].Table, it[0].FieldID]}
        .groupTuple()

    // write out include list fields list for each table to extract
    INCLUDE_CH = INCLUDE(TABLES_FIELDS_INCLUDE_CH)
    
    // Download the encodings file
    ENCODING_CH = Channel
        .fromPath(params.encoding)
    
    // convert to tab and R code
    R_CH = CONVERT(UKB_ENC_CH, ENCODING_CH, INCLUDE_CH)
    
    // source the R code and write out the data.frame as RDS
    RDS_CH = RDS(R_CH)
    
    // copy to duckdb
    DUCK_CH = DUCK(RDS_CH)
    
}

/* Stage release files */
process STAGE {
    tag "Staging data"
    queue "staging"
    
    input:
    path enc
    path key
    
    output:
    tuple path('ukb.enc'), path('ukb.key')
    
    script:
    """
    cp $enc ukb.enc
    cp $key ukb.key
    """
}

/* Unpack the data file */
process UNPACK {
    tag "Unpacking data"
    
    cpus = 1
    memory = 4.GB
    time = '30m'
    
    input:
    tuple path(enc), path(key)
    
    output:
    path('ukb.enc_ukb')
    
    script:
    """
    ukbunpack $enc $key
    """
}

/* Convert the data file to CSV */
process DOCS {
    tag "Extracting HTML dictionary"
    
    cpus = 1
    memory = 2.GB
    time = '6h'

    scratch true
    stageInMode 'copy'
    stageOutMode 'copy'
    
    input:
    path(enc_ukb)
    
    output:
    tuple path('ukb.html'), path('fields.ukb')
    
    script:
    """
    dconvert $enc_ukb docs
    """
}


/* Parse categories and fields from showcase dictionary */
process DICTIONARY {
    tag "Categorising fields"
    
    executor 'local'
    module 'igmm/apps/R/4.1.0'
    
    input:
    path(showcase)
    
    output:
    tuple path("Tables_Fields.csv"), path("Dictionary.tsv")
    
    script:
    """
    #!/usr/bin/env Rscript
    library(dplyr)
    library(readr)
    library(tidyr)
    library(snakecase)

    # Convert showcase into a table of categories and fields

    showcase <- read_tsv("$showcase")

    showcase_fields <-
    showcase |>
        separate(Path, sep=' > ', into=c('Cat1', 'Cat2'), extra='drop') |>
        mutate(Table=to_any_case(Cat2, case='big_camel', sep_in='[^A-Za-z]')) |>
        select(-Cat1, -Cat2) |>
        select(Table, FieldID, Field, everything())

    write_tsv(showcase_fields, 'Dictionary.tsv')

    write_csv(select(showcase_fields, Table, FieldID), 'Tables_Fields.csv')
    """
}


/* Write out a list of fields for an include list to convert */
process INCLUDE {
    tag "${table}"
    
    executor 'local'

    input:
    tuple val(table), val(fields)

    output:
    path("${table}.fields")

    script:
    """
    echo "${fields.join('\n')}" > ${table}.fields
    """

}

/* Convert included fields in the data file to R */
process CONVERT {
    tag "Converting ${include.simpleName}"
    
    cpus = 1
    memory = 2.GB
    time = '6h'

    scratch true
    stageInMode 'copy'
    stageOutMode 'copy'
    
    // ignore dconvert exiting with status 1 (caused by "ROSETTA Error"),
    // which can be safely ignored according to UKB Access
    errorStrategy 'ignore'
    
    input:
    path(enc_ukb)
    path(encoding)
    each path(include)
    
    output:
    tuple path("${include.simpleName}.r"), path("${include.simpleName}.tab")
    
    script:
    """
    dconvert $enc_ukb r -o${include.simpleName} -e${encoding} -i${include}
    """
}

/* Process tab data files using the associated R script 
   Write out as .rds */
process RDS {
    tag "Running ${rsource.simpleName}"
    
    module 'igmm/apps/R/4.1.0'
    
    cpus = 3
    memory = 48.GB
    time = '2h'
    
    scratch true
    stageInMode 'copy'
    stageOutMode 'copy'
    
    input:
    tuple path(rsource), path(tab)
    
    output:
    path("${rsource.simpleName}.rds")
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # load the data file
    bd <- read.table("${tab}", header=TRUE, sep="\\t")
    
    # source the R script, but skip line3 that loads the data file
    # with a hardcoded full path
    rsource_lines <- readLines("${rsource}", warn=FALSE)
    rsource_exprs <- str2expression(rsource_lines[-3])
    source(exprs=rsource_exprs)
    
    # write out the data.frame an RDS file
    saveRDS(bd, "${rsource.simpleName}.rds")
    """
}

/* Copy the RDS to DuckDB */
process DUCK {
    tag "DuckDB ${rds.simpleName}"
    
    module 'igmm/apps/R/4.1.0'
    
    cpus = 2
    memory = 24.GB
    time = '30m'

    scratch true
    stageInMode 'copy'
    stageOutMode 'copy'
    
    input:
    path(rds)
    
    output:
    path("${rds.simpleName}.duckdb")
    
    script:
    """
    #!/usr/bin/env Rscript
    library(DBI)
    
    con = dbConnect(duckdb::duckdb(), dbdir="${rds.simpleName}.duckdb", read_only=FALSE)
    
    Table <- readRDS("${rds}")
    
    dbWriteTable(con, "${rds.simpleName}", Table)
    dbDisconnect(con, shutdown=TRUE)
    """
}
