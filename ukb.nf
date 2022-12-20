nextflow.preview.dsl=2

/* Process UK Biobank data release into a DuckDB database */

params.enc = "ukb12345.enc"
params.key = "k1234r12345.key"
params.showcase = "https://biobank.ndph.ox.ac.uk/~bbdatan/Data_Dictionary_Showcase.tsv"
params.encoding = "https://biobank.ndph.ox.ac.uk/ukb/ukb/utilx/encoding.dat"

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
    
    // publish the data dictionary
    PUB(UKB_DOCS_CH, ENC_CH)

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
    // also pass in the dictionary for additional data munging
    R_DICTIONARY_CH = R_CH
        .combine(DICTIONARY_CH)
    RDS_CH = RDS(R_DICTIONARY_CH)
    
    // collect RDS files into a single list
    RDS_COLLECT_CH = RDS_CH
        .collect()
    
    // copy to duckdb
    DUCK_CH = DUCK(RDS_COLLECT_CH, DICTIONARY_CH, ENC_CH)
    
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

/* Publish the HTML data dictionary */
process PUB {
    tag "Publish dictionary"
    
    publishDir 'release', mode: 'copy'
    
    executor 'local'
    
    input:
    tuple path('ukb.html'), path('fields.ukb')
    path(enc)
    
    output:
    path("${enc.simpleName}.html")
    
    script:
    """
    cp ukb.html "${enc.simpleName}.html"
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
    dconvert $enc_ukb r -o${include.simpleName} -e${encoding} -i${include} || true
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
    tuple path(rsource), path(tab), path(dictionary)
    
    output:
    path("${rsource.simpleName}.rds")
    
    script:
    """
    #!/usr/bin/env Rscript
    library(readr)
    library(dplyr)
    library(stringr)
    library(lubridate)
    
    # load the data file
    bd <- read.table("${tab}", header=TRUE, sep="\\t")
    
    # load the dictionary
    dictionary <- read_tsv("${dictionary}")
    
    # source the R script, but skip line3 that loads the data file
    # with a hardcoded full path
    rsource_lines <- readLines("${rsource}", warn=FALSE)
    rsource_exprs <- str2expression(rsource_lines[-3])
    source(exprs=rsource_exprs)
    
    # find fields that are bulk text identifiers and convert them to factors
    bulk_field_ids <- dictionary |> filter(ItemType == "Bulk") |> mutate(FieldID=paste0("f.", FieldID, ".")) |> pull(FieldID)
    bulk_fields <- paste0("f.", bulk_field_ids,  ".")
    bulk_field_list <- unlist(lapply(bulk_fields, function(field) str_subset(colnames(bd), fixed(field))))
    
    # find fields that are times and convert them to timestamps
    time_field_ids <- dictionary |> filter(ValueType == "Time") |> pull(FieldID)
    time_fields <- paste0("f.", time_field_ids,  ".")
    time_field_list <- unlist(lapply(time_fields, function(field) str_subset(colnames(bd), fixed(field))))
    
    # only process columns that are still character
    bulk_fields_mu <- bd |> select(one_of(bulk_field_list)) |> select(where(is.character)) |> colnames()
    time_fields_mu <- bd |> select(one_of(time_field_list)) |> select(where(is.character)) |> colnames()
    
    bd <- bd |> mutate(across(one_of(bulk_fields_mu), factor)) |>
                mutate(across(one_of(time_fields_mu), ymd_hms))
    
    # write out the data.frame an RDS file
    saveRDS(bd, "${rsource.simpleName}.rds")
    """
}

/* Copy the RDSs to DuckDB */
process DUCK {
    tag "DuckDB"
    
    publishDir 'release', mode: 'copy'
    module 'igmm/apps/R/4.1.0'
    
    cpus = 4
    memory = 64.GB
    time = '1h'

    scratch true
    stageInMode 'symlink'
    stageOutMode 'copy'
    
    input:
    path(rds_collection)
    path(dictionary)
    path(enc)
    
    output:
    path("${enc.simpleName}.duckdb")
    
    script:
    """
    #!/usr/bin/env Rscript
    
    library(DBI)
    library(dbplyr)
    library(readr)
    library(stringr)
    
    # parse out list of RDS files containing table data
    rds_paths <- "${rds_collection}"
    rds_path_list <- str_split(rds_paths, pattern=' ')[[1]]
    
    # open the database collection
    con = dbConnect(duckdb::duckdb(), dbdir="${enc.simpleName}.duckdb", read_only=FALSE)
    
    # load the dictionary
    dictionary <- read_tsv("${dictionary}")
    
    dbWriteTable(con, 'Dictionary', dictionary)
    
    # process each table
    for(rds in rds_path_list) {
        Table <- readRDS(rds)
        TableName <- str_split(rds, pattern='\\\\.')[[1]][1]
        dbWriteTable(con, TableName, Table)
        rm(Table); gc()
        index_query <- build_sql("CREATE UNIQUE INDEX ", ident(TableName), "_feid_idx ON ", ident(TableName), " (\\"f.eid\\")", con = con)
        dbExecute(con, index_query)
    }
    
    dbDisconnect(con, shutdown=TRUE)
    """
}
