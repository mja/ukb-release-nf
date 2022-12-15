nextflow.enable.dsl=2

/* Process fields based on data dictionary categories */

params.showcase = "https://biobank.ndph.ox.ac.uk/~bbdatan/Data_Dictionary_Showcase.tsv"
params.fields = "fields2.ukb"

workflow {
    /* Download and process the showcase data dictionary */
    SHOWCASE_CH = Channel
        .fromPath(params.showcase)

    SHOWCASE_FIELDS_CH = DICTIONARY(SHOWCASE_CH)

    TABLES_FIELDS_CH = SHOWCASE_FIELDS_CH
        .flatten()
        .first()
        .splitCsv(header: true)
        .map {['FieldID': it.FieldID, 'Table': it.Table]}

    DICTIONARY_CH = SHOWCASE_FIELDS_CH.flatten().last()

    /* process the fields contained in the data release */
    FIELDS_CH = Channel
        .fromPath(params.fields)
        .splitCsv(header: ['FieldID'])

    // filter the tables/fields list to fields that are in the data release
    TABLES_FIELDS_INCLUDE_CH = TABLES_FIELDS_CH
        .cross(FIELDS_CH)
        .map {[it[0].Table, it[0].FieldID]}
        .groupTuple()

    // write out include list fields list for each table to extract
    INCLUDE_CH = INCLUDE(TABLES_FIELDS_INCLUDE_CH)

}

/* Parse categories and fields from showcase dictionary */
process DICTIONARY {
    tag "Categorising fields"
    
    input:
    path(showcase)
    
    output:
    tuple path("Tables_Fields.csv"), path("Dictionary.tsv")
    
    script:
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
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

    input:
    tuple val(table), val(fields)

    output:
    tuple val(table), path("${table}.fields")

    script:
    """
    echo "${fields.join('\n')}" > ${table}.fields
    """

}