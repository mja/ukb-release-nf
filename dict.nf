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

    DICTIONARY_CH = SHOWCASE_FIELDS_CH.flatten().last()

    /* process the fields contained in the data release */
    FIELDS_CH = Channel
        .fromPath(params.fields)
        .splitText()

    // example of filtering by field
    TABLES_FIELDS_CH
        .filter {it.FieldID == '3'}
        .view()

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