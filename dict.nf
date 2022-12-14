nextflow.enable.dsl=2

/* Process fields based on data dictionary categories */

params.showcase = "https://biobank.ndph.ox.ac.uk/~bbdatan/Data_Dictionary_Showcase.tsv"

workflow {
    SHOWCASE_URL_CH = Channel
        .value(params.showcase)

    SHOWCASE_CH = SHOWCASE(SHOWCASE_URL_CH)

    TABLE_FIELDS_CH = DICTIONARY(SHOWCASE_CH)

}

/* Download showcase data dictionary */
process SHOWCASE {
    tag "Downloading dictionary"
    
    input:
    val url
    
    output:
    path("Data_Dictionary_Showcase.tsv")
    
    script:
    """
    wget $url
    """
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