nextflow.preview.dsl=2

/* Process UK Biobank data release into a DuckDB database */

params.enc = "ukb12345.enc"
params.key = "k1234r12345.key" 

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
    
    // convert to csv
    //CSV_CH = CONVERT(UKB_ENC_CH)
    
    // copy to duckdb
    //DUCK_CH = DUCK(CSV_CH)
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
    
    input:
    path(enc_ukb)
    
    output:
    tuple path('ukb.html'), path('fields.ukb')
    
    script:
    """
    dconvert $enc_ukb docs
    """
}

/* Convert the data file to CSV */
process CONVERT {
    tag "Converting data"
    
    cpus = 1
    memory = 2.GB
    time = '6h'
    
    input:
    path(enc_ukb)
    
    output:
    path('ukb.csv')
    
    script:
    """
    dconvert $enc_ukb csv
    """
}

/* Copy the CSV to DuckDB */
process DUCK {
    tag "Creating DuckDB"
    
    cpus = 8
    memory = 64.GB // check job id 26361907
    time = '6h'
    
    input:
    path(csv)
    
    output:
    path('ukb.duckdb')
    
    script:
    """
    duckdb ukb.duckdb "CREATE TABLE ukb AS SELECT * FROM read_csv_auto('$csv');"
    """
}