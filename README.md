# Nextflow workflows for processing UK Biobank releases

Process data releases from [UK Biobank](https://www.ukbiobank.ac.uk/).

## Data

Download your `.enc` encoded phenotype file and `.key` keyfile from the [AMS Portal](https://bbams.ndph.ox.ac.uk/ams/).

# d*U*c*K*d*B*: UKB Release → DuckDB Nextflow Workflow

Process new [UKB releases](https://biobank.ndph.ox.ac.uk/showcase/index.cgi) into tables in [DuckDB](https://duckdb.org/), organised by category with a built-in data dictionary.

## Initial setup

Download UKB programs
```
mkdir bin
wget -nd -P bin biobank.ndph.ox.ac.uk/ukb/util/ukbconv
wget -nd -P bin biobank.ndph.ox.ac.uk/ukb/util/dconvert
wget -nd -P bin biobank.ndph.ox.ac.uk/ukb/util/ukbunpack
for bin in bin/*; do
    chmod 755 $bin
done
```

Install required R libraries
```
Rscript -e "install.packages(c('dplyr', 'readr', 'tidyr', 'stringr', 'lubridate', 'snakecase', 'remotes'))"
Rscript -e "remotes::install_version('duckdb', '0.7.1-1')"
```

## Running the workflow

### Configuration

`process.beforeScript` should load the version of R where the DuckDB library has been installed. See [nf-core/configs](https://github.com/nf-core/configs) for possible information on your system. 

### Execution

Run the workflow on a UKB release download and key file.

```sh
nextflow run duck.nf \
-c custom.config \
-resume \
--enc ukb12345.enc \
--key k1234r12345.key
```



## Access the database

The database is written out to the `release` directory.

### R

```r
install.packages('duckdb')
```

Work with the database tables using [dbplyr](https://dbplyr.tidyverse.org)

```r
library(DBI)
library(dplyr)

con = dbConnect(duckdb::duckdb(), dbdir="release/ukb12345.duckdb", read_only=TRUE)

# read in the the data dictionary
dictionary <- tbl(con, 'Dictionary')

# get a list of tables
dbListTables(con)
##  [1] "AlgorithmicallyDefinedOutcomes" "BaselineCharacteristics"       
##  [3] "BiologicalSampling"             "BloodAssays"                   
##  [5] "CancerRegister"                 "CognitiveFunction"             
##  [7] "CognitiveFunctionOnline"        "DeathRegister"                 
##  [9] "Dictionary"                     "DietByHourRecall"              
## [11] "DigestiveHealth"                "ExperienceOfPain"              
## [13] "FirstOccurrences"               "Genotypes"                     
## [15] "HospitalInpatient"              "Imaging"                       
## [17] "LocalEnvironment"               "MentalHealth"                  
## [19] "OngoingCharacteristics"         "PhysicalActivityMeasurement"   
## [21] "PhysicalMeasures"               "PrimaryCare"                   
## [23] "ProceduralMetrics"              "Recruitment"                   
## [25] "SampleInventory"                "Touchscreen"                   
## [27] "UrineAssays"                    "VerbalInterview"               
## [29] "WorkEnvironment"     

# find which table a field is
dictionary |>
    filter(FieldID == 4598) |>
    select(Table, FieldID, Field)
## # Source:   lazy query [?? x 3]
## # Database: DuckDB 0.7.1 [R 4.1.2/release/ukb12345.duckdb]
##   Table       FieldID Field                          
##   <chr>         <dbl> <chr>                          
## 1 Touchscreen    4598 Ever depressed for a whole week

# get the Touchscreen table
touchscreen <- tbl(con, 'Touchscreen')

# work with the table using normal dplyr functions
touchscreen |>
    count(f.4598.0.0)
## # Source:   lazy query [?? x 2]
## # Database: DuckDB 0.7.1 [R 4.1.2/release/ukb12345.duckdb]
##   f.4598.0.0                n
##   <fct>                 <dbl>
## 1 NA                   329735
## 2 Yes                   89351
## 3 No                    78777
## 4 Do not know            3876
## 5 Prefer not to answer    650

# close the database
dbDisconnect(con, shutdown=TRUE)
```

### Python

```sh
pip install duckdb
```

Use the [Python API](https://github.com/duckdb/duckdb/blob/master/examples/python/duckdb-python.py) [[Reference](https://duckdb.org/docs/api/python/reference/)]
```python
import duckdb
import pandas as pd

# connect to the database
con = duckdb.connect(database='release/ukb12345.duckdb', read_only=True)

# get the Dictionary table
dictionary = con.table('Dictionary')

# find the table a particular field is in
dictionary.filter('FieldID == 4598').to_df()[['Table', 'FieldID', 'Field']]
##          Table  FieldID                            Field
## 0  Touchscreen   4598.0  Ever depressed for a whole week

# get the touchscreen table
touchscreen = con.table('Touchscreen')

# summarise the data
touchscreen.aggregate('"f.4598.0.0", count("f.4598.0.0")').to_df()
##              f.4598.0.0  count("f.4598.0.0")
## 0                   NaN                    0
## 1                    No                78777
## 2                   Yes                89351
## 3           Do not know                 3876
## 4  Prefer not to answer                  650

# close the database
con.close()
```
