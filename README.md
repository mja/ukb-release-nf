# d*U*c*K*d*B*: UKB Release → DuckDB Nextflow Workflow

Process new [UKB releases](https://biobank.ndph.ox.ac.uk/showcase/index.cgi) into tables in [DuckDB](https://duckdb.org/), organised by category with a built-in data dictionary.

## Initial setup

Download UKB and database programs
```
mkdir bin
wget -nd -P bin biobank.ndph.ox.ac.uk/ukb/util/ukbconv
wget -nd -P bin biobank.ndph.ox.ac.uk/ukb/util/dconvert
wget -nd -P bin biobank.ndph.ox.ac.uk/ukb/util/ukbunpack
for bin in bin/*; do
    chmod 755 $bin
done
cp /exports/igmm/eddie/GenScotDepression/local/bin/duckdb-0.7.1 bin/duckdb
```

Install required R libraries
```
module load igmm/apps/R/4.1.0
Rscript -e "install.packages(c('dplyr', 'readr', 'tidyr', 'stringr', 'lubridate', 'snakecase', 'remotes'))"
Rscript -e "remotes::install_version('duckdb', '0.7.1-1')"
```

## Running the workflow

The workflow is run from an interactive session. Use a screen manager so that you can disconnect from the session while the workflow runs. 

Connect to Eddie and start a [tmux](https://www.redhat.com/sysadmin/introduction-tmux-linux) session
```sh
ssh UUN@eddie.ecdf.ed.ac.uk
tmux new-session -s ukb
```
- start an interactive session
```sh
qlogin -l h_vmem=8G
```
Navigate back to the `ukb-release-nf` directory

Load the UGE and Nextflow modules
```sh
module load uge
module load igmm/apps/nextflow/20.04.1
```

Run the workflow on a UKB release download and key file:
```sh
nextflow run ukb.nf \
-c eddie.config \
-work-dir /exports/eddie/scratch/$USER/ukb-release \
-resume \
--enc /exports/igmm/datastore/GenScotDepression/data/ukb/release/ukb670429/ukb670429.enc \
--key /exports/igmm/datastore/GenScotDepression/data/ukb/release/k4844-keys/k4844r670429.key
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

con = dbConnect(duckdb::duckdb(), dbdir="release/ukb670429.duckdb", read_only=TRUE)

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
## # Database: DuckDB 0.6.1
## #   4.1.0/release/ukb670429.duckdb]
##   Table       FieldID Field                          
##   <chr>         <dbl> <chr>                          
## 1 Touchscreen    4598 Ever depressed for a whole week

# get the Touchscreen table
touchscreen <- tbl(con, 'Touchscreen')

# work with the table using normal dplyr functions
touchscreen |>
    count(f.4598.0.0)
## # Source:   lazy query [?? x 2]
## # Database: DuckDB 0.6.1 [madams23@Linux 3.10.0-1160.71.1.el7.x86_64:R
## #   4.1.0/release/ukb670429.duckdb]
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
con = duckdb.connect(database='release/ukb670429.duckdb', read_only=True)

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
