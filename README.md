# UKB Release Nextflow Workflow

Process new UKB releases.

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
cp /exports/igmm/eddie/GenScotDepression/local/bin/duckdb bin/duckdb
```

Install required R libraries
```
module load igmm/apps/R/4.1.0
Rscript -e "install.packages(c('dplyr', 'readr', 'tidyr', 'stringr', 'snakecase', 'duckdb'))"
```

## Running the workflow

The workflow is run from an interactive session. Use a screen manager so that you can disconnect from the session while the workflow runs. 

- Connect to Eddie and start a [tmux](https://www.redhat.com/sysadmin/introduction-tmux-linux) session
```sh
ssh UUN@eddie.ecdf.ed.ac.uk
tmux new-session -s ukb
```
- start an interactive session
```sh
qlogin -l h_vmem=8G
```
- navigate back to the `ukb-release-nf` directory

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

