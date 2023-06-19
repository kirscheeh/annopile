# Read PileUp with Custom Variant Annotation

Simplistic pileups of specified regions with basic CNV calling and coverage plotting [soon].

## Set-Up
## Prerequisities
- python and the following packages: pysam, pandas, openpyxl
- docker
- current version of annovar databases for hg19 

### Config Files
To run this small gathering of scripts, first customize the [config file](path.conf) containing the paths to several locations.
- DATA_PATH: parent directory containing subfolders of data (not mandatory)
- PARENT_ANNOVAR: path to folder containing hg19 annovar database (mandatory) --> inside, a folder hg19_annovar_db is needed containing the database files
- REPO_PATH: path to this reporitory (mandatory)
- RESULT_PATH: path to output directory (mandatory)

If you want to automatically generate the SIO.csv file, which is the main file containing information about the data you want to process, you can use [the file gathering script](scripts/find_files.py)- This will be only usable if you adapt it to your local infrastructure. It is important, however, tthat the [samples of interest file](input/SOI.csv) is a comma-separated list containing path/to/bam/file.bam and sample_name. 

The [bed file](input/GOI.bed) is a simple 4-column bed file containing the regions you are interest in.

### Docker
Run the following command to load the annovar docker image:

    docker load -i /path/to/repo/annovar_docker.tar

## Usage
Simply call

    bash scripts/createSummary.sh

and find the results in your deisgnated results folder.

### Variant Calls
Since the pileup is far more sensitive than default variant callers, only SNVs are reported in the VCF sheet that occure with an allele frequency of 0.1. Due to the difficulty of handling MNVs and InDels, this script does not include them in the variant call and annotation.

## TODO
- add visualizations for coverage
- add annovar container
