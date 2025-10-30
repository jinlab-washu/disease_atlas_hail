# Hail annotation pipeline
Hail pipeline used to prepare VCF file to export into elastic search database used for the WashU Disease Atlas https://washudiseaseatlas.org

## Requirements
```
Docker image: dreammaerd/hail_vep_gnomad:0.2.120_104_0.6.4_fixJson

Docker image contains the following key software
Python 3.8.10
Hail 0.2.120
Java 1.8
ensembl-vep-release-104.3
```

## High level summary of Pipeline steps
```
Step 1: Break multi-allelics to bi-allelics and use minimal representation of alleles
Step 2: Count of alleles by population
Step 3: Annotate using Variant Effect Predictor (VEP) with LOFTEE plugin
Step 4: Reformat for export to Hail Table and VCF file
Step 5: Reformat and export to Elastic Search
```

## Inputs and Outputs
**Input 1:** bgzipped VCF created by GATK. Other VCF files created by other variant callers has not been tested  
**Input 2:** Tab separated meta file. See example_meta.tsv for format.

**Output:** Hail table - individual genotypes removed, VEP functional annotation, various summary counts/metrics, flattened structure for export to Elastic Search.

## Example: Annotation and output to hail table
```
# Start up docker
LSF_DOCKER_PORTS="$JPORT:$JPORT" bsub -Is -G compute-jin810-t3 -q subscription \
-n 16 -R 'select[port8888=1] span[hosts=1] rusage[mem=128GB]' \
-a 'docker(dreammaerd/hail_vep_gnomad:0.2.120_104_0.6.4_fixJson)' /bin/bash

#Run annotation command
python3 hail_annotate_pipeline.py -i input.vcf.gz \
-m meta_file.tsv \
-o output.ht

```

## Example: Upload hail table to Elastic Search database
```
python populate_data.py -i input.ht
```
