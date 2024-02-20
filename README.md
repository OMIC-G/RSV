# RSV genomic surveillance in Galicia

The main goal of genomic surveillance of RSV is to report the geographical spread and temporal patterns of RSV clades and to detect resistance mutations. The objective is to contribute to reducing the burden of the seasonal RSV.

Surveillance of respiratory viruses has a multidisciplinary character, and must integrate clinical, epidemiological and virological data in an organized way to report to the ECDC by TESSY (The European Surveillance System).

As a part of the RELECOV (Red de Laboratorios de secuenciación de SARS-CoV-2) and for coordination with CNM-ISCIII (Centro Nacional de Microbiología, Instituto de Salud Carlos III), several activities are part of the routine of a sequencing laboratory.

Here we present an integrated RSV analysis pipeline including :
1. Sequencing and bioinformatic analysis
2. Submission to GISAID of consensus sequences
3. Report of clades to autonomous community Health Authorities (SERGAS)

This is the result of a continuous collaborative effort of the OMIC-G network (Red de Laboratorios para la aplicación de Ómicas a la Microbiología Clínica en Galicia)

## Dependencies

The following software was used for this pipeline.

* BBSplit from [BBMap](https://sourceforge.net/projects/bbmap/)
* [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)
* [samtools](https://www.htslib.org/)
* [ivar](https://andersen-lab.github.io/ivar/html/manualpage.html)
* [nextclade](https://github.com/nextstrain/nextclade)
* [multiqc](https://multiqc.info/)
* [rsvCLI](https://gisaid.org/)


## Sample files and folder structure

The starting point of this analysis are .fastq files stored in the folder *fastq*

The following folders are used through the whole analysis:

* *fastq*: with the original .fastq.gz files
* *gisaid*: for uploading the final .fasta and .csv files to GISAID
* *mid*: temporary files created in the process and .fa files of each sample
* *out*: final files created in the analysis: report, fasta, csv, ...
* *refs*: .fasta reference files for each subtype (A and B) and primer bed files

First time, before alignment, you have to index the fasta references with bwa:

```
bwa-mem2 index -p RSVA RSVA.fasta
bwa-mem2 index -p RSVB RSVB.fasta
```


## 1. RSV typing: RSVa, RSVb

We use the `bbsplit` utility from _BBmap_ in order to align the fastq sequences against each one of the 2 references (RSVa and RSVb)

```
FILE="sample_name"
bbsplit.sh build=1 \
    in=fastq/${FILE}_L001_R1_001.fastq.gz \
    in2=fastq/${FILE}_L001_R2_001.fastq.gz \
    ref_a=refs/RSVA.fasta out_a=fastq/pre/A${FILE}.fastq.gz \
    ref_b=refs/RSVB.fasta out_b=fastq/pre/B${FILE}.fastq.gz
```

As a result, we get one .fastq.gz file for each sample and each type, stored in the folder _fastq/pre_ and named like this:
`${TYPE}${SAMPLE}.fastq.gz`

After that, we align the sequences in each obtained fast.gz file with `bwa-mem2` and study the coverage of each one with `samtools`. We only keep the ones with a coverage >80%. Usually one for each sample, except in case of coinfections.
 

```
FILE=${TYPE}${SAMPLE}

# Aligning:
bwa-mem2 mem -t 8 RSV${TYPE} fastq/pre/${FILE}.fastq.gz > mid/${FILE}.sam

# Converting from .sam to .bam:
samtools view -bS mid/${FILE}.sam -o mid/${FILE}.bam

# Sorting:
samtools sort mid/${FILE}.bam -o mid/${FILE}.sorted.bam

# Triming primers:
ivar trim -i mid/${FILE}.sorted.bam -b refs/RSV${TYPE}.primer.bed -p mid/${FILE}.trimmed -e -m 32

# Re-sorting:
samtools sort mid/${FILE}.trimmed.bam -o mid/${FILE}.sorted.bam

# Choosing the best between RSV A/B
samtools index mid/${FILE}.sorted.bam
if [[ $FILE == A* ]];then
    cover=$(samtools coverage -H -r MN078114.1:4900-5800 mid/${FILE}.sorted.bam|cut -f6)
else
    cover=$(samtools coverage -H -r ON729320.1:4900-5800 mid/${FILE}.sorted.bam|cut -f6)
fi

# Cleaning:
if (( $(echo "$cover < 80" |bc -l) )); then
    rm mid/${FILE}*
else
    rm mid/${FILE}.sam mid/${FILE}.bam mid/${FILE}.trimmed.bam
fi
```

As a result, in the _mid_ folder we get all the .bam files of the samples with this name structure:
`${TYPE}${SAMPLE}.sorted.bam`


## 2. Generation of .fasta consensus of each sample

We use `samtools` and `ivar` in order to get the .fasta file of the consensus of each sample, using the previously obtained .sorted.bam files

```
FILE=${TYPE}${SAMPLE}
samtools mpileup -aa -A -d 0 -Q 0 mid/${FILE}.sorted.bam \
    | ivar consensus -t 0 -p mid/${FILE} -i ${FILE}
```

As a result, in the _mid_ folder we get the .fa files with the same name structure:
`${TYPE}${SAMPLE}.fa`

Besides, we create a single .fasta file for each subtype of all samples. We'll need these files later to query nextclade. 
```
cat mid/A*.fa > out/a.fasta
cat mid/B*.fa > out/b.fasta
```

## 3. Sumarize quality and coverage data

We obtain coverage data at 10X and 100X for each sample and store the results in a .csv file called _coverage.csv_. Later, we will use this file as part of the final report.

```
FILE=${TYPE}${SAMPLE}
samtools coverage -H mid/${FILE}.sorted.bam
samtools mpileup mid/${FILE}.sorted.bam
```


## 4. Geting information of each sample from nextclade

We use the `nextclade CLI` from Nextstrain to obtain data for each sample. For this, we need to provide the .fasta file of all the samples of each type, so we query Nextstrain for each different reference

```
nextclade run --input-dataset=rsv_a --output-csv=out/nextclade_a.csv out/a.fasta
nextclade run --input-dataset=rsv_b --output-csv=out/nextclade_b.csv out/b.fasta
```

The results, stored as .csv files in the _out_ folder, also will be used as a part of the final report.


## 5. Variant calling

We get data for minority variants, so we can make histogram plots in the final report based on this data.

```
FILE=${TYPE}${SAMPLE}
samtools mpileup -A -d 0 --reference refs/RSV${TYPE}.fasta \
    -Q 0 mid/${FILE}.sorted.bam | ivar variants -p out/tsv/${FILE} \
    -t 0.03 -r refs/RSV${TYPE}.fasta -m 10 
```


## 6. Final report

The final report is made with an R script, summarizing all data obtained previously plus quality data obtained with `FastQC`, `Qualimap` and `multiqc` tools.

```
samtools stats ${FILE} > mid/${FILE}.stats

fastqc fastq/A*.fastq.gz --outdir=mid/fastqc/
fastqc fastq/B*.fastq.gz --outdir=mid/fastqc/

qualimap multi-bamqc -d mid/qualimap_list_a.txt \
    -gff genemap_a.gff -outdir out/qualimap_A/ -r mid/A*.sorted.bam
qualimap multi-bamqc -d mid/qualimap_list_b.txt \
    -gff genemap_b.gff -outdir out/qualimap_B/ -r mid/B*.sorted.bam

multiqc --force -o out -n "multiqc.html" -i "Report" \
    -b "<a href='report.html'>Global Report</a>" mid
```


## RSV consensus sequences submission to GISAID

The objective of this part is to ease the process to send the metadata and .fasta files to *GISAID* database. It's done with an R script and with the functionality of the `rsvCLI` application from GISAID. 

In order to use the `rsvCLI` command line utility, a GISAID user and a client_id it's needed.


### Metadata and fasta files

You need to complete the `gisaid_rsv_template.csv` with the metadata of the samples (date, patient age, gender, Lab ID, ...). In our case, we use an .ods file from the LIS and process it with R to fill the template.csv.

The .fasta file with the sequence of each sample is done by concatenating the individual .fa files in the _mid_ folder. The names of each sequence in the final .fasta file, should match the names in the template.csv. Again, we use R to do the work.


### Uploading to GISAID

The files generated in the previous step (.fasta and .csv) are uploaded to GISAID with the `rsvCLI` utility. 

```
rsvCLI upload --username XXXX --password YYYY --clientid ZZZZ \
    --metadata metadata.csv --fasta sequences.fasta \
    --dateformat YYYYMMDD --log result.log
```

In the _result.log_ generated you'll have the accession_id of each sample uploaded to GISAID.


## Final report to SAÚDE PÚBLICA DE GALICIA

In this final step, an .xlsx file is created with the data requested by Saúde Pública de Galicia. 
Again, we use R to write the .xlsx file including:
- accession_id numbers from the rsvCLI log
- patient data and clade from the .ods from the LIS
- GISAID sample names from the .csv sent to GISAID


