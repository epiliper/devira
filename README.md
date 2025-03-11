<p align="center">
<img src="./img/adar_logo_dark.png#gh-dark-mode-only" width=80%>
<img src="./img/adar_logo_light.png#gh-light-mode-only" width=80%>
<img src="./img/adar_workflow.png" width=80%>
</p>


## About
ADAR is a pipeline for reference-guided denovo assembly of viral genomes from short-read sequencing data, tested with both shotgun and amplicon sequencing approaches. It is inspired by the Broad Institute's `assemble_denovo` workflow.

Right now, ADAR is intended primarily for assembly of respiratory viruses, including:

- Seasonal coronaviruses
- Parainfluenza
- SARS-CoV2
- Influenzas A and B
- Enteroviruses (including rhinoviruses)
## Workflow
1. ADAR takes in raw FASTQ files specified in a user-created samplesheet, performs read trimming and QC reporting, and assigns reads to specific taxon bins with [kraken2](https://github.com/DerrickWood/kraken2). 
2. The reads associated with each bin are then assembled into contigs with [megahit](https://github.com/voutcn/megahit), and the contigs are compared to FASTA references in a user-supplied database; for all references above average nucleotide identity and coverage thresholds, ADAR will use the references to guide and refine contig arrangement into scaffolds. 
3. Scaffold gaps undergo gap-filling with 1) sequencing reads and, if any gaps are left, 2) reference sequence. 
4. To generate a consensus genome, sequencing reads are realigned to the scaffolds with [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2), and the alignment is used to call consensus with [ivar consensus](https://github.com/andersen-lab/ivar). This is performed twice to get longer assemblies with less bias towards the chosen reference genome.


## Instructions
0. Ensure the Docker desktop client is running. [Install it](https://docs.docker.com/get-started/get-docker/) if you haven't already.
1. Arrange all input fastqs in their own directory. It's recommended to have single-end and paired-end libraries in separate directories.
2. From this repository, download `bin/fastq_dir_to_samplesheet.py` script. Run the script as follows:   

    ```bash 
    python3 fastq_dir_to_samplesheet.py $PATH_TO_FASTQ_DIR \\
        -r1 $R1_SUFFIX -r2 $R2_SUFFIX \\
        ${SAMPLESHEET_NAME}.tsv
    ```
    where `$R1_SUFFIX` and `$R2_SUFFIX` are the unique suffixes between the paired-end files of each sample, e.g. 
    ```bash
    -r1 _1.fastq.gz -r2 _2.fastq.gz
    ```
    If creating a samplesheet for single-end libraries, just use `-r1`.

3. Run the pipeline with the samplesheet as input.
    ```bash
    nextflow run epiliper/adar -r main -latest \\
        --input ${SAMPLESHEET_NAME}.tsv --output $OUTPUT_DIR \\
        -profile docker
    ```
4. Once the pipeline is done, consensus genomes can be found in `$OUTPUT_DIR/final_files/final_assemblies`. For reports of any genomes/samples that failed assembly QC thresholds, see files in `$OUTPUT_DIR/fail`.

You can test this with out with the pipeline's example data, located in `fastq_example`.

## Notes
After trimmed reads are classified with Kraken2 early in the pipeline, reads associated with a specific taxon are extracted into a new file, and this is repeated for all taxon IDs listed in a TSV file input to the pipeline. The default taxon ID list can be found in `assets/taxids.tsv`. Reads not related to any of the taxon IDs in the list, or unclassified reads period, are not used in downstream assembly. 

> [!IMPORTANT]   
> If you want to make your own ID list, it's recommended to keep taxon IDs at species-level or higher, to avoid over-stratification of reads and lots of unnecessary file/generation processing.

---

The reference database used for selecting sequences to guide scaffolding is the same as in our reference-based assembly pipeline, [revica-strm](https://github.com/greninger-lab/revica-strm). It's comprised of multiple representatives of a variety of respiratory virus species such as enterovirus, seasonal coronavirus, SARS-CoV2, parainfluenza, measles, influenza, and more. Inspect `assets/ref.fa` if curious. If you intend to use your own database, ensure the fasta headers are structured as follows:   

```ACCESSION<SPACE>REF_TAG<SPACE>SAMPLE_HEADER```   

where REF_TAG should be unique to a species-specific segment/genome. Take these entries for Flu A segments PB1 and NS1, and an enterovirus genome, for example:

```
>NC_007364.1 fluA_NS1 Influenza A virus (A/goose/Guangdong/1/1996(H5N1)) segment 8, complete sequence
>NC_007375.1 fluA_PB1 Influenza A virus (A/Korea/426/1968(H2N2)) segment 2, complete sequence
>AF406813.1 EV Porcine enterovirus 8 strain V13, complete genome
```
You can specify your own reference database with `--refs $REF_DB`.

The kraken2 database we use is generated solely from the fasta sequences in `assets/ref.fa`. To add your own sequences, you'll have to recreate the database. See [our instructions for how to do so](creating_kraken2_db.md). You can specify a path to your Kraken2 database with `--kraken2_db $DB`.
