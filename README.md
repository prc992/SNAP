# SNAP Pipeline Documentation

## Input Data Options
There are two ways to provide input files to the pipeline:  
1. A **directory** containing the files to be processed.  
2. A **spreadsheet** with basic metadata about the files.  

The input files can be of three types:  
- **FASTA files**  
- **Raw BAM files** (generated immediately after alignment with a reference genome)  
- **Processed BAM files** (expected to be sorted, deduplicated, and filtered for unique reads)  

Below is a detailed explanation of how to use each method.  

---

## Providing Input as FASTA Files
If your sample files are in **FASTA format**, use the `--sample_dir_fasta` parameter. The pipeline assumes that each sample is stored in a separate directory, and all files within that directory belong to the corresponding sample.  

### Example Directory Structure:

```
Sample_folder/
├── Sample1/
│   ├── Sample1_ABC_123_1.fasta
│   ├── Sample1_ABC_123_2.fasta
├── Sample2/
│   ├── Sample2_ABC_123_1.fasta
│   ├── Sample2_ABC_123_2.fasta
├── Sample3/
│   ├── Sample3_ABC_123.fasta
```

### Automatic Spreadsheet Generation
When using this method, the pipeline will generate a **spreadsheet** (`CSV format`) with the following structure: 

sampleId, enrichment_mark, read1, read2

- The pipeline supports both **paired-end** and **single-end** FASTA files.
- You can specify which **enrichment mark** should be calculated using the `--enrichment_mark` parameter.

### Enrichment Marks
By default, the SNAP pipeline supports the following enrichment marks:
- **H3K4me3**
- **H3K27ac**
- **MeDIP**

The reference files required for this calculation are stored in:

ref_files/enrichment_states/enrichment_mark/

If you need to provide **custom enrichment mark files**, place the **on-target** and **off-target** BED files in a local directory and name it after the desired histone mark. The expected file names are:
- `off.target.filt.bed`
- `on.target.filt.bed`

Then, set:
- `--enrichment_mark <custom_mark_name>`
- `--enrichment_states_ref <path_to_directory>`

📌 **Note:** If using the **directory-based** input method, all samples will have the **same enrichment mark** applied. If you need to assign different histone marks to individual samples within the same execution, you must provide a **spreadsheet** instead.

---

## Providing Input via Spreadsheet (FASTA)
To provide input using a spreadsheet, use the `--samplesheetfasta` parameter. This expects a **CSV file** with the following structure:

sampleId, enrichment_mark, read1, read2

- The **enrichment_mark** field can be left blank if no enrichment mark calculation is required.
- The **read2** field can be left blank for **single-end** FASTA files.

---

## Providing Input as BAM Files
If your data is already aligned (BAM format), you can use one of the following options:
- **Directory-based input:** `--sample_dir_bam`
- **Spreadsheet-based input:** `--samplesheetBams`

### Example Directory Structure for BAM Files:

```
Sample_folder/
│── Sample1/
│   ├── Sample1_ABC_123.bam
│── Sample2/
│   ├── Sample2_ABC_123.bam
│── Sample3/
│   ├── Sample3_ABC_123.bam
```

### Providing Input via Spreadsheet (BAM)
To provide BAM files via a **spreadsheet**, use the `--samplesheetBams` parameter. The expected CSV format is:

sampleId, enrichment_mark, bam

- The **enrichment_mark** field is optional.

📌 **Automatic Spreadsheet Generation:**
Whenever a directory is provided using `--sample_dir_fasta` or `--sample_dir_bam`, the pipeline will generate a **spreadsheet** in the **output directory**. The file will be named:
- `snap-samplesheet-bam-date_and_time.csv`
- `snap-samplesheet-fasta-date_and_time.csv`

---

## Handling Pre-Processed BAM Files
If your BAM files are **already processed** (sorted, deduplicated, and filtered for unique reads), you can **skip redundant processing steps** by setting:

–deduped_bam true

This should be used in combination with either `--samplesheetBams` or `--sample_dir_bam`.

📌 **Skipping Steps:** When `--deduped_bam true` is set, the pipeline will **bypass**:
- BAM sorting
- Library complexity calculation
- Unique read filtering
- Quality filtering
- `samtools` stats calculation
- Duplicate removal

This reduces processing time if your BAM files have already undergone these steps.

---

## Reference Genome Selection
Specify the reference genome version using: 

–genome <reference_genome>

The SNAP pipeline is pre-configured to support:
- **hg19**
- **hg38**

For these references, no additional downloads are required—the necessary FASTA files, **blacklist regions**, and SNP files for sample identification are already available in:
ref_files/genome/genome_paths.csv

If using a **custom reference genome**, provide a **local version** of this file and specify:
- `--genomeInfoPaths <full_path_to_genome_paths.csv>`
- `--genome <custom_genome_name>` (must match the **Genome** field in the CSV file)

---

## Pipeline Execution Phases
The output directory where results will be saved is specified with:
–outputFolder <directory_path>

(Default: `analysis_folder`)

The SNAP pipeline consists of **six main processing stages**. By default, all stages are executed, but you can stop the pipeline at a specific stage using the `--until` parameter:

```
PREPROCESSING
DOWNLOAD_REFERENCES
ALIGNMENT
BAM_PROCESSING
FRAGMENTS_PROCESSING
BAM_SIGNAL_PROCESSING
```
---

## Excluding Blacklisted Genomic Regions
By default, the SNAP pipeline **removes genomic regions known to cause artifacts**, reducing false positives and biases in the analysis.

To **disable** this step, set:

–exclude_dac_regions false

---

## END Motif and GC Content Analysis
The pipeline generates an **END motif analysis file** at:motifs/bp_motif.bed

By default, this analysis considers **4-mers**, but this value can be adjusted using: motifs/bp_motif.bed

---

## Pileup Reports
The SNAP pipeline generates **pileup reports** in two formats:
1. **MultiQC report** (found in `reports/multiqc/`)
2. **IGV session file** (stored in `reports/igv_session/`, viewable in IGV Browser)

By default, both reports use regions defined in: ref_files/pileup_report/test_housekeeping.bed

To use **custom regions**, specify:–genes_pileup_report <path_to_custom_bed_file>
