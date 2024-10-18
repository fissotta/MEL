# Installation

```bash
# Create a new environment and install metaSNV
mamba create --name metaSNV -c bioconda -c conda-forge 'metasnv>=2.0.1'
```

[metaSNV GitHub Repository](https://github.com/metasnv-tool/metaSNV/tree/master)

```bash
# Activate the metaSNV environment
mamba activate metaSNV
```

## FIX 1

Download `htslib-1.11` from [here](https://sourceforge.net/projects/samtools/files/):

```bash
# Navigate to the htslib-1.11 directory and build
cd htslib-1.11 && make
```

Download `samtools-1.11` from [here](https://sourceforge.net/projects/samtools/files/). Note: **Do not use the latest version** because `qaTools` requires a library called `sam.h` that has been renamed in later versions.

```bash
# Navigate to the samtools-1.11 directory and build
cd samtools-1.11 && make

# Clone qaTools repository
git clone https://github.com/CosteaPaul/qaTools.git
```

Edit `qaTools` Makefile as described in the [qaTools GitHub](https://github.com/CosteaPaul/qaTools) so that `HTSLIB` and `SAMTOOLS` point to the `htslib-1.11` and `samtools-1.11` paths.

```bash
# Export library path
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/path/to/htslib-1.11"

# Navigate to qaTools directory and build
cd qaTools && make
```

### Modify `metaSNV.py`:
Update:
```python
cmd = ['{}/src/qaTools/qaCompute'.format(basedir),
```
to:
```python
cmd = ['/path/to/qaTools/qaCompute',
```

Located in: `/path/to/anaconda3/share/metasnv-2.0.4/metaSNV.py`

## FIX 2

### Move Line in `metaSNV_subpopr.R`
Find the script `metaSNV_subpopr.R` and move line 365 above line 358.

**Before:**
```r
bpParam <- MulticoreParam(workers = min(N.CORES,length(species)),
                          jobname = "subpopr",
                          stop.on.error = FALSE,
                          threshold = "DEBUG",
                          log = TRUE,
                          progressbar = printProgressBar,
                          logdir = paste0(OUT.DIR,"/threadLogs"))
dir.create(paste0(OUT.DIR,"/threadLogs"), recursive = T, showWarnings = FALSE)
```

**After:**
```r
dir.create(paste0(OUT.DIR,"/threadLogs"), recursive = T, showWarnings = FALSE)
bpParam <- MulticoreParam(workers = min(N.CORES,length(species)),
                          jobname = "subpopr",
                          stop.on.error = FALSE,
                          threshold = "DEBUG",
                          log = TRUE,
                          progressbar = printProgressBar,
                          logdir = paste0(OUT.DIR,"/threadLogs"))
```

## Input Preparation

### Join with 100 Ns

```bash
# Joining contigs
for f in *fna; do perl ../contigs_joiner.pl $f > ${f%.fna}_J.fna; done

# Replace ".fasta" or ".fna" extensions
sed -i 's/\.fasta//g' *_J.fasta
sed -i 's/\.fna//g' *_J.fna
```

### CoverM Genome Coverage

```bash
coverm genome --methods length covered_bases count mean relative_abundance rpkm --min-read-aligned-length 45 --min-read-percent-identity 97 --min-covered-fraction 0 --discard-unmapped --output-format sparse --genome-fasta-files references/NC_009767.fna --genome-fasta-extension fna --bam-file-cache-directory NC_009767_MGTA --threads 90 --output-file NC_009767_MGTA.txt -1 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto10_GGTCACGA-GTATTATG_L001_R1_001_TRIM.fq -2 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto10_GGTCACGA-GTATTATG_L001_R2_001_TRIM.fq
```

### CoverM Contigs

```bash
coverm contig --methods length covered_bases covered_fraction reads_per_base count mean trimmed_mean rpkm tpm variance --include-secondary --min-read-aligned-length 45 --min-read-percent-identity 97 --min-covered-fraction 0 --discard-unmapped --output-format sparse --sharded --reference references --bam-file-cache-directory bams_coverm --threads 60 --output-file AFP_ICO_FCA_coverm.txt -1 /media/WALLROSE/FASTQ_MG/LSM1_2019_1.fq -2 /media/WALLROSE/FASTQ_MG/LSM1_2019_2.fq ...
```

### Check Headers

```bash
# Extract and modify SAM header
samtools view -H coverm-genome.S16_2020_1.fq.gz.bam > header.sam && nano header.sam
```

If the header does not match the reference:
```bash
# Replace SAM header in all BAM files
parallel --jobs 20 'samtools reheader header.sam {} > {}.rehead' ::: *bam
rm *bam && rename 's/.rehead//' *rehead && rm header.sam
```

### Sort and Index BAM Files

```bash
# Sort BAM files
for f in *bam; do samtools sort -@ 30 -o ${f%.bam}.sorted.bam $f; done

# Index BAM files
for f in *sorted.bam; do samtools index -@ 30 $f; done
```

### Prodigal GFF Conversion

```bash
# Generate GFF and protein files using Prodigal
for file in *.fasta; do prodigal -i "$file" -o "${file%.fasta}.gff" -f gff -a "${file%.fasta}_proteins.faa"; done

# Clean up GFF files
sed -i '/Model/d' *gff
sed -i '/Sequence/d' *gff
```

### Convert GTF to GFF3

```bash
# Convert GTF to GFF3
parallel --jobs 4 'agat_convert_sp_gxf2gxf.pl -g {} -o {}.gff3' ::: *gtf
```

## Execution of metaSNV

```bash
# Activate environment
mamba activate metaSNV

# Index the reference
samtools faidx NC_010175.fna
```

### Running metaSNV

```bash
# Example run
metaSNV.py --threads 96 NC_010175_metaSNV NC_010175_bam.list references/NC_010175.fna --db_ann metaSNV_anntotations/NC_010175_metaSNV_anntotations.txt

# MetaSNV Filtering
metaSNV_Filtering.py --n_threads 30 NC_010175_metaSNV
```

### Fix LB

```bash
# Filter using AWK script
awk '{ delete_line=0; for (i=2; i<=NF; i++) if (($i+0) > 0 && ($i+0) < 0.05) { delete_line=1; break } } !delete_line'  S3_2020_Vamb_5.filtered.freq > S3_2020_Vamb_5.filtered_LB.freq

# Run DistDiv
metaSNV_DistDiv.py --n_threads 30 --dist --div --divNS --matched --filt NC_010175_metaSNV/filtered/pop/

# Subpopulation Analysis
metaSNV_subpopr.R -i NC_010175_metaSNV -p 12
```

### Generate GFF

```bash
# Generate GFF files using Prodigal
for f in *fna; do prodigal -i $f -f gff -o ${f%.fna}.gff -a ${f%.fna}.faa; done
```

### MetaSNV Scripts

```bash
# Run metaSNV pipeline steps
for f in *fna; do metaSNV.py --threads 96 ${f%.fna}_metaSNV bam_${f%.fna}_list_bams.txt $f; done

for f in *_metaSNV; do metaSNV_Filtering.py --n_threads 96 $f; done

# DistDiv analysis without GFF
for f in *_metaSNV; do metaSNV_DistDiv.py --n_threads 90 --dist --div --matched --filt $f/filtered/pop/; done
```
