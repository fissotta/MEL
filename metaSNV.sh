####################
### Installation ###
#####################

mamba create --name metaSNV -c bioconda -c conda-forge 'metasnv>=2.0.1'

https://github.com/metasnv-tool/metaSNV/tree/master

mamba activate metaSNV

#############
### FIX 1 ###
#############

download htslib-1.11 @ https://sourceforge.net/projects/samtools/files/

cd htslib-1.11 & make

#download samtools-1.11 @ https://sourceforge.net/projects/samtools/files/ (if you download or clone the last version, it won't work because qaTools requires a library called sam.h that's renamed in later versions)

cd samtools-1.11 & make

git clone https://github.com/CosteaPaul/qaTools.git

#edit qaTools' Makefile as described in https://github.com/CosteaPaul/qaTools so that HTSLIB and SAMTOOLS point at the htslib-1.11 and samtools-1.11 paths

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/path/to/htslib-1.11"

cd qaTools & make

#change

 cmd = ['{}/src/qaTools/qaCompute'.format(basedir),
#to

    cmd = ['/path/to/qaTools/qaCompute',
in /path/to/anaconda3/share/metasnv-2.0.4/metaSNV.py



#############
### FIX 2 ###
#############

# have solved this problem.
# Find the location of the script "metaSNV_subpopr.R" and move the line 365 above the line 358:
# Before:

 ~/mambaforge/envs/metaSNV/share/metasnv-2.0.4/metaSNV_subpopr.R
 
bpParam <- MulticoreParam(workers = min(N.CORES,length(species)),
                          jobname = "subpopr",
                          stop.on.error = FALSE,
                          threshold = "DEBUG",
                          log = TRUE,
                          progressbar = printProgressBar,
                          logdir = paste0(OUT.DIR,"/threadLogs"))
dir.create(paste0(OUT.DIR,"/threadLogs"), recursive = T, showWarnings = FALSE)

# After:

dir.create(paste0(OUT.DIR,"/threadLogs"), recursive = T, showWarnings = FALSE)
bpParam <- MulticoreParam(workers = min(N.CORES,length(species)),
                          jobname = "subpopr",
                          stop.on.error = FALSE,
                          threshold = "DEBUG",
                          log = TRUE,
                          progressbar = printProgressBar,
                          logdir = paste0(OUT.DIR,"/threadLogs"))


######################### 
### Input preparation ### 
######################### 

#JOIN with 100 Ns
for f in *fna; do perl ../contigs_joiner.pl $f > ${f%.fna}\_J.fna; done

#FR with coverm
coverm genome --methods length covered_bases count mean relative_abundance rpkm --min-read-aligned-length 45 --min-read-percent-identity 97 --min-covered-fraction 0 --discard-unmapped --output-format sparse --genome-fasta-files references/NC_009767.fna --genome-fasta-extension fna --bam-file-cache-directory NC_009767_MGTA --threads 90 --output-file NC_009767_MGTA.txt -1 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto10_GGTCACGA-GTATTATG_L001_R1_001_TRIM.fq -2 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto10_GGTCACGA-GTATTATG_L001_R2_001_TRIM.fq -1 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto11_48_dna_GAATCCGTGG-TCACATGAGA_L001_R1_001_TRIM.fq -2 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto11_48_dna_GAATCCGTGG-TCACATGAGA_L001_R2_001_TRIM.fq -1 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto11_62_dna_ACGTGAACAT-ACTAGCCGTG_L001_R1_001_TRIM.fq -2 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto11_62_dna_ACGTGAACAT-ACTAGCCGTG_L001_R2_001_TRIM.fq -1 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto11_CTGCTTCC-GATAGATC_L001_R1_001_TRIM.fq -2 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto11_CTGCTTCC-GATAGATC_L001_R2_001_TRIM.fq -1 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto12_TCATCCTT-AGCGAGCT_L001_R1_001_TRIM.fq -2 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto12_TCATCCTT-AGCGAGCT_L001_R2_001_TRIM.fq -1 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto13_AGGTTATA-CAGTTCCG_L001_R1_001_TRIM.fq -2 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto13_AGGTTATA-CAGTTCCG_L001_R2_001_TRIM.fq -1 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto1_CAATTAAC-CGAGATAT_L001_R1_001_TRIM.fq -2 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto1_CAATTAAC-CGAGATAT_L001_R2_001_TRIM.fq -1 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto2_TGGCCGGT-TAGAGCGC_L001_R1_001_TRIM.fq -2 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto2_TGGCCGGT-TAGAGCGC_L001_R2_001_TRIM.fq -1 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto3_AGTACTCC-AACCTGTT_L001_R1_001_TRIM.fq -2 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto3_AGTACTCC-AACCTGTT_L001_R2_001_TRIM.fq -1 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto4_GACGTCTT-GGTTCACC_L001_R1_001_TRIM.fq -2 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto4_GACGTCTT-GGTTCACC_L001_R2_001_TRIM.fq -1 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto5_TGCGAGAC-CATTGTTG_L001_R1_001_TRIM.fq -2 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto5_TGCGAGAC-CATTGTTG_L001_R2_001_TRIM.fq -1 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto6_CATAGAGT-TGCCACCA_L001_R1_001_TRIM.fq -2 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto6_CATAGAGT-TGCCACCA_L001_R2_001_TRIM.fq -1 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto7_ACAGGCGC-CTCTGCCT_L001_R1_001_TRIM.fq -2 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto7_ACAGGCGC-CTCTGCCT_L001_R2_001_TRIM.fq -1 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto8_GTGAATAT-TCTCATTC_L001_R1_001_TRIM.fq -2 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto8_GTGAATAT-TCTCATTC_L001_R2_001_TRIM.fq -1 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto9_AACTGTAG-ACGCCGCA_L001_R1_001_TRIM.fq -2 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto9_AACTGTAG-ACGCCGCA_L001_R2_001_TRIM.fq -1 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto_3_45_dna_CTCGAAGGAA-CGAGGCCTAT_L001_R1_001_TRIM.fq -2 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto_3_45_dna_CTCGAAGGAA-CGAGGCCTAT_L001_R2_001_TRIM.fq -1 ../DATA_PD/MGTA/TRIMMED/PAIRED/T60_TAT2020_1.fq -2 ../DATA_PD/MGTA/TRIMMED/PAIRED/T60_TAT2020_2.fq -1 ../DATA_PD/MGTA/TRIMMED/PAIRED/T82_TAT2020_1.fq -2 ../DATA_PD/MGTA/TRIMMED/PAIRED/T82_TAT2020_2.fq

#check headers
samtools view -H coverm-genome.Pto10_GGTCACGA-GTATTATG_L001_R1_001_TRIM.fq.bam > header.sam

#si no coincide con el header de la referencia

#modificar contenido de header.sam
# @HD     VN:1.6  SO:coordinate
# @PG     ID:minimap2     PN:minimap2     VN:2.27-r1193   CL:minimap2 --split-prefix /tmp/coverm-minimap2-split-indexXA4ri3 -a -x sr -t 90 /tmp/coverm-mapping-index.KNZFrjjFlqZI/coverm-concatenated-fasta7HKXpt ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto10_GGTCACGA-GTATTATG_L001_R1_001_TRIM.fq ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto10_GGTCACGA-GTATTATG_L001_R2_001_TRIM.fq
# @PG     ID:samtools     PN:samtools     PP:minimap2     VN:1.18 CL:samtools sort -T /tmp/coverm_fifo.JzjpONg0eR84/coverm-make-samtools-sortynM5CN -l0 -@ 89
# @PG     ID:samtools.1   PN:samtools     PP:samtools     VN:1.18 CL:samtools view -F4 -@ 89 -b -o NC_010175_MGTA/coverm-genome.Pto10_GGTCACGA-GTATTATG_L001_R1_001_TRIM.fq.bam
# @PG     ID:samtools.2   PN:samtools     PP:samtools.1   VN:1.20 CL:samtools view -H coverm-genome.Pto10_GGTCACGA-GTATTATG_L001_R1_001_TRIM.fq.bam
# @SQ     SN:NC_010175    LN:5258541

parallel --jobs 20 'samtools reheader header.sam {} > {}.rehead' ::: *bam

rm *bam && rename 's/.rehead//' *rehead

#sort
for f in *bam; do samtools sort -@ 30 -o ${f%.bam}.sorted.bam $f; done

#index
for f in *sorted.bam; do samtools index -@ 30 $f; done

#annotations gtf to gff3
parallel --jobs 4 ' agat_convert_sp_gxf2gxf.pl -g {} -o {}.gff3' ::: *gtf

#edit gff3 solo dejar 1 comment tline
nano gff2metaSNV_annotation.py

##gff-version 3
NC_010175       RefSeq  gene    336     1775    .       +       .       ID=CAUR_RS00005;gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;gene_id=CAUR_RS00005;locus_tag=CAUR_RS00005;old_locus_tag=Caur_0001;transcript_id=""
NC_010175       AGAT    mRNA    336     1775    .       +       .       ID=unassigned_transcript_1;Parent=CAUR_RS00005;Ontology_term=GO:0006270,GO:0006275,GO:0003677,GO:0003688,GO:0005524;exon_number=1;gbkey=CDS;gene=dnaA;gene_id=CAUR_>NC_010175       AGAT    exon    336     1775    .       +       .       ID=agat-exon-61;Parent=unassigned_transcript_1;Ontology_term=GO:0006270,GO:0006275,GO:0003677,GO:0003688,GO:0005524;exon_number=1;gbkey=CDS;gene=dnaA;gene_id=CAUR_>NC_010175       Protein Homology        CDS     336     1775    .       +       0       ID=agat-cds-1;Parent=unassigned_transcript_1;Ontology_term=GO:0006270,GO:0006275,GO:0003677,GO:0003688,GO:0005524;exon_number=1;gbkey=CDS;gene=dnaA>NC_010175       Protein Homology        start_codon     336     338     .       +       0       ID=agat-start_codon-1;Parent=unassigned_transcript_1;Ontology_term=GO:0006270,GO:0006275,GO:0003677,GO:0003688,GO:0005524;exon_number=1;gbk>NC_010175       Protein Homology        stop_codon      1773    1775    .       +       0       ID=agat-stop_codon-1;Parent=unassigned_transcript_1;Ontology_term=GO:0006270,GO:0006275,GO:0003677,GO:0003688,GO:0005524;exon_number=1;gbke>NC_010175

#copiar
cp ~/mambaforge/envs/metaSNV/share/metasnv-2.0.4/src/gff2metaSNV_annotation.py .

# editar lineas en gff2...py
# these paths need to to be adjusted
input_gff            = r"/media/issotta/SSDWORK/loutraki/reference_genomes/ncbi_dataset/data/dbann/NC_010175.gff3"
output_folder        = r"/media/issotta/SSDWORK/loutraki/reference_genomes/ncbi_dataset/data/dbann/NC_010175"


#########################
### Execution metaSNV ###
#########################

#activcate env
mamba activate metaSNV


#faidx references
samtools faidx NC_010175.fna

##contenido file list_bams.txt
#rutas a los bams sorted & con index creado sin lineas en blanco o da error
/media/antiguo_stark/Dilanaz/TEST_metaSNV/VAN18202NCBI_S8-2019_sorted.bam
/media/antiguo_stark/Dilanaz/TEST_metaSNV/VAN18202NCBI_S82020S_sorted.bam
/media/antiguo_stark/Dilanaz/TEST_metaSNV/VAN18202NCBI_S82020_sorted.bam


#run metaSNV
metaSNV.py --threads 30  NC_010175_metaSNV NC_010175_bam.list references/NC_010175.fna --db_ann metaSNV_anntotations/NC_010175_metaSNV_anntotations.txt

metaSNV_Filtering.py --n_threads 30 NC_010175_metaSNV

metaSNV_DistDiv.py --n_threads 30 --dist --div --divNS --matched --filt NC_010175_metaSNV/filtered/pop/

metaSNV_subpopr.R -i NC_010175_metaSNV -p 12
