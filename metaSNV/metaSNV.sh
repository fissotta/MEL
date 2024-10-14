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
sed -i 's/\.fasta//g' *_J.fasta
#FR with coverm
coverm genome --methods length covered_bases count mean relative_abundance rpkm --min-read-aligned-length 45 --min-read-percent-identity 97 --min-covered-fraction 0 --discard-unmapped --output-format sparse --genome-fasta-files references/NC_009767.fna --genome-fasta-extension fna --bam-file-cache-directory NC_009767_MGTA --threads 90 --output-file NC_009767_MGTA.txt -1 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto10_GGTCACGA-GTATTATG_L001_R1_001_TRIM.fq -2 ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto10_GGTCACGA-GTATTATG_L001_R2_001_TRIM.fq

#check headers
samtools view -H coverm-genome.S16_2020_1.fq.gz.bam > header.sam && nano header.sam

#si no coincide con el header de la referencia

#modificar contenido de header.sam
# @HD     VN:1.6  SO:coordinate
# @PG     ID:minimap2     PN:minimap2     VN:2.27-r1193   CL:minimap2 --split-prefix /tmp/coverm-minimap2-split-indexXA4ri3 -a -x sr -t 90 /tmp/coverm-mapping-index.KNZFrjjFlqZI/coverm-concatenated-fasta7HKXpt ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto10_GGTCACGA-GTATTATG_L001_R1_001_TRIM.fq ../DATA_PD/MGTA/TRIMMED/PAIRED/Pto10_GGTCACGA-GTATTATG_L001_R2_001_TRIM.fq
# @PG     ID:samtools     PN:samtools     PP:minimap2     VN:1.18 CL:samtools sort -T /tmp/coverm_fifo.JzjpONg0eR84/coverm-make-samtools-sortynM5CN -l0 -@ 89
# @PG     ID:samtools.1   PN:samtools     PP:samtools     VN:1.18 CL:samtools view -F4 -@ 89 -b -o NC_010175_MGTA/coverm-genome.Pto10_GGTCACGA-GTATTATG_L001_R1_001_TRIM.fq.bam
# @PG     ID:samtools.2   PN:samtools     PP:samtools.1   VN:1.20 CL:samtools view -H coverm-genome.Pto10_GGTCACGA-GTATTATG_L001_R1_001_TRIM.fq.bam
# @SQ     SN:NC_010175    LN:5258541

parallel --jobs 20 'samtools reheader header.sam {} > {}.rehead' ::: *bam

rm *bam && rename 's/.rehead//' *rehead && rm header.sam

#sort
for f in *bam; do samtools sort -@ 30 -o ${f%.bam}.sorted.bam $f; done

#index
for f in *sorted.bam; do samtools index -@ 30 $f; done

#prodigal
for file in *.fasta; do prodigal -i "$file" -o "${file%.fasta}.gff" -f gff -a "${file%.fasta}_proteins.faa"; done
sed -i '/Model/d' *gff
sed -i '/Sequence/d' *gff

#annotations gtf to gff3
parallel --jobs 4 ' agat_convert_sp_gxf2gxf.pl -g {} -o {}.gff3' ::: *gtf

#edit gff3 solo dejar 1 comment tline
nano gff2metaSNV_annotation.py

#example
samtools view -H coverm-genome.S16_2020_1.fq.gz.bam > header.sam && nano header.sam
parallel --jobs 10 'samtools reheader header.sam {} > {}.rehead' ::: *bam
rm *bam && rename 's/.rehead//' *rehead && rm header.sam && cd ../ && ls

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
metaSNV.py --threads 96  NC_010175_metaSNV NC_010175_bam.list references/NC_010175.fna --db_ann metaSNV_anntotations/NC_010175_metaSNV_anntotations.txt

metaSNV_Filtering.py --n_threads 30 NC_010175_metaSNV

###LB FIX##
awk '{ delete_line=0; for (i=2; i<=NF; i++) if (($i+0) > 0 && ($i+0) < 0.05) { delete_line=1; break } } !delete_line'  S3_2020_Vamb_5.filtered.freq > S3_2020_Vamb_5.filtered_LB.freq

metaSNV_DistDiv.py --n_threads 30 --dist --div --divNS --matched --filt NC_010175_metaSNV/filtered/pop/

metaSNV_subpopr.R -i NC_010175_metaSNV -p 12




#########################3
#########################3


#coverM
for f in *fna; do echo coverm genome --methods length covered_bases count mean relative_abundance rpkm --min-read-aligned-length 45 --min-read-percent-identity 97 --min-covered-fraction 0 --discard-unmapped --output-format sparse --genome-fasta-files $f --genome-fasta-extension fna --bam-file-cache-directory bam_${f%.fna} --threads 96 --output-file ${f%.fna}.tsv -1 ../FASTQ_MG/S4_2019_1.fq.gz -2 ../FASTQ_MG/S4_2019_2.fq.gz -1 ../FASTQ_MG/S6_2019_1.fq.gz -2 ../FASTQ_MG/S6_2019_2.fq.gz -1 ../FASTQ_MG/S7_2020_1.fq.gz -2 ../FASTQ_MG/S7_2020_2.fq.gz -1 ../FASTQ_MG/S8_2019_1.fq.gz -2 ../FASTQ_MG/S8_2019_2.fq.gz -1 ../FASTQ_MG/S16_2020_1.fq.gz -2 ../FASTQ_MG/S16_2020_2.fq.gz -1 ../FASTQ_MG/S17_2020_1.fq.gz -2 ../FASTQ_MG/S17_2020_2.fq.gz; done

#get gff
for f in *fna; do prodigal -i $f -f gff -o ${f%.fna}.gff -a ${f%.fna}.faa; done

#fix header
samtools view -H coverm-genome.S16_2020_1.fq.gz.bam > header.sam && nano header.sam
parallel --jobs 10 'samtools reheader header.sam {} > {}.rehead' ::: *bam
rm *bam && rename 's/.rehead//' *rehead && rm header.sam && cd ../ && ls

#sort
for f in bam_S*/*bam; do samtools sort -@ 96 -o ${f%.bam}.sorted.bam $f; done

#index
for f in bam_S*/*sorted.bam; do samtools index -@ 96 $f; done

#make list
for f in bam*; do ls -d $PWD/$f/*sorted.bam > $f\_list_bams.txt; done

#run shit 1
for f in *fna; do metaSNV.py --threads 96 ${f%.fna}\_metaSNV bam_${f%.fna}\_list_bams.txt $f; done

#run shit 2
for f in *_metaSNV; do metaSNV_Filtering.py --n_threads 96 $f; done

#run shit 3 (whitot nd/ds now so no gff)
#metaSNV_DistDiv.py --n_threads 90 --dist --div --divNS --matched --filt S4_2019_Vamb_2_J_metaSNV/filtered/pop/
for f in *_metaSNV; do metaSNV_DistDiv.py --n_threads 90 --dist --div --matched --filt $f/filtered/pop/; done
