http://bioinformatics.intec.ugent.be/kmarchal/EVORhA/

perl ~/SCRIPTS/singleline.pl NC_010175.fna > NC_010175.fasta

sed '/^>/! s/[^ACTG]/N/g' NC_010175.fasta > NC_010175_ACTG.fasta
