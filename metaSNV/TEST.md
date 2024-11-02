#Sort BAM files
for f in *bam; do samtools sort -@ 10 -o ${f%.bam}.sorted.bam $f; done

#extract references from bams
for f in *bam; do samtools view -H $f; done | grep SN | sort -u | sed 's/.*SN:\([^\t]*\).*/\1/' > referencias_extract.list

#extract by reference
for f in *sorted.bam; do while read -r l; do echo "samtools view -h -b $f $l > ${f%.sorted.bam}__${l}_filtered.bam"; done < referencias_extract.list; done > split_bams.sh

#move splited
mkdir splited_bam && mv *__* splited_bam

#Extract SAM headers from BAM files
for f in *bam; do samtools view -H $f > ${f%.bam}.head; done

# Modified multiple headers
# Example usage: python filter_sn_lines.py /path/to/files/

import os
import sys

def filter_sn_lines_inplace(directory):
    # Iterate through each file in the directory
    for filename in os.listdir(directory):
        if filename.endswith(".head"):
            filepath = os.path.join(directory, filename)
            target_sn = filename.split("__")[-1].replace("_filtered.head", "")  # Extract SN target based on filename
            
            with open(filepath, 'r') as file:
                lines = file.readlines()

            # Filter lines: Keep @HD and @PG lines, and only the matching @SQ line
            new_lines = [line for line in lines if line.startswith("@HD") or line.startswith("@PG") or 
                         (line.startswith("@SQ") and f"SN:{target_sn}" in line)]

            # Overwrite the original file with the filtered lines
            with open(filepath, 'w') as file:
                file.writelines(new_lines)
            
            print(f"Processed and modified {filepath}")

# Pass the directory path as an argument when running the script
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python filter_sn_lines.py /path/to/files/")
    else:
        filter_sn_lines_inplace(sys.argv[1])


# Reheader BAM files using modified headers
for f in *head; do samtools reheader $f ${f%.head}\.bam > ${f%.head}\.bam_reheaded; done; 

for f in *bam; do samtools view -h -o ${f%bam}sam $f; done
        
#Sort BAM files
for f in *bam; do samtools sort -@ 10 -o ${f%.bam}.sorted.bam $f; done

# Index BAM files
for f in *.bam; do samtools index -@ 30 $f; done
