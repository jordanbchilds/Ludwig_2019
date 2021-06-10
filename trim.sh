
  ## Trim quality and adaptors based on experiment: edit categories.txt, eg. SRP149536 is scATAC-seq runs ##

readarray -t types < categories.txt
for i in types; do
# trim reads for low quality. REMEMBER: change hisat2 input to _trimmed.fastq for alignment
echo "trimming ${i}";
cutadapt -j 0 -q 30 -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -A GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -o fastq/${i}_2_trimmed.fastq -p fastq/${i}_2_trimmed.fastq fastq/${rt}_1.fastq fastq/${rt}_2.fastq;
# ^ NOT reverse complement of adaptors
done


## Split into nuclear sequences and mitochondrial sequences

