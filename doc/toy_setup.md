Instructions for setting up a toy sized project.

```
PROJECT_DIR=/tmp/mytoy
mkdir $PROJECT_DIR
mkdir $PROJECT_DIR/data
echo "fq1" > $PROJECT_DIR/data/fq1
echo "fq2" > $PROJECT_DIR/data/fq2
./shim.py config create --config-file $PROJECT_DIR/config.yaml --project-name 'Toy' --project-dir $PROJECT_DIR --reference-name 'Toy genome' --reference-gff-path $PROJECT_DIR/data/reference.gff --reference-fasta-path $PROJECT_DIR/data/reference.fasta  --bwa-path /work/software/bwa/bin/bwa --samtools /work/software/samtools/bin/samtools --htseq-count-path /work/software/htseq/bin/htseq-count --R-path=/home/dacb/anaconda3/bin/R
./shim.py sample add --config-file $PROJECT_DIR/config.yaml --name sample1 --group group1 --fastq-files $PROJECT_DIR/data/fq1 --description 'first sample' 
./shim.py sample add --config-file $PROJECT_DIR/config.yaml --name sample2 --group group2 --fastq-files $PROJECT_DIR/data/fq1 $PROJECT_DIR/data/fq2 --description 'second sample'
```
