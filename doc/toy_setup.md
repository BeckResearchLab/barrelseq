Instructions for setting up a toy sized project.

```
PROJECT_DIR=/tmp/$USER/mytoy
mkdir -p $PROJECT_DIR
mkdir $PROJECT_DIR/data
cp barrelseq/data/* $PROJECT_DIR/data
./shim.py config create --config-file $PROJECT_DIR/config.yaml --project-name 'Toy' --project-dir $PROJECT_DIR --reference-name 'Toy genome' --reference-gff-path $PROJECT_DIR/data/reference.gff --reference-fasta-path $PROJECT_DIR/data/reference.fna  --bwa-path /work/software/bwa/bin/bwa --samtools /work/software/samtools/bin/samtools --htseq-count-path /work/software/htseq/bin/htseq-count --R-path=/home/dacb/anaconda3/bin/R
./shim.py sample add --config-file $PROJECT_DIR/config.yaml --name sample1 --group group1 --fastq-files $PROJECT_DIR/data/sample1.fq --description 'first sample' 
./shim.py sample add --config-file $PROJECT_DIR/config.yaml --name sample2 --group group1 --fastq-files $PROJECT_DIR/data/sample2.fq $PROJECT_DIR/data/fq2 --description 'second sample'
```
