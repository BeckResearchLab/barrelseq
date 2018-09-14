# barrelseq
Bacteria &amp; Archaea, Rna-seq, Resequenceing (EL) for isolate and simple community analysis pipeline
--

Dependencies:
* pyyaml
* biopython
* bwa
* samtools
* htseq
* R
    * deseq2

Until `setup.py` is complete, you shuld run the tool with `shim.py` in this directory. E.g.
```
./shim.py --help
./shim.py config create --config-file /tmp/b --project-name 'Test' --project-dir /tmp --reference-name /tmp/a --reference-gff-path /tmp/a --reference-fasta-path /tmp/a  --bwa-path /work/software/bwa/bin/bwa --samtools /work/software/samtools/bin/samtools --htseq-count-path /work/software/htseq/bin/htseq-count --R-path=/home/dacb/anaconda3/bin/R --pair-ended
./shim.py config edit --config-file /tmp/b --project-name 'TestA'
```
