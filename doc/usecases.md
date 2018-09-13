Use cases:
* Generate default configuration file `generate-config`
    * Arguments:
        * output file
        * path to bwa - file must exist
        * path to samtools - file must exist
        * path to htseq-count - file must exist
        * path to R - file must exist
        * reference human readable name
        * path to reference gff - file must exist
        * path to reference nt fasta - file must exist
        * pair-ended (T/F)
        * pass-through options for bwa mem
        * pass-through options for htseq-count
        * pass-through options for samtools view (sam to bam)
        * pass-through options for samtools sort (sort bam)
        * pass-through options for samtools index (index the sorted bam)
        *
        *
        *
* Validate configuration file `validate-config`
    * Arguments:
        * input config file - file must exist
* Add a sample to the sample set `add-sample`
    * Arguments:
        * input config file - file must exist
        * sample name (uMax_run1) - must start with letter, no dashes or periods (use underscores)
        * sample group (uMax) - must start with letter, no dashes or periods (use underscores)
        * path to fastq(s) - file must exist
        * 
        *
        *
        *
        *
        *
* Run alignments and compute read attribution `run`
    * Arguments:
        * input config file
        * processes to use
        *
        *
        *
        *
* Perform statistical analysis and extract table `stats-analysis`
    * for specified subset of samples or all
    * Arguments:
        * input config file
        *
        *
        *
        *
        * (optional) regroup: create a new group with samples separated by ,
        *
* Extract data table `extract-data`
    * for specified subset of samples of all
	* as TPM, raw counts, RPKM
    * Arguments:
        * input config file
* Sample vs. sample parity plots (optional)
