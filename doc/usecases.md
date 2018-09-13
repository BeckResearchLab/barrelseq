Use cases:
* Generate default configuration file `generate-config`
    * Arguments:
        * output file - must be writable
        * path to bwa - file must exist & be readable
        * path to samtools - file must exist & be readable
        * path to htseq-count - file must exist & be readable
        * path to R - file must exist & be readable
        * path to project directory - directory must exist & be writable
        * reference human readable name
        * path to reference gff - file must exist & be readable
        * path to reference nt fasta - file must exist & be readable
        * pair-ended (T/F)
        * pass-through options for bwa mem
        * pass-through options for htseq-count
        * pass-through options for samtools view (sam to bam)
        * pass-through options for samtools sort (sort bam)
        * pass-through options for samtools index (index the sorted bam)
        *
        *
        *
        *
        *
        *
* Validate configuration file `validate-config`
    * Arguments:
        * input config file - file must exist & be readable
* Add a sample to the sample set `add-sample`
    * compute md5sum for each fastq
    * sample\_info will have separate date fields for added, last\_run, last\_analyzed
    * Arguments:
        * input config file - file must exist & be readable
        * sample name (uMax\_run1) - must start with letter, no dashes or periods (use underscores)
        * sample group (uMax) - must start with letter, no dashes or periods (use underscores)
        * path to fastq(s) - file(s) must exist & be readable
        * full sample description - must not contain tabs or line breaks
        *
* Run alignments and compute read attribution `run`
    * can run the commands directly or make bash files for inspection or other uses
    * Arguments:
        * input config file - file must exist & be readable
        * number of processes to use - must be integer 1>=1, default=1, not passed to bwa
        * name of samples from sample\_info to run - default is all samples with empty last\_run date
        * save as script - default is F, when T shell scripts are written with the commands, nothing is run
        * preserve intermediate results - default is F, when T commands that remove files, e.g. SAM, are not generated
        *
        *
        *
        *
        *
* Perform statistical analysis and extract table `stats-analysis`
    * for specified subset of samples or all
    * Arguments:
        * input config file
        * output name prefix, default = {A}\_vs\_{B}
        * groupA - group must exist in sample info, must have a last run date for all samples in group
        * groupB - group must exist in sample info, must have a last run date for all samples in group
        * (optional) regroupA - create a new group with samples separated by ',', must have last run date for all samples in listed samples, e.g. --regroupA EPS=uMax1,uMax2 would define a new temporary group EPS with samples uMax1, uMax2
        * (optional) regroupB - create a new group with samples separated by ',', must have last run date for all samples in listed samples
        * shrinkage estimator - default is 'apeglm', string must be one of 'normal', 'apeglm', 'ashr'
        * output transformation - default is 'vst', string must be one of 'vst' or 'rlog'
        * produce figures - default is T
        *
* Extract data table `extract-data`
    * for specified subset of samples of all
	* as TPM, raw counts, RPKM
    * Arguments:
        * input config file
        * output name
        * output values, default='TPM', must be string from list 'TPM', 'raw', 'RPKM'
        * output type, default='tsv', must be string from list 'tsv', 'csv', 'xls'
        * sample name list, default is None, otherwise comma separated string of sample names, must all have last run date
        *
        *
        *
        *
* Sample vs. sample parity plots (optional) - possibly deprecated
