import argparse
import os

import pandas as pd

import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from barrelseq import config
from barrelseq import sample

# calculate tpms
# need gene length and other gene info
# need genbank annotation file merged in
# need genbank file added to the config

# Should this be added to the config somehow...? I use this list
# to filter several times so I pulled it out here instead of passing
# it around between functions
GENE_INFO_COLS = ['product','type','gene_symbol','locus',
                  'start_coord','end_coord','note','translation'
]

def make_locus_info_dict(gb_file):
    '''
    Given a genbank file, load it into a genbank object, parse out
    the relevant gene information and return as a dictionary of dictionaries
    { locus_tag: {
        gene_symbol: foo,
        product: bar,
        ...
        }
    }
    '''
    # load genbank file into object
    gb_rec = SeqIO.parse(gb_file, "genbank").__next__()

    # initialize a dictionary where the locus tag is the key. The values
    # will also be dictionaries of gene info
    locus_dict = {}
    # for every feature
    for feature in gb_rec.features:
        # start gene info dict for this specific locus_tag
        locus_info = {}
        # only consider features with a locus tag (ignores things like
        # "regulatory" or "repeat region")
        if 'locus_tag' in feature.qualifiers:
            locus_tag = feature.qualifiers['locus_tag'][0]

            locus_info['start_coord'] = feature.location.start.real
            locus_info['end_coord'] = feature.location.end.real
            locus_info['type'] = feature.type
                
            # only consider CDS features
            if feature.type == 'CDS':
                # extract relevant info about genes
                locus_info['locus'] = gb_rec.id
                locus_info['gene_symbol'] = "" if 'gene' not in feature.qualifiers else feature.qualifiers['gene'][0]
                locus_info['product'] = feature.qualifiers['product']
                locus_info['note'] = feature.qualifiers['note']
                locus_info['translation'] = str(feature.translate(gb_rec).seq)
                
            # add this info to the bigger dictionary
            locus_dict[locus_tag] = locus_info

    return locus_dict

def try_get_locus(locus_tag, col_name, locus_dict):
    '''
    Since not every locus tag is a CDS, ignore the ones that don't
    have extra product/gene symbol info available
    '''
    try:
        return locus_dict[locus_tag][col_name]
    except(KeyError):
        # return an empty string for cols that weren't avaiable in the locus_dict
        return ""

def add_gene_info_to_master(master_df, gb_file):
    '''
    Given the master df with raw counts that have already been extracted,
    now merge in other gene info data from the genbank file.
    '''

    # load the genbank file
    locus_dict = make_locus_info_dict(gb_file)

    # add columns to master_df
    for col in GENE_INFO_COLS:
        master_df[col] = master_df['locus_tag'].apply(lambda x: try_get_locus(x,col,locus_dict))
        #master_df[col] = master_df['locus_tag'].apply(lambda x: locus_dict[x][col])

    # also, add gene_len as a column
    master_df['gene_len'] = master_df['end_coord'] - master_df['start_coord'] +1

def calculate_tpm(original_df, sample_cols):
    '''
    Given a dataframe of genes by experimental samples containing read counts,
    calculate the transcripts per million (TPM)
    
    from: https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
    1.) Divide the read counts by the length of each gene in kilobases. This 
        gives you reads per kilobase (RPK).
    2.) Count up all the RPK values in a sample and divide this number by 
        1,000,000. This is your “per million” scaling factor.
    3.) Divide the RPK values by the “per million” scaling factor. This gives 
        you TPM.
    '''
    # copy the original df
    df = original_df.copy()
    
    # loop through the sample columns to create each RPK then TPM column
    rpk_cols = []
    tpm_cols = []
    for col in sample_cols:
        # make new column strings and add to respective list
        rpk_col = col + "_rpk"
        rpk_cols.append(rpk_col)
        tpm_col = col + "_tpm"
        tpm_cols.append(tpm_col)
        
        # calculate reads per kilobase by dividing each column's by the gene length
        # in kilobases (len/1000)
        df[rpk_col] = df[col]/(df['gene_len']/1000)
        
        # calculate rpk scale factor by summing column and divding by 1,000,000
        rpk_scale_factor = df[rpk_col].sum()/float(1000000)
        
        # divde each rpk val by the scale factor
        df[tpm_col] = df[rpk_col]/rpk_scale_factor
    
    # initialze a new TPM only df with the gene info columns
    tpm_df = df[['locus_tag']+GENE_INFO_COLS+tpm_cols]
    
    # Right now, the columns all have "_tpm" appended to them... maybe this 
    # is a bad idea. Happy to change it if that's annoying but maybe it will 
    # be a helpful reminder that this is in TPM space?
    return tpm_df, tpm_cols

def extract(args):
    '''
    After running htseq-count on a set of specified RNA-seq samples, 
    combine all of their summary count files into one master dataframe
    of gene rows and experiment sample columns.
    '''

    # validate config and retrieve specified samples (use all by default)
    cfg = config.validate(args)
    samples = sample.validate(cfg)

    os.chdir(os.path.join(cfg.project_dir, "workspace"))

    # If sample names were provided from the command line, check if
    # it exists in the config
    if args.samples is not None:
        for smp in args.samples:
            # raise error if sample name provided was not in the config
            if smp not in samples['name'].values:
                raise KeyError(f"Invalid sample name", f"The sample {smp} is not contained in the list of valid samples.")
        
        # filter to only the samples passed in the args
        samples = samples.loc[samples['name'].isin(args.samples)]

    # initialize empty df to collect experiment sample columns
    master_df = None
    # initialize a row counter for catching changes in upstream software
    rows = 0

    # for every sample
    for smp in samples['name'].values:
        # load up its summary count file
        summary_file = os.path.join(smp, smp + '.summary.dat')
        
        # If this is the first sample being added, initiialize df
        if master_df is None:
            master_df = pd.read_csv(summary_file, sep='\t', names=['locus_tag', smp])
            rows = master_df.shape[0]
        
        # otherwise, merge this new sample in with the master df
        else:
            # make a temp df with the new experiment file
            tmp_df = pd.read_csv(summary_file, sep='\t', names=['locus_tag', smp])
            
            # check if this experiment has a different number of rows in the summary
            # (if so, something probably changed in upstream counting software)
            if rows != tmp_df.shape[0]:
                raise ValueError("Data frame size mismatch", f"The number of rows in the master data frame is {rows} and the number of rows in the data frame for smp {smp} contains {tmp_df.shape[0]}")
            
            # merge in the new experiment, joining on the gene names
            master_df = master_df.merge(tmp_df, how='outer', on='locus_tag')
            
            # check again for a row count mismatch, just to be sure
            if rows != master_df.shape[0]:
                raise(ValueError, "Data frame size mismatch", f"The number of rows in the old master data frame is {rows} and the number of rows in the merged data frame including smp {smp} contains {master_df.shape[0]}")

    # drop locus rows like "__no_feature", "__ambiguous", "__too_low_quality" from df
    actual_loci = [x for x in master_df['locus_tag'] if not x.startswith("__")]
    master_df = master_df.loc[master_df['locus_tag'].isin(actual_loci)]

    # add gene info to master df
    add_gene_info_to_master(master_df, cfg.reference_gb_path)

    # convert to tpm if requested
    if args.values == "TPM":
        master_df, tpm_cols = calculate_tpm(master_df, samples['name'].values)

    # save final merged df a tsv
    print(f'saving output to {args.output}')
    master_df.to_csv(args.output, sep='\t', header=True, index=False)

    return
