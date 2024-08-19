#!/usr/bin/python3

"""
Script to add the number of derived alleles per sample to the gerp output file. Written to be run per window, 
by providing contig (or scaffold/chromosome) name, start and end position (1-based) on the command line.

It is based on Verena Kutschera's script for the same purpose but with major improvements in runtime and memory-usage. 

Author: Jesús Castrejón

(tested with polars==1.5)

Usage (gunzipped vcf):
    python3 gerp_derived_alleles.py -gz -g gerp_file -v vcf_file_gzipped -c contig -s start_position -e end_position -o out_file_name 

Usage (unzipped vcf):
    python3 gerp_derived_alleles.py -g gerp_file -v vcf_file_unzipped -c contig -s start_position -e end_position -o out_file_name 
"""

import polars as pl
import gzip as gz
from datetime import datetime
import argparse

def read_gerp_windows(gerpFile, chrom, start, end):
    schema={'column_1':pl.String,'column_2':pl.Int64,'column_3':pl.String,'column_4':pl.Float64}
    names=['#CHROM','POS','ancestral_state','gerp_score']
    cols = ['column_1', 'column_2', 'column_3', 'column_4',]

    # Reads (lazy) gerp file and filters positions
    df = pl.scan_csv(gerpFile,separator='\t',has_header=False, schema=schema).rename(dict(zip(cols,names)))
    df = df.filter((pl.col("#CHROM") == chrom) & (pl.col("POS").is_between(start, end) ) ).collect()

    # Add a row with "NaN" in case the window is missing from the file
    if len(df)==0:
        df = pl.DataFrame({"#CHROM": chrom, 'POS':start,'ancestral_state':float("nan"),'gerp_score':float("nan")})
    return df

def read_vcf_windows(vcf_path, chrom, start, end, gzip=True):
    def get_names(vcf_path, gzip,):
        ifile = gz.open(vcf_path, "rt") if gzip else open(vcf_path, "rt")
        with ifile:
              for line in ifile:
                if line.startswith("#CHROM"):
                    vcf_names = [x for x in line.split('\t')]
                    break
        ifile.close()
        return [x.split('\n')[0] if '\n' in x else x for x in vcf_names]

    names = get_names(vcf_path,gzip=gzip,) 
    inds = names[names.index('FORMAT')+1:]

    if(gzip):
        # Reads (in-memory) vcf file
        with gz.open(vcf_path, "rt") as f:   
            df = pl.read_csv(f,comment_prefix='#', separator="\t", has_header=False,)
        f.close()
        # Filters by position
        df = df.filter((pl.col("column_1") == chrom) & (pl.col("column_2").is_between(start, end) ) )
    else:
        # Reads (lazy) vcf file
        df = pl.scan_csv(vcf_path,comment_prefix='#', separator="\t", has_header=False,)
        # Filters by position
        df = df.filter((pl.col("column_1") == chrom) & (pl.col("column_2").is_between(start, end) ) ).collect()
    
    cols = df.collect_schema().names()
    df = df.rename(dict(zip(cols,names)))

    # Fix the genotype column to only contain data until the first ":" and drop unused columns
    df = df.with_columns( pl.col(x).str.split(':').list.first() for x in inds).drop(["ID","QUAL","FILTER","INFO","FORMAT"])
    # Add a row with "NaN" in case the window is missing from the file
    if len(df)==0:
        d1 = {"#CHROM": chrom, 'POS':start,}
        d1.update({c:float("nan") for c in names[2:]})
        df = pl.DataFrame(d1)
    return df

def count_derived_alleles(gerp_path,vcf_path,chrom, start, end, gzip=False):
    # Loads GERP and VCF files into dataframes
    df_gerp = read_gerp_windows(gerp_path, chrom, start, end)
    df_vcf  = read_vcf_windows(vcf_path, chrom, start, end, gzip=gzip)

    names = df_vcf.columns
    inds = names[names.index('ALT')+1:]

    # Join GERP and VCFs for each window, and drops null rows
    df = df_gerp.join(df_vcf, on=['#CHROM','POS'], how='left').drop_nulls(pl.col('REF'))
    del df_gerp, df_vcf #Drops unused DFs from memory
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "GERP output file and VCF file successfully joined for chromosome", chrom, "from position", start, "to position", end)

    ### Creates temp column find number of alleles different to ancestral allele
    df = df.with_columns(pl.when(pl.col('ancestral_state')==pl.col('REF')).then(pl.lit('1')).otherwise(pl.lit('0')).alias('num_ancestral'))
    df = df.with_columns(
            pl.when(
                pl.col(x).str.contains(r'\.'))
                .then(float("nan"))
                .when(pl.col('ancestral_state')=='N')
                .then(float("nan"))
                .when( pl.col(x).is_null() )
                .then(float("nan"))
                .otherwise(  pl.col(x).str.count_matches(pl.col('num_ancestral'))).alias(x+'_no_derived_alleles')
            for x in inds
        )
    # Drops unused columns and exits function
    [df.drop_in_place(x) for x in inds + ['num_ancestral','REF','ALT']]
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "GERP output file and VCF file successfully processed for chromosome", chrom, "from position", start, "to position", end)
    return df

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A Python script to add the number of derived alleles per sample to the gerp output file.')
    parser.add_argument('-g','--gerp', help='GERP Dataframe',required=True)   # gerp output file with ancestral state and gerp score per site
    parser.add_argument('-v','--vcf', help='VCF Dataframe',required=True)    # VCF file to be merged
    parser.add_argument('-c','--contig', help='Contig to process',required=True) # contig to be processed
    parser.add_argument('-s','--start', help='Start position for window',required=True)  # start position for window (1-based)
    parser.add_argument('-e','--end', help='End position for window',required=True)    # end position for window
    parser.add_argument('-o','--output', help='Output file path',required=True) # output file path
    parser.add_argument('-gz','--gzip', help='Boolean to indicate whether vcf file is gunzip compressed or not (default False)',action='store_true')  

    args = vars(parser.parse_args())
    # Call the function
    df = count_derived_alleles(args['gerp'], args['vcf'], args['contig'], int(args['start']), int(args['end']), gzip=args['gzip'])

    # Write the window to file for further processing    
    df.write_csv(args['output'], separator='\t', null_value='NaN',)
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "Result for chromosome", args['contig'], "from position", int(args['start']), "to position", int(args['end']), "successfully written to", args['output'])