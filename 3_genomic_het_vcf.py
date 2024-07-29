#!/usr/bin/env python3

# Author: Jesus Castrejon
# Script to count call rate and Heterogocity from VCF file

# Outputs to Std Output (see -h)

import polars as pl
import argparse
import sys

def get_names(vcf_path):
    with open(vcf_path, "rt") as ifile:
        for line in ifile:
            if line.startswith("#CHROM"):
                vcf_names = [x for x in line.split('\t')]
                break
    ifile.close()
    ll=[]
    for name in vcf_names:
        if '#' in name:
            name = name.split('#')[1]
        elif '\n' in name:
            name = name.split('\n')[0]
        ll.append(name)
    return ll



def parse_dataframe(myfile):

    names = get_names(myfile)

    df = pl.scan_csv(myfile, comment_prefix='#', separator="\t", has_header=False,)
    cols = df.collect_schema().names()
    df = df.rename(dict(zip(cols,names)))

    idx = names.index('FORMAT')
    inds = names[idx+1:]

    df = df.select(pl.col(x) for x in inds).with_columns( pl.col(x).str.split(':').list.first() for x in inds)
    dfu = df.unpivot()
    dfu = dfu.with_columns(
            pl.when(pl.col("value") == './.')
            .then(0)
            .when(pl.col("value").str.split('/').list.n_unique() > 1)
            .then(2)
            .otherwise(1).alias("value")).group_by(
            ['variable','value']).len().collect()

    dfu = dfu.pivot("value", index="variable", values="len").sort('variable').rename({'2':'#het','1':'#hom','0':'#NAs','variable':'individuals'})
    dfu = dfu.with_columns((pl.col('#het') + pl.col('#hom')).alias('#call'))
    return dfu

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A Python script to count heterogocity and call rate from a VCF dataframe')
    parser.add_argument('filename', help='The path of the VCF file')           # positional argument
    args = vars(parser.parse_args())

    df = parse_dataframe(args['filename'], )
    df.write_csv(sys.stdout, separator='\t')