import polars as pl
import argparse

def parse_piawka_het(filename):
	# Name output files
	fname = input_filename.split('/')[-1]
	wdir = '/'.join(input_filename.split('/')[:-1])
	if '/' in input_filename:
		wdir +=  '/'
	output_het_filename = wdir   +  'genomic_het_table.txt'

	names = ['locus','nSites','pop1','pop2','nUsed ','metric','value','numerator','denominator',]
	cols = ['column_1', 'column_2', 'column_3', 'column_4', 'column_5','column_6', 'column_7','column_8','column_9','column_10','column_11']

	# Read data into a Lazy DataFrame
	df = pl.scan_csv(input_filename,separator="\t",has_header=False).rename(dict(zip(cols,names)))
	# Filters per metric
	df_het  = df.filter((pl.col('metric')=='het_pixy')).select(['pop1','numerator','denominator',])

	#--------------------------------------------------
	# HET Table
	df_het = df_het.group_by("pop1").agg(pl.col("numerator").sum(), pl.col("denominator").sum(),).sort('pop1')
	df_het = df_het.with_columns((pl.col('numerator') / pl.col('denominator')).alias('het')).collect()#.rename({"numerator":'diffs','denominator':'comps'})
	df_het.write_csv(output_het_filename,separator='\t',)
	del df_het, df
	return


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='A Python script to reduce a pixy HET dataframe.')
	parser.add_argument('filename', help='The path of the dataframe')           # positional argument  
	args = vars(parser.parse_args())
	input_filename = args['filename']
	parse_piawka_het(args['filename'],)