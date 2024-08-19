# Population genomics

A collection of very efficient scripts written in python using Polars.

## Installation

Install the required dependencies with:

```
pip install -r requirements.txt
```

## Description and Usage

The scrips provided are:

-  01_genomic_piawka_pi_dxy_fst.py: script to summarise the output of PIAWK (PI). It produces as output tables and matrix for genomic pi, dxy, and fst statistics.

-  02_genomic_piawka_het.py: script to summarise the output of PIAWK (HET). It produces as output a table with genomic_het statistics.

-  03_genomic_het_vcf.py: script to count heterogocity and call rate directly from a VCF file.

- 04_gerp_derived_alleles.py: script to add the number of derived alleles per sample to the gerp output file.

For further information call the help option. For example:

```
python 04_gerp_derived_alleles.py -h
```

## Contact

jesus.castrejon.fig@gmail.com