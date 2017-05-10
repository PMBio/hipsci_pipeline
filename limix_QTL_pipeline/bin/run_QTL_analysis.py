#!/usr/bin/env python2
import argparse

def get_args():
    parser = argparse.ArgumentParser(description='QTL pipeline')
    parser.add_argument('-geno_prefix','--geno_prefix', required=True)
    parser.add_argument('-pheno_file','--pheno_file', required=False)
    parser.add_argument('-out_file','--out_file', required=True)
    parser.add_argument('-cis_window_kb','--cis_window_kb', required=True)
    args = parser.parse_args()

    geno_prefix = args.geno_prefix
    pheno_file = args.pheno_file
    out_file = args.out_file
    cis_window_kb = args.cis_window_kb

    return geno_prefix, pheno_file, out_file, cis_window_kb