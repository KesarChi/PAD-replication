'''
Author: Kesar-Chi && momiji0613@163.com
Date: 2024-06-07 21:09:57
LastEditors: Kesar-Chi && momiji0613@163.com
LastEditTime: 2026-07-13 13:45:57
FilePath: /BIGCode/utils/pgen2txt.py
Description: Institute: Beijing Institute Of Genomics, CAS.

Copyright (c) 2024 by ${momiji0613@163.com}, All Rights Reserved. 
'''


import pandas as pd
import subprocess
from multiprocessing import Pool
import os
from glob import glob
import numpy as np
import pgenlib as pg
import argparse
import warnings
warnings.filterwarnings('ignore')


def run_cmd(args):
    subprocess.run(args, shell=True, check=True)
    

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--geno', type=str, help='pgen file')
    parser.add_argument('-o', '--opath', type=str, help='output path')
    parser.add_argument('-s', '--snp', type=str, help='snp list file, no header')
    parser.add_argument('-b', '--burden', action='store_true', help='save burden score')
    parser.add_argument('-n', '--name', type=str, help='burden name')
    return parser.parse_args()


def run_plink(args):
    cmd = f"plink2 --pfile {args.geno} --extract {args.snp} --make-pgen --out {args.opath}/{os.path.basename(args.snp).split('.')[0]}"
    run_cmd(cmd)
    return f"{args.opath}/{os.path.basename(args.snp).split('.')[0]}"
    
    
def parse_geno(snp_geno, opath):
    pvar = pg.PvarReader(bytes(f"{snp_geno}.pvar", 'utf-8'))
    pgen = pg.PgenReader(bytes(f"{snp_geno}.pgen", 'utf-8'))
    ids = pd.read_csv(f"{snp_geno}.psam", sep=r'\s+')
    # keep first col, keep dataframe format, not series
    ids = ids.iloc[:, [0]]
    ids.columns = ['IID']
    snps, samples = pgen.get_variant_ct(), pgen.get_raw_sample_ct()
    res = pd.DataFrame({'eid':ids['IID']})
    for i in range(snps):
        snpid = pvar.get_variant_id(i).decode('utf-8')
        out = np.array(np.empty(samples), dtype=np.int8)
        pgen.read(i, out)
        res[snpid] = out
    res.to_csv(f"{snp_geno}.txt", sep='\t', index=False)
    print(f"Genotype file has been saved to {snp_geno}.txt")
    pvar.close()
    pgen.close()
    return res


def clean_opath(snp_geno):
    run_cmd(f"rm {snp_geno}.p*")


if __name__ == '__main__':
    args = parse_args()
    if not os.path.exists(args.opath):
        os.makedirs(args.opath)
    
    snp_geno = run_plink(args)
    res = parse_geno(snp_geno, args.opath)
    clean_opath(snp_geno)
    
    if args.burden:
        res = res.replace(-9, 0)
        res = res.fillna(0)
        res[args.name] = res.drop(columns=['eid']).sum(axis=1)
        res[['eid', args.name]].to_csv(f"{args.opath}/{args.name}-burden.txt", sep='\t', index=False)
        print(f"Burden score file has been saved to {args.opath}/{args.name}-burden.txt")

