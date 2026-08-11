# export rare variants burden
```bash
chmod 777 ./pgen2txt

# exist: chr1.psam, chr1.pgen, chr1.pvar
./pgen2txt -g ukb_wes_imp_chr1 -s ../RVAS/pcsk9.snp -o ./ -b -n PCSK9
./pgen2txt -g ukb_wes_imp_chr6 -s ../RVAS/LPA.snp -o ./ -b -n LPA
```

# replicate
```bash
chmod 777 ./replicate

./replicate --help

./replicate -g LPA-burden.txt,PCSK9-burden.txt -p lpa.score,ldl.sscore -w ./ --assoc-file cov.txt.gz --phenos PAD,AAA,AA --covariates age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,smoking -res PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --keep keep.eid --markername lpa,ldl_a --rvas --score --quant
```

# data format
- First col of each file is used to merge.
- Files of -g and -p have 2 col: eid and burden / eid and score
- --assoc-file has covariates and phenotypes
- variables of -res are used to adjust population struction for genetic score
- --keep: one col, no header, contains sample id
- --markername: biomarkers corresponding to two files of -p or -g
- --rvas, --score, --quant: perform which analysis
