# export rare variants burden
```bash
python pgen2txt.py -g chr1 -o ./ -s keep1.snpid -b -n PCSK9i
python pgen2txt.py -g chr6 -o ./ -s keep6.snpid -b -n LPAi
```

# replicate
```bash
python replicate.py -g LPA-burden.txt,PCSK9-burden.txt -p lpa.score,ldl.sscore -w ./ --assoc-file cov.txt.gz --phenos PAD,AAA,AA --covariates age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,smoking -res PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --keep keep.eid --markername lpa,ldl_a --rvas --score --quant
```
