'''
Author: Kesar-Chi && momiji0613@163.com
Date: 2026-07-13 11:45:17
LastEditors: Kesar-Chi && momiji0613@163.com
LastEditTime: 2026-07-13 19:39:27
FilePath: /BIGCode/project/GeneHunt/LPA/LPA-paper2/replicate.py
Description: Institute: Beijing Institute Of Genomics, CAS.

Copyright (c) 2026 by ${momiji0613@163.com}, All Rights Reserved. 
'''
import pandas as pd
import numpy as np
import statsmodels.api as sm
import numpy as np
import os
import argparse
import warnings
warnings.filterwarnings('ignore')


def parse_args():
    parser = argparse.ArgumentParser(description="Replicate the results of LPA-paper2")
    parser.add_argument('-g', "--geno", help="burden files, comma split", type=str)
    parser.add_argument('-p', "--prs", help="score files, comma split, first is for lpa, second is for ldl-c, two cols, eid and score, score col needs name", type=str)
    parser.add_argument('-w', "--workdir", required=True, help="workdir", type=str)
    parser.add_argument("--assoc-file", required=True, help="tab split, pheno+covariates")
    parser.add_argument("--phenos", default="", help="逗号分隔的二元表型")
    parser.add_argument("--covariates", default="age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10", type=str)
    parser.add_argument('-res', "--residual", default="PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10", type=str)
    parser.add_argument("--keep", help="keep file, one col, no header", type=str)
    parser.add_argument("--markername", default="Lpa,ldl", help="col names for lpa and ldl-c levels", type=str)
    parser.add_argument("--rvas", default=False, action="store_true", help="perform rvas?")
    parser.add_argument("--score", default=False, action="store_true", help="perform score analysis?")
    parser.add_argument("--quant", default=False, action="store_true", help="perform quant analysis?")
    return parser.parse_args()


def read_and_merge(arg, names, merge_key="eid"):
    dfs = []
    files = arg.split(',')
    for f in files:
        df = read_table(f)
        dfs.append(df)

    lpaf, ldlf = dfs[0], dfs[1]
    lpaf.columns = ['eid', names[0]]
    ldlf.columns = ['eid', names[1]]
    lpaf[merge_key] = lpaf[merge_key].astype(str)
    ldlf[merge_key] = ldlf[merge_key].astype(str)
    merged_df = pd.merge(lpaf, ldlf, on=merge_key, how="inner")

    return merged_df


def catx_y_counts(data, y, x, cov):
    if isinstance(x, str):
        x = [x]
    features = x + cov if cov is not None else x
    test = data[features + [y]].dropna().reset_index(drop=True)
    for i in x:
        # tmp = test.loc[test[i]==1, y].value_counts().to_dict()
        # num1, num2 = tmp[1], tmp[0]
        print(f"{i}:{test.loc[test[i]==1, y].value_counts().to_dict()}")
    return {i:test.loc[test[i]==1, y].value_counts().to_dict() for i in x}
    
    
def parse_case_control(res):
    parsed = {}
    for k, v in res.items():
        n_case = 0
        n_control = 0
        for key, count in v.items():
            if int(key) == 1:
                n_case += count
            elif int(key) == 0:
                n_control += count
        parsed[k] = {"N_case": n_case, "N_control": n_control}
    return parsed


def read_table(file):
    ext = os.path.basename(file).lower()
    if ext.endswith(".gz"):
        compression = "gzip"
        ext = ext[:-3]
    else:
        compression = None

    if ext.endswith((".txt", ".tsv")):
        sep = "\t"
    elif ext.endswith(".csv"):
        sep = ","
    else:
        sep = "\t"

    return pd.read_csv(file, sep=sep, compression=compression)


def get_resid(data, y, x):
    data = data.reset_index(drop=True)
    data = data[x+[y]]
    data = data.astype(float)
    X = sm.add_constant(data[x])
    model = sm.OLS(data[y], X)
    result = model.fit().resid
    result += data[y].mean()
    return result


def logit_model(data, y, x, cov):
    if isinstance(x, str): 
        x = [x]
    features = x + cov if cov is not None else x
    df = data[features + [y]].dropna().reset_index(drop=True)
    X = df[features]
    X = sm.add_constant(X)
    Y = df[y]
    # logistic regression
    model = sm.Logit(Y, X.astype(float))
    result = model.fit()
    return result


def RVAS(df, phenos, covariates, masks, root='./'):
    df['Both carrier'] = (df[masks[0]] == 1) & (df[masks[1]] == 1)
    df[f'only {masks[0]}'] = (df[masks[0]] == 1) & (df[masks[1]] == 0)
    df[f'only {masks[1]}'] = (df[masks[0]] == 0) & (df[masks[1]] == 1)
    df['Non carrier'] = (df[masks[0]] == 0) & (df[masks[1]] == 0)

    groups = ['Non carrier', f'only {masks[1]}', f'only {masks[0]}', 'Both carrier']
    print("Carrier counts:")
    for g in groups:
        print(f"{g}: {df[g].value_counts().to_dict()}")
    print()
    
    analysis(df, groups, covariates, phenos, f"{root}/Result1-RVAS.tsv", "Non carrier")


def Score(df, phenos, covariates, masks, root='./'):
    print("Score analysis:")
    df[f'{masks[0]}_q'] = pd.qcut(df[masks[0]], q=2, labels=[1, 2])
    df[f'{masks[1]}_q'] = pd.qcut(df[masks[1]], q=2, labels=[1, 2])
    df = pd.get_dummies(df, columns=[f'{masks[0]}_q', f'{masks[1]}_q'])
    df[f'{masks[0]}-low-{masks[1]}-low'] = df[f'{masks[0]}_q_1'] & df[f'{masks[1]}_q_1']
    df[f'{masks[0]}-low-{masks[1]}-high'] = df[f'{masks[0]}_q_1'] & df[f'{masks[1]}_q_2']
    df[f'{masks[0]}-high-{masks[1]}-low'] = df[f'{masks[0]}_q_2'] & df[f'{masks[1]}_q_1']
    df[f'{masks[0]}-high-{masks[1]}-high'] = df[f'{masks[0]}_q_2'] & df[f'{masks[1]}_q_2']  
    
    groups = [f'{masks[0]}-high-{masks[1]}-high', f'{masks[0]}-high-{masks[1]}-low', f'{masks[0]}-low-{masks[1]}-high', f'{masks[0]}-low-{masks[1]}-low']
    print("Carrier counts:")
    for g in groups:
        print(f"{g}: {df[g].value_counts().to_dict()}")
    print()
    
    analysis(df, groups, covariates, phenos, f"{root}/Result2-ScoreMedian.tsv", f'{masks[0]}-high-{masks[1]}-high')


def Quant(df, pheno, covariates, output, pair1=['Lpas', 'lpa'], pair2=['LDLs', 'ldl']):
    print(f"Quant analysis for {pheno}:")
    res_lpa_percent = []
    res_ldl_percent = []
    res_beta = []
    res_beta_lpal = []
    res_beta_pcsk9l = []
    res_diff_lpa = []
    res_diff_ldl = []
    res_p = []
    res_p_lpa = []
    res_p_ldl = []
    res_se = []
    res_se_lpa = []
    res_se_ldl = []

    for i in range(46):
        top1_lpa = df[pair1[0]].quantile(0.5 + 0.01*i)
        bot1_lpa = df[pair1[0]].quantile(0.5 - 0.01*i)
        top1_lpa_rows = df[df[pair1[0]] >= top1_lpa].index
        bot1_lpa_rows = df[df[pair1[0]] <= bot1_lpa].index
        df['top-lpa'] = 0
        df.loc[top1_lpa_rows, 'top-lpa'] = 1
        df['bot-lpa'] = 0
        df.loc[bot1_lpa_rows, 'bot-lpa'] = 1
        for j in range(46):
            res_lpa_percent.append(0.5 - 0.01*i)
            res_ldl_percent.append(0.5 - 0.01*j)
            top1_pcsk9 = df[pair2[0]].quantile(0.5 + 0.01*j)
            bot1_pcsk9 = df[pair2[0]].quantile(0.5 - 0.01*j)
            top1_pcsk9_rows = df[df[pair2[0]] >= top1_pcsk9].index
            bot1_pcsk9_rows = df[df[pair2[0]] <= bot1_pcsk9].index
            df['top-pcsk9'] = 0
            df.loc[top1_pcsk9_rows, 'top-pcsk9'] = 1
            df['bot-pcsk9'] = 0
            df.loc[bot1_pcsk9_rows, 'bot-pcsk9'] = 1
            
            df['T-LPA-T-PCSK9'] = df['top-lpa'] * df['top-pcsk9']
            df['T-LPA-B-PCSK9'] = df['top-lpa'] * df['bot-pcsk9']
            df['B-LPA-T-PCSK9'] = df['bot-lpa'] * df['top-pcsk9']
            df['B-LPA-B-PCSK9'] = df['bot-lpa'] * df['bot-pcsk9']
            newtest = df[(df['T-LPA-T-PCSK9']==1) | (df['T-LPA-B-PCSK9']==1) | (df['B-LPA-T-PCSK9']==1) | (df['B-LPA-B-PCSK9']==1)].reset_index(drop=True)
            res_diff_lpa.append(newtest.loc[newtest['top-lpa']==1, pair1[1]].mean() - newtest.loc[newtest['bot-lpa']==1, pair1[1]].mean())
            res_diff_ldl.append(newtest.loc[newtest['top-pcsk9']==1, pair2[1]].mean() - newtest.loc[newtest['bot-pcsk9']==1, pair2[1]].mean())
            
            model = logit_model(newtest, pheno, ['B-LPA-B-PCSK9', 'B-LPA-T-PCSK9', 'T-LPA-B-PCSK9'], covariates)
            res_beta.append(model.params['B-LPA-B-PCSK9'])
            res_beta_lpal.append(model.params['B-LPA-T-PCSK9'])
            res_beta_pcsk9l.append(model.params['T-LPA-B-PCSK9'])
            res_se.append(model.bse['B-LPA-B-PCSK9'])
            res_se_lpa.append(model.bse['B-LPA-T-PCSK9'])
            res_se_ldl.append(model.bse['T-LPA-B-PCSK9'])
            res_p.append(model.pvalues['B-LPA-B-PCSK9'])
            res_p_lpa.append(model.pvalues['B-LPA-T-PCSK9'])
            res_p_ldl.append(model.pvalues['T-LPA-B-PCSK9'])
            
    sv = pd.DataFrame({'lpa-percent': res_lpa_percent, 'ldl-percent': res_ldl_percent, 'lpa-diff': res_diff_lpa, 'ldl-diff': res_diff_ldl, 'beta-BB': res_beta, 'beta-lpaB': res_beta_lpal, 'beta-pcsk9B': res_beta_pcsk9l, 'p-BB': res_p, 'p-lpaB': res_p_lpa, 'p-pcsk9B': res_p_ldl, 'se-BB': res_se, 'se-lpaB': res_se_lpa, 'se-pcsk9B': res_se_ldl})
    sv.to_csv(output, index=False, sep='\t')
    print(f"Save results to {output}")


def analysis(df, groups, covariates, phenos, output, ref):
    results = []
    for p in phenos:
        print(f"============={p} =============")
        cnt = parse_case_control(catx_y_counts(df, p, groups, covariates))
        model_vars = groups[1:]
        result = logit_model(df, y=p, x=model_vars, cov=covariates)
        conf = result.conf_int()
        for _, var in enumerate(groups):
            row = {"Trait": p, "Vars": var, "Ncase": cnt[var]["N_case"], "Ncontrol": cnt[var]["N_control"], "OR": np.nan, "Pval": np.nan, "L-OR": np.nan, "H-OR": np.nan}
            if var == ref:
                row["OR"] = 1
            else:
                row["OR"] = np.exp(result.params[var])
                row["Pval"] = result.pvalues[var]
                row["L-OR"] = np.exp(conf.loc[var, 0])
                row["H-OR"] = np.exp(conf.loc[var, 1])
            results.append(row)

    results = pd.DataFrame(results)
    results.to_csv(f"{output}", sep="\t", index=False)
    print(f"Save results to {output}")


if __name__ == "__main__":
    args = parse_args()
    cov = read_table(args.assoc_file)
    cov.rename(columns={cov.columns[0]: "eid"}, inplace=True)
    cov['eid'] = cov['eid'].astype(str)
    
    phenos = args.phenos.split(',')
    cov_list = args.covariates.split(',')
    residual_list = args.residual.split(',')
    
    if args.keep:
        cov = cov[cov['eid'].isin(pd.read_csv(args.keep, header=None)[0].astype(str))]
        print(f"Keep {cov.shape[0]} samples")
    
    if args.rvas:
        maskr = ['LPAi', 'PCSK9i']
        geno = read_and_merge(args.geno, maskr)
        geno.rename(columns={geno.columns[0]: "eid"}, inplace=True)
        geno['eid'] = geno['eid'].astype(str)
        geno['LPAi'] = geno['LPAi'].clip(0, 1)
        geno['PCSK9i'] = geno['PCSK9i'].clip(0, 1)
        geno = geno.merge(cov, on="eid", how="inner")
        RVAS(geno, phenos, cov_list, maskr, args.workdir)
        
    if args.score:
        masks = ['Lpas', 'LDLs']
        score = read_and_merge(args.prs, masks)
        score.rename(columns={score.columns[0]: "eid"}, inplace=True)
        score['eid'] = score['eid'].astype(str)
        score = score.merge(cov, on="eid", how="inner")
        score[f'{masks[0]}'] = get_resid(score, masks[0], residual_list)
        score[f'{masks[1]}'] = get_resid(score, masks[1], residual_list)
        Score(score, phenos, cov_list, masks, args.workdir)
        if args.quant:
            mname = args.markername.split(',')
            p1, p2 = [masks[0], mname[0]], [masks[1], mname[1]]
            for p in phenos:
                Quant(score, p, cov_list, f"{args.workdir}/Result3-{p}.tsv", p1, p2)
