
import sys
import scipy.stats as stats

#input_fp = sys.argv[1]
input_fp = 'TF_SNP_Motif_count'
# perform Fisher's excat test
with open('Fisher_excat_test.txt' ,'w') as out:
    out.write('\t'.join(str(s) for s in ['TF_name','odds ratio', 'p-value', 'Dis_SNP', 'Dis_Motif','Con_SNP', 'Con_Motif'])+'\n')
    with open(input_fp, 'U') as f:
        while True:
            tf = f.readline().strip()
            if not tf: break
            Dis_SNP = int(f.readline().strip())
            Dis_Motif = int(f.readline().strip())
            Con_SNP = int(f.readline().strip())
            Con_Motif = int(f.readline().strip())
            s_oddsratio, s_pvalue = stats.fisher_exact([[Dis_SNP, Dis_Motif],[Con_SNP, Con_Motif]])
            print tf, s_oddsratio, s_pvalue, Dis_SNP, Dis_Motif,Con_SNP, Con_Motif
            out.write('\t'.join(str(s) for s in [tf, s_oddsratio, s_pvalue, Dis_SNP, Dis_Motif,Con_SNP, Con_Motif])+'\n')
            