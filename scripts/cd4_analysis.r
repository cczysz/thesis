files_dir = "/group/stranger-lab/nicolel/mRNA_expression/CEL_files/CD4/"
out_dir = paste(files_dir,'out',sep='')

# Normalization files
phenotype_file = "CD4.Samples.POSTQC.ImmVarFinal.txt"
raw_exp_rfile = "cd4_cau_raw.RData"
norm_exp_rfile = "cd4_cau_expr.RData"

# PEER files
peer_factors_file = "cd4_cau_factors.txt"
residuals_file = "cd4_lm.peer_residuals.txt"

# DE files
probeinfo = "/home/t.cczysz/HuGeneProbeInfo.csv"
