#install.packages("/nas40t2/egamazon/matrixeqtl/MatrixEQTL_1.6.2.tar.gz", repos = NULL, type="source")
#library(MatrixEQTL)

#.libPaths("/userhome/genegateRPackages/lib.3.1.1")
library(MatrixEQTL)
#source('/nas40t2/pevans/eQTL/MatrixEQTL/R/Matrix_eQTL_engine.R')

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR ; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
# SNP_file_name = "/group/stranger-lab/GTEx/Thyroid_dosage_euroids.txt"; 
#snps_location_file_name = "/nas40t2/egamazon/gtex/processed_egamazon/snp_location_info";
# snps_location_file_name = "/group/stranger-lab/gtex/genotypes/OMNI_arrays/GTEx_dosage_maflt5_names_snppositons.txt";
# SNP_file_name= "/group/stranger-lab/nicolel/mRNA_expression/CEL_files/CD4/cd4_expr.txt"
# snps_location_file_name = "/group/stranger-lab/nicolel/genotypes/both_eur.cd4.allchr.dosage.filt.eqtl.tsv"
snps_location_file_name = "/group/stranger-lab/nicolel/genotypes/snpspos_snphits_cd4"
SNP_file_name = "/group/stranger-lab/nicolel/genotypes/snppos_eur.cd4.allchr.dosage.filt.eqtl.mod.tsv"
# Gene expression file name
# expression_file_name = "/group/stranger-lab/GTEx/Thyroid_euroids_renormexpr_20p.txt";
expression_file_name = "/group/stranger-lab/nicolel/mRNA_expression/CEL_files/CD4/cd4_expr.txt"
gene_location_file_name = "/group/stranger-lab/nicolel/mRNA_expression/CEL_files/CD4/for_eqtl/genepos_cd4_mod"

# Covariates file name
# Set to character() for no covariates
covariates_file_name =  character();

# Output file name
# output_file_name_cis = "/group/stranger-lab/Thyroid_gtex_ciseqtl_residuals_20p.txt";
output_file_name_cis = "/group/stranger-lab/nicolel/mRNA_expression/CEL_Files/CD4/cd4_ciseqtl.txt"
# output_file_name_tra = "/group/stranger-lab/Thyroid_gtex_transeqtl_residuals_20p.txt";
output_file_name_tra = "/group/stranger-lab/nicolel/mRNA_expression/CEL_Files/CD4/cd4_treqtl.txt"

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1e-3;
pvOutputThreshold_tra = 0;

#trixEQTL Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

cisDist = 1e6;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t"; # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1; # one row of column labels
snps$fileSkipColumns = 1; # one column of row labels
snps$fileSliceSize = 2000; # read file in pieces of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t"; # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1; # one row of column labels
gene$fileSkipColumns = 1; # one column of row labels
gene$fileSliceSize = 2000; # read file in pieces of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t"; # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1; # one row of column labels
cvrt$fileSkipColumns = 1; # one column of row labels
cvrt$fileSliceSize = 2000; # read file in one piece
if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

snpspos = data.frame(snpspos[,-4])

me = Matrix_eQTL_main(
        snps = snps,
        gene = gene,
        cvrt = cvrt,
        output_file_name = output_file_name_tra,
        pvOutputThreshold = pvOutputThreshold_tra,
        useModel = useModel,
        errorCovariance = errorCovariance,
        verbose = TRUE,
        output_file_name.cis = output_file_name_cis,
        pvOutputThreshold.cis = pvOutputThreshold_cis,
        snpspos = snpspos,
        genepos = genepos,
        cisDist = cisDist,
        pvalue.hist = "qqplot"
)

## Plot the Q-Q plot of local and distant p-values

pdf('/group/stranger-lab/nicolel/CEL_Files/CD4/qqplot_residuals_20p.pdf')
plot(me)
dev.off()
