library(affy)
library(BiocInstaller)
biocLite('maskBAD')
library(maskBAD)

cd4_dir = "/group/stranger-lab/nicolel/mRNA_expression/CEL_files/CD4"
cd14_dir = "/group/stranger-lab/nicolel/mRNA_expression/CEL_files/CD14"

cd4_phen = read.table(file=paste(cd4_dir,'CD4.Samples.Caucasian.txt',sep='/'),header=T,as.is=T,sep=',',fill=T)

cd14_phen = read.table(file=paste(cd14_dir,'CD14.Caucasian.txt',sep='/'),header=T,as.is=T,sep=',',fill=T)

cd14_files = paste(cd14_dir,cd14_phen$FileName,sep='/')
cd4_files = paste(cd4_dir,cd4_phen$File.Name,sep='/')

all_files = c(cd14_files,cd4_files)

# cel_data = ReadAffy(filenames=all_files)
# save(list=c('cel_data'),file="/home/t.cczysz/all_cel_data.Robj")
load('/home/t.cczysz/all_cel_data.Robj')

groups = c(numeric(length(cd14_files))+1,numeric(length(cd4_files))+2)

exmask <- mask(affy=cel_data,ind=groups,useExpr=F)
save(list=c('exmask'),file="/home/t.cczysz/mask.Robj")

pdf('scores_hist.pdf')
hist(exmask$probes$quality.scores)
dev.off()


