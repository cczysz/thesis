library(oligo)
setwd('/group/stranger-lab/moliva/ImmVar/')
load("Robjects/phen.Robj")

load.cel.files <- function(population,cell_type) {

	if ( cell_type == "CD14" ) {
		
		## Load EUR CD14 CEL files
		cell_type_markers="CD14+16-Mono";
	} 
	else if ( cell_type == "CD4" ) {
		## Load EUR CD4 CEL files
		cell_type_markers="CD3-CD14+CD16-";
	} 
	else {
		stop("Cell type needs to be CD14 or CD4")
	}
	if ( !population%in%c("Caucasian","African-American","Asian") ) {
		stop("Population needs to be Caucasian African-American or Asian")
	}

	files_sufix=as.character(phen[phen$Race == population & phen$CellType == cell_type_markers,"FileName"])
	filenames=unlist(lapply(files_sufix,function(x) paste(c("data/mRNAexp/CEL_files",cell_type,x),collapse="/")));
	raw_oligo=read.celfiles(filenames=filenames,phenoData=AnnotatedDataFrame(phen[phen$Race == population & phen$CellType == cell_type_markers,]))
	return(raw_oligo)
}

backcorrect.normalize.probe.level <- function(ExpressionFeatureSet) {

	bgCorrected <- backgroundCorrect(ExpressionFeatureSet);
	normalized2 <- normalize(bgCorrected, method="quantile");
	return(normalized2);
}
