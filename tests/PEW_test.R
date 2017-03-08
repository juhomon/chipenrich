sourceDir = function(path) {
	for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
		if(trace) cat(nm,":")           
		source(file.path(path, nm), ...)
		if(trace) cat("\n")
	}
}

sourceDir("~/GitHub/chipenrich-weights/R/")

chipenrich(peaks = "./ENCODEpeaks/wgEncodeAwgTfbsHaibGm12878Mef2csc13268V0416101UniPk.narrowPeak",
		   out_name = 'PEWtest',
		   out_path = ".", 
		   genome = "hg19", 
		   genesets = "GOBP", 
		   locusdef = "nearest_tss", 
		   method = "polyenrich_weighted",
		   qc_plots = F, 
		   max_geneset_size = 2000, 
		   n_cores = 1)