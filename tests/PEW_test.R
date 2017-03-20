sourceDir <- function(path, trace = TRUE, ...) {
	for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
		if(trace) cat(nm,":")
		source(file.path(path, nm), ...)
		if(trace) cat("\n")
	}
}

sourceDir("~/GitHub/chipenrich-weights/R/")

chipenrich(peaks = "./ENCODEpeaks/wgEncodeAwgTfbsHaibGm12878Pax5c20Pcr1xUniPk.narrowPeak",
		   out_name = 'PEWtest_wgEncodeAwgTfbsHaibGm12878Pax5c20Pcr1xUniPk',
		   out_path = ".", 
		   genome = "hg19", 
		   genesets = "GOBP", 
		   locusdef = "nearest_tss", 
		   method = "polyenrich_weighted",
		   qc_plots = F, 
		   max_geneset_size = 2000, 
		   n_cores = 1)

justpeaks(peaks = "./ENCODEpeaks/wgEncodeAwgTfbsHaibH1hescJundV0416102UniPk.narrowPeak",
		  out_name = 'PEWtest_wgEncodeAwgTfbsHaibH1hescJundV0416102UniPk',
		  out_path = ".", 
		  genome = "hg19", 
		  genesets = "GOBP",
		  locusdef = "nearest_tss",
		  method = "polyenrich_weighted",
		  qc_plots = F, 
		  max_geneset_size = 2000, 
		  n_cores = 1)

justenrich(ppg = read.table("./PEWtest_peaks-per-gene.tab",header=T),
		   out_name = 'PEWtest_wgEncodeAwgTfbsHaibGm12878Pax5c20Pcr1xUniPk',
		   out_path = ".", 
		   genome = "hg19", 
		   genesets = "GOBP", 
		   locusdef = "nearest_tss", 
		   method = "polyenrich_weighted",
		   qc_plots = F, 
		   max_geneset_size = 2000, 
		   n_cores = 1)

justenrich(ppg_none,
		   out_name = 'PEWtest_none_wgEncodeAwgTfbsHaibGm12878Pax5c20Pcr1xUniPk',
		   out_path = ".", 
		   genome = "hg19", 
		   genesets = "GOBP", 
		   locusdef = "nearest_tss", 
		   method = "polyenrich_weighted",
		   qc_plots = F, 
		   max_geneset_size = 2000, 
		   n_cores = 1)
justenrich(ppg_log,
		   out_name = 'PEWtest_log_wgEncodeAwgTfbsHaibGm12878Pax5c20Pcr1xUniPk',
		   out_path = ".", 
		   genome = "hg19", 
		   genesets = "GOBP", 
		   locusdef = "nearest_tss", 
		   method = "polyenrich_weighted",
		   qc_plots = F, 
		   max_geneset_size = 2000, 
		   n_cores = 1)
justenrich(ppg_loglog,
		   out_name = 'PEWtest_loglog_wgEncodeAwgTfbsHaibGm12878Pax5c20Pcr1xUniPk',
		   out_path = ".", 
		   genome = "hg19", 
		   genesets = "GOBP", 
		   locusdef = "nearest_tss", 
		   method = "polyenrich_weighted",
		   qc_plots = F, 
		   max_geneset_size = 2000, 
		   n_cores = 1)
