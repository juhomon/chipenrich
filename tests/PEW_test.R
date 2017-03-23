sourceDir <- function(path, trace = TRUE, ...) {
	for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
		if(trace) cat(nm,":")
		source(file.path(path, nm), ...)
		if(trace) cat("\n")
	}
}

devtools::load_all("~/GitHub/chipenrich-weights/")

sourceDir("~/GitHub/chipenrich-weights/R/")

peaks = read.table("./PEWtest_wgEncodeAwgTfbsHaibH1hescJundV0416102UniPk_peaks.tab",header=T)
peaks_log = peaks; peaks_log$peak_weight = log(peaks_log$peak_weight)
peaks_loglog = peaks_log; peaks_loglog$peak_weight = log(peaks_loglog$peak_weight)

ppg_log = calc_genes_peak_weight(peaks_log, ppg)

chipenrich(peaks = "./ENCODEpeaks/wgEncodeAwgTfbsHaibH1hescJundV0416102UniPk.narrowPeak",
		   out_name = 'PEWtest_wgEncodeAwgTfbsHaibH1hescJundV0416102UniPk.narrowPeak',
		   out_path = ".", 
		   genome = "hg19", 
		   genesets = "GOBP", 
		   locusdef = "nearest_tss", 
		   method = "polyenrich_weighted",
		   qc_plots = F, 
		   max_geneset_size = 2000, 
		   n_cores = 1)

chipenrich(peaks = "./ENCODEpeaks/wgEncodeAwgTfbsSydhK562JundIggrabUniPk.narrowPeak",
		   out_name = 'PE_wgEncodeAwgTfbsSydhK562JundIggrabUniPk.narrowPeak',
		   out_path = ".", 
		   genome = "hg19", 
		   genesets = "GOBP", 
		   locusdef = "nearest_tss", 
		   method = "countenrich_fast",
		   qc_plots = F, 
		   max_geneset_size = 2000, 
		   n_cores = 2)

chipenrich(peaks = "./ENCODEpeaks/wgEncodeAwgTfbsSydhK562JundIggrabUniPk.narrowPeak",
		   out_name = 'PEW_none_wgEncodeAwgTfbsSydhK562JundIggrabUniPk.narrowPeak',
		   out_path = ".", 
		   genome = "hg19", 
		   genesets = "GOBP", 
		   locusdef = "nearest_tss", 
		   method = "polyenrich_weighted",
		   qc_plots = F, 
		   max_geneset_size = 2000, 
		   n_cores = 2)



justpeaks(peaks = "./ENCODEpeaks/wgEncodeAwgTfbsSydhK562JundIggrabUniPk.narrowPeak",
		  out_name = 'PE_wgEncodeAwgTfbsSydhK562JundIggrabUniPk.narrowPeak',
		  out_path = ".", 
		  genome = "hg19", 
		  genesets = "GOBP",
		  locusdef = "nearest_tss",
		  method = "polyenrich_weighted",
		  qc_plots = F, 
		  max_geneset_size = 2000, 
		  n_cores = 1)

justenrich(ppg,
		   out_name = 'PEWtest_PE_wgEncodeAwgTfbsHaibGm12878Pax5c20Pcr1xUniPk',
		   out_path = ".", 
		   genome = "hg19", 
		   genesets = "GOBP", 
		   locusdef = "nearest_tss", 
		   method = "countenrich_fast",
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

plot_spline_counts = function(ppg, variable) {
	gpw = ppg[order(ppg$log10_length),]
	
	fitspl = gam(as.formula(sprintf("%s~s(log10_length,m=3,bs='cr')",variable)),data=gpw,family="nb")
	as.numeric(predict(fitspl, gpw, type="terms"))->gpw$spline
	
	spline_table = matrix(0,nrow=ceiling(nrow(gpw)/25),ncol=2)
	for (i in 1:(nrow(spline_table)-1)) {
		spline_table[i,] = c(mean(gpw$log10_length[(25*(i-1)+1):(25*i)]),
							 mean(gpw[(25*(i-1)+1):(25*i),variable]))
	}
	spline_table[nrow(spline_table),] = c(mean(gpw$log10_length[(25*(nrow(spline_table)-1)+1):nrow(gpw)]),
										  mean(gpw[(25*(nrow(spline_table)-1)+1):nrow(gpw),variable]))
	
	plot(spline_table[,1],spline_table[,2],xlab="log10 locus length", ylab="Y: counts",pch=19,cex=0.3)
	lines(gpw$log10_length, gpw$spline - mean(gpw$spline) + mean(gpw[,variable]) ,col = "red", lwd=2)
	
}

plot_spline_counts_nobins = function(ppg, variable) {
	gpw = ppg[order(ppg$log10_length),]
	
	fitspl = gam(as.formula(sprintf("%s~s(log10_length,bs='cr')",variable)),data=gpw,family="nb")
	as.numeric(predict(fitspl, gpw, type="terms"))->gpw$spline
	
	plot(gpw$log10_length,gpw[,variable],xlab="log10 locus length", ylab="Y: counts",pch=19,cex=0.3)
	lines(gpw$log10_length, gpw$spline - mean(gpw$spline) + mean(gpw[,variable]) ,col = "red", lwd=2)
	
}


