polyenrich = function(
	peaks,
	out_name = "polyenrich",
	out_path = getwd(),
	genome = "hg19",
	genesets = c(
		'GOBP',
		'GOCC',
		'GOMF'),
	locusdef = "nearest_tss",
	method = c('polyenrich','polyenrich_weighted'),
	weight_col = NA,
	use_mappability = F,
	mappa_file = NULL,
	read_length = 36,
	qc_plots = T,
	max_geneset_size = 2000,
	n_cores = 1
)