SUPPORTED_METHODS = list(
  'chipenrich' = "test_gam",
  'chipenrich_fast' = "test_gam_fast",
  'fet' = "test_fisher_exact",
  'broadenrich' = "test_gam_ratio",
  'countenrich' = "test_gam_nb",
  'countenrich_fast' = "test_gam_nb_fast",
  'chipapprox' = "test_approx"
)

HIDDEN_METHODS = list(
  'binomial' = "test_binomial",
  'broadenrich_splineless' = "test_gam_ratio_splineless",
  'polyenrich_weighted' = "test_pew_fast"
)

METHOD_NAMES = list(
	'fet' = "Fisher's Exact Test",
	'chipenrich' = "ChIP-Enrich",
	'broadenrich' = "Broad-Enrich",
	'chipapprox' = "ChIP-Enrich Approximate",
	'broadenrich_splineless' = "Broad-Enrich Splineless",
	'countenrich' = "Poly-Enrich",
	'chipenrich_fast' = "ChIP-Enrich Fast",
    'countenrich_fast' = "Poly-Enrich Fast",
    'polyenrich_weighted' = "Poly-Enrich Weighted"
)
