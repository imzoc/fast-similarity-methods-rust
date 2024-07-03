SMALL_DB = "tests/inputs/sequences_3.csv"
BIG_DB = "tests/inputs/sequences_2001.csv"
SMALL_OUT = "tests/outputs/data_3_"
BIG_OUT = "tests/outputs/data_2001_"
GEN_COMP = cargo run --bin generate_comparison_data --
MINIMIZER = minimizer_l2_norm


smol_l2norm:
	$(GEN_COMP) -s $(SMALL_DB) -e l2_norm -o $(SMALL_OUT)l2norm.csv

big_l2norm:
	$(GEN_COMP) -s $(BIG_DB) -e l2_norm -o $(BIG_OUT)l2norm.csv

smol_minimizer:
	$(GEN_COMP) -s $(SMALL_DB) -e $(MINIMIZER) -o $(SMALL_OUT)minimizer.csv

big_minimizer:
	$(GEN_COMP) -s $(BIG_DB) -e $(MINIMIZER) -o $(BIG_OUT)minimizer.csv

smol_strobemer:
	$(GEN_COMP) -s $(SMALL_DB) -e strobemer -o $(SMALL_OUT)strobemer.csv

smol_cosine:
	$(GEN_COMP) -s $(SMALL_DB) -e cosine_similarity -o $(SMALL_OUT)cosine_similarity.csv

big_cosine:
	$(GEN_COMP) -s $(BIG_DB) -e cosine_similarity -o $(BIG_OUT)cosine_similarity.csv