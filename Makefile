SMALL_DB = "tests/inputs/sequences_3.csv"
BIG_DB = "tests/inputs/sequences_2001.csv"
SMALL_OUT = "tests/outputs/data_3_"
BIG_OUT = "tests/outputs/data_2001_"
GEN_COMP = cargo run --bin generate_comparison_data --
MINIMIZER = minimizer_l2_norm


small_kmer_eucidean_distance:
	$(GEN_COMP) -s $(SMALL_DB) -e kmer_euclidean_distance -o $(SMALL_OUT)KED_k6_step3.csv --k 6 --step 3

big_kmer_eucidean_distance:
	$(GEN_COMP) -s $(BIG_DB) -e kmer_euclidean_distance -o $(BIG_OUT)KED_k6_step3.csv --k 6 --step 3

small_kmer_cosine_similarity:
	$(GEN_COMP) -s $(SMALL_DB) -e kmer_cosine_similarity -o $(SMALL_OUT)KCS_k6_step3.csv --k 6 --step 3

big_kmer_cosine_similarity:
	$(GEN_COMP) -s $(BIG_DB) -e kmer_cosine_similarity -o $(BIG_OUT)KCS_k6_step3.csv --k 6 --step 3

small_minimizer_euclidean_distance:
	$(GEN_COMP) -s $(SMALL_DB) -e minimizer_euclidean_distance -o $(SMALL_OUT)MED_k5_window20_step10.csv -m 20 --step 10 --k 5

big_minimizer_euclidean_distance:
	$(GEN_COMP) -s $(BIG_DB) -e minimizer_euclidean_distance -o $(BIG_OUT)MED_k5_window20_step10.csv -m 20 --step 10 --k 5

small_strobemer_euclidean_distance:
	$(GEN_COMP) -s $(SMALL_DB) -e strobemer_euclidean_distance -o $(SMALL_OUT)SED_order3_strobe4_gap2_window10.csv --strobemer-order 3 --strobe-length 4 --strobe-window-gap 2 --strobe-window-length 10 

big_strobemer_euclidean_distance:
	$(GEN_COMP) -s $(BIG_DB) -e strobemer_euclidean_distance -o $(BIG_OUT)SED_order3_strobe4_gap2_window10.csv --strobemer-order 3 --strobe-length 4 --strobe-window-gap 2 --strobe-window-length 10 
