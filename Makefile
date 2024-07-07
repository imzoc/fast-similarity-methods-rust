GEN_COMP = cargo run --bin generate_comparison_data --

small_kmer_eucidean_distance:
	DB_SIZE = 3
	DB = tests/inputs/sequences_$(DB_SIZE).csv
	METHOD = kmer_euclidean_distance
	K = 6
	STEP = 3
	OUTPATH = tests/outputs/data_$(DB_SIZE)/$(METHOD)/$(K)/$(STEP)/data.csv
	mkdir -p $(OUTPATH)
	$(GEN_COMP) -s $(SMALL_DB) -e METHOD -o $(SMALL_OUT)KED_k6_step3.csv --k 6 --step 3
