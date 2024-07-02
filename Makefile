SEQ_3 = "tests/inputs/sequences_3.csv"
SEQ_1K = "tests/inputs/sequences_1000.csv"
GEN_COMP = cargo run --bin generate_comparison_data --
MINIMIZER = minimizer_l2_norm

smol:
	$(GEN_COMP) -s $(SEQ_3) -o tests/outputs/out3.csv

big:
	$(GEN_COMP) -s $(SEQ_1K) -o tests/outputs/out1k.csv

smol_minimizer:
	$(GEN_COMP) -s $(SEQ_1K) -e $(MINIMIZER) -o tests/outputs/out1k_minimizer.csv
