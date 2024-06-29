SEQ_3 = "tests/inputs/sequences_3.csv"
SEQ_1K = "tests/inputs/sequences_1000.csv"
GEN_COMP = cargo run --bin generate_comparison_data --

smol:
	$(GEN_COMP) -s $(SEQ_3) -o tests/outputs/out3.csv

big:
	$(GEN_COMP) -s $(SEQ_1K) -o tests/outputs/out1k.csv
