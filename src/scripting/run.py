import os
import itertools

# Define the parameter values
db_size = 2001
method = "strobemer"
sim_func = "jaccard_similarity"
order = 3
strobe_length = 7
strobe_window_gap = 5
strobe_window_length = 15
step = 4

outpath = (f"../../tests/outputs/"
    f"data_{db_size}/"
    f"{method}/"
    f"{sim_func}/"
    f"order_{order}/"
    f"strobe_length_{strobe_length}/"
    f"strobe_window_gap_{strobe_window_gap}/"
    f"strobe_window_length_{strobe_window_length}/"
    f"step_{step}/"
)

print(f"Creating output directory...")
os.system(f"mkdir -p {outpath}")
command = (
    "cargo run --bin generate_comparison_data -- "
    f"-i ../../tests/inputs/sequences_{db_size}.csv "
    f"-o {outpath}data.csv"
    f"--representation-method {method} "
    f"--similarity-method {sim_func} "
    f"--order {order} "
    f"--strobe-length {strobe_length} "
    f"--strobe-window-gap {strobe_window_gap} "
    f"--strobe-window-length {strobe_window_length} "
    f"--step {step} "
)
print(f"Executing: {command}")
os.system(command)