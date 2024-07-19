import os
import itertools

# Define the parameter values
A_values = [2001, 3]  # database size
B_values = ["strobemer"]  # representation method
X_values = ["euclidean_distance", "jaccard_similarity"] # distance function
C_values = [2, 3]  # order
D_values = [7, 8]  # strobe_length
E_values = [5]  # strobe_window_gap
F_values = [10, 15]  # strobe_window_length
G_values = [2, 4]  # step

# Generate all combinations of the parameters
combinations = list(itertools.product(A_values, B_values, X_values, C_values, D_values, E_values, F_values, G_values))

# Base command
base_command = "cargo run --bin generate_comparison_data -- "

for combo in combinations:
    A, B, X, C, D, E, F, G = combo
    outpath = f"../../tests/outputs/data_{A}/{B}/{X}/order_{C}/strobe_length_{D}/strobe_window_gap_{E}/strobe_window_length_{F}/step_{G}/"

    print(f"Creating output directory...")
    os.system(f"mkdir -p {outpath}")
    command = (
        f"{base_command}"
        f"-s ../../tests/inputs/sequences_{A}.csv "
        f"-r {B} "
        f"-d {X} "
        f"--order {C} "
        f"--strobe-length {D} "
        f"--strobe-window-gap {E} "
        f"--strobe-window-length {F} "
        f"--step {G} "
        f"-o {outpath}data.csv"
    )
    print(f"Executing: {command}")
    os.system(command)