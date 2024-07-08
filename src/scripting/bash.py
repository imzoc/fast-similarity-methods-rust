import os
import itertools

# Define the parameter values
A_values = [3, 2001]  # database size
B_values = ["strobemer_euclidean_distance"]  # estimation method
C_values = [2, 3]  # order
D_values = [7, 8]  # strobe_length
E_values = [5]  # strobe_window_gap
F_values = [10, 15]  # strobe_window_length
G_values = [2, 4]  # step

# Generate all combinations of the parameters
combinations = list(itertools.product(A_values, B_values, C_values, D_values, E_values, F_values, G_values))

# Base command
base_command = "cargo run --bin generate_comparison_data -- "

for combo in combinations:
    A, B, C, D, E, F, G = combo
    command = (
        f"{base_command}"
        f"-s tests/inputs/sequences_{A}.csv "
        f"-e {B} "
        f"--order {C} "
        f"--strobe_length {D} "
        f"--strobe_window_gap {E} "
        f"--strobe_window_length {F} "
        f"--step {G} "
        f"-o tests/outputs/data_{A}/{B}/order_{C}/strobe_length_{D}/strobe_window_gap_{E}/strobe_window_length_{F}/step_{G}/data.csv"
    )
    print(f"Executing: {command}")
    os.system(command)