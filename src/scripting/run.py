import os
import itertools

def run_gapmer(
    db_size,
    sim_method,
    k,
    gaps,
    step
):
    outpath = (f"../../tests/outputs/"
        f"data_{db_size}/"
        f"gapmer/{sim_method}/"
        f"k_{k}/"
        f"gaps_{gaps}/"
        f"step_{step}/"
    )
    print(f"Creating output directory...")
    os.system(f"mkdir -p {outpath}")
    command = (
        "cargo run --bin generate_comparison_data -- "
        f"-i ../../tests/inputs/sequences_{db_size}.csv "
        f"-o {outpath}data.csv "
        f"--representation-method gapmer "
        f"--similarity-method {sim_method} "
        f"-k {k} "
        f"--gaps {gaps} "
        f"--step {step} "
    )
    print(f"Executing: {command}")
    os.system(command)

def strobemer1():
    db_size = 2001
    sim_method = "jaccard_similarity"
    order = 3
    strobe_length = 7
    strobe_window_gap = 5
    strobe_window_length = 15
    step = 4

    run_strobemer(
        db_size=db_size,
        sim_method=sim_method,
        order=order,
        strobe_length=strobe_length,
        strobe_window_gap=strobe_window_gap,
        strobe_window_length=strobe_window_length,
        step=step
    )

def gapmer1():
    db_size = 2001
    sim_method = "jaccard_similarity"
    k = 10
    gaps = 1
    step = 4

    run_gapmer(
        db_size=db_size,
        sim_method=sim_method,
        k=k,
        gaps=gaps,
        step=step
    )

def run_alignment(db_size):
    outpath = (f"../../tests/outputs/data_{db_size}/alignment/")
    print(f"Creating output directory...")
    os.system(f"mkdir -p {outpath}")
    command = (
        "cargo run --bin generate_comparison_data -- "
        f"-i ../../tests/inputs/sequences_{db_size}.csv "
        f"-o {outpath}data.csv "
        f"--representation-method direct_alignment "
        f"--similarity-method a "
    )
    print(f"Executing: {command}")
    os.system(command)


def run_kmer(
    db_size,
    sim_method,
    k,
    step
):
    outpath = (f"../../tests/outputs/"
        f"data_{db_size}/"
        f"kmer/{sim_method}/"
        f"k_{k}/"
        f"step_{step}/"
    )
    print(f"Creating output directory...")
    os.system(f"mkdir -p {outpath}")
    command = (
        "cargo run --bin generate_comparison_data -- "
        f"-i ../../tests/inputs/sequences_{db_size}.csv "
        f"-o {outpath}data.csv "
        f"--representation-method kmer "
        f"--similarity-method {sim_method} "
        f"-k {k} "
        f"--step {step} "
    )
    print(f"Executing: {command}")
    os.system(command)

def kmers(k, step):
    run_kmer(
        db_size=2001,
        sim_method="jaccard_similarity",
        k=k,
        step=step
    )


def run_spaced_kmer(
    db_size,
    sim_method,
    k,
    spaces,
    step
):
    outpath = (f"../../tests/outputs/"
        f"data_{db_size}/"
        f"spaced_kmer/{sim_method}/"
        f"k_{k}/"
        f"spaces_{spaces}/"
        f"step_{step}/"
    )
    print(f"Creating output directory...")
    os.system(f"mkdir -p {outpath}")
    command = (
        "cargo run --bin generate_comparison_data -- "
        f"-i ../../tests/inputs/sequences_{db_size}.csv "
        f"-o {outpath}data.csv "
        f"--representation-method spaced_kmer "
        f"--similarity-method {sim_method} "
        f"-k {k} "
        f"--spaces {spaces} "
        f"--step {step} "
    )
    print(f"Executing: {command}")
    os.system(command)

def spaced_kmer():
    for k in [5, 6]:
        for spaces in [2]:
            run_spaced_kmer(2001, "jaccard_similarity", k, spaces, 1)



def run_strobemer(
    db_size,
    sim_method,
    order,
    strobe_length,
    strobe_window_gap,
    strobe_window_length,
    step
):
    outpath = (f"../../tests/outputs/"
        f"data_{db_size}/"
        f"strobemer/{sim_method}/"
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
        f"-o {outpath}data.csv "
        f"--representation-method strobemer "
        f"--similarity-method {sim_method} "
        f"--order {order} "
        f"--strobe-length {strobe_length} "
        f"--strobe-window-gap {strobe_window_gap} "
        f"--strobe-window-length {strobe_window_length} "
        f"--step {step} "
    )
    print(f"Executing: {command}")
    os.system(command)

def strobemer():
    for order in [2]:
        for strobe_length in [5]:
            for strobe_window_gap in [0]:
                for strobe_window_length in [40]:
                    for step in [2, 3]:
                        run_strobemer(2001, "jaccard_similarity", order, strobe_length, strobe_window_gap, strobe_window_length, step)



def run_minimizer(
    db_size,
    sim_method,
    k,
    w,
    step
):
    outpath = (f"../../tests/outputs/"
        f"data_{db_size}/"
        f"kmer/{sim_method}/"
        f"k_{k}/"
        f"w_{w}/"
        f"step_{step}/"
    )
    print(f"Creating output directory...")
    os.system(f"mkdir -p {outpath}")
    command = (
        "cargo run --bin generate_comparison_data -- "
        f"-i ../../tests/inputs/sequences_{db_size}.csv "
        f"-o {outpath}data.csv "
        f"--representation-method kmer "
        f"--similarity-method {sim_method} "
        f"-k {k} "
        f"--minimizer-window-length {w} "
        f"--step {step} "
    )
    print(f"Executing: {command}")
    os.system(command)
