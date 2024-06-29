use anyhow::Result;
use clap::Parser;
use csv::{ReaderBuilder, WriterBuilder};
use serde::Deserialize;
use std::collections::HashMap;

#[derive(Debug, Parser)]
#[command(about, author, version)]
/// Similarity Methods
struct Args {
    #[arg(short, long)]
    sequences: String,

    #[arg(short, long, default_value = "tests/outputs/data.csv")]
    outfile: String,

    #[arg(short, long, default_value = "l2_norm")]
    estimation_method: String,

    #[arg(short, long, default_values_t = [2, 3, 4], num_args(1..=3))] // What does num_args() do?
    k_range: Vec<usize>,

    #[arg(short, long, default_value_t = 10)]
    w: usize,
}

#[allow(unused_imports)]
use similarity_methods::utils::sequence;
use similarity_methods::utils::tensor;

/* Using serde to parse CSV data. */
#[derive(Debug, Deserialize)]
struct DatabaseRecord {
    base_sequence: String,
    modified_sequence: String,
    edit_distance: usize,
}

impl DatabaseRecord {
    fn generate_kmer_parameters(&self, k: &usize, w: &usize) -> tensor::Parameters {
        tensor::Parameters{
            k:k.clone(),
            w:w.clone(),
            base_sequence: self.base_sequence.chars().collect(),
            modified_sequence: self.modified_sequence.chars().collect()
        }
    }
}

/* Custom struct to store data mapped to k and edit distance */
struct KAndEditDistanceHashMap {
    data: HashMap<(usize, usize), f64>,
}

impl KAndEditDistanceHashMap {
    // Create a new MyStruct with an empty HashMap
    fn new() -> KAndEditDistanceHashMap {
        KAndEditDistanceHashMap {
            data: HashMap::new(),
        }
    }

    fn insert(&mut self, k: usize, edit_distance: usize, value: f64) {
        self.data.insert((k, edit_distance), value);
    }

    fn update(&mut self, k: usize, edit_distance: usize, value: f64) {
        let key = (k, edit_distance);
        if let Some(cur) = self.data.get_mut(&key) {
            *cur += value;
        } else {
            self.insert(k, edit_distance, value);
        }
    }

    fn get(&self, k: usize, edit_distance: usize) -> &f64 {
        self.data.get(&(k, edit_distance)).unwrap()
    }
}

// --------------------------------------------------
fn main() {
    if let Err(e) = run(Args::parse()) {
        eprintln!("{e}");
        std::process::exit(1);
    }
}

// --------------------------------------------------
// See this repo's README file for pseudocode
fn run(args: Args) -> Result<()> {
    // Rolling magnitude and variability data.
    let mut edit_distance_sums = KAndEditDistanceHashMap::new();
    let mut edit_distance_squared_sums = KAndEditDistanceHashMap::new();
    let mut edit_distance_counts = KAndEditDistanceHashMap::new();
    let mut edit_distances = vec![];
    let mut rdr = ReaderBuilder::new().from_path(&args.sequences)?;

    for result in rdr.deserialize() {
        let record: DatabaseRecord = result?;
        edit_distances.push(record.edit_distance);

        for k in &args.k_range {
            let params = record.generate_kmer_parameters(k, &args.w);
            let estimated_distance: f64 = match args.estimation_method.as_str() {
                "l2_norm" => tensor::l2norm(params),
                "cosine_similarity" => tensor::cosine_similarity(params),
                "minimizer_l2_norm" => tensor::minimizer_l2_norm(params),
                "strobemer" => tensor::strobemer(params),
                _ => unreachable!(),
            };

            edit_distance_sums.update(
                *k,
                record.edit_distance,
                estimated_distance,
            );

            edit_distance_squared_sums.update(
                *k,
                record.edit_distance,
                estimated_distance * estimated_distance,
            );

            edit_distance_counts.update(*k, record.edit_distance, 1.);
        }
    }

    // Compute mean and confidence intervals from rolling data.
    let mut lower_bounds = KAndEditDistanceHashMap::new();
    let mut upper_bounds = KAndEditDistanceHashMap::new();
    let mut means = KAndEditDistanceHashMap::new();
    for k in &args.k_range {
        for edit_distance in edit_distances.clone() {
            let sum = edit_distance_sums.get(*k, edit_distance);
            let squared_sum =
                edit_distance_squared_sums.get(*k, edit_distance);
            let count = edit_distance_counts.get(*k, edit_distance);

            let mean = sum / count;
            let variance = (squared_sum - mean * sum) / count;
            let mean_se = (variance / count).sqrt();

            means.update(*k, edit_distance, mean);
            lower_bounds.update(*k, edit_distance, mean - mean_se);
            upper_bounds.update(*k, edit_distance, mean + mean_se);
        }
    };

    // Create the header row for comparison data file
    let mut header_row = vec!["edit distance".to_string()];
    for k in &args.k_range {
        header_row.push(format!("mean (k={})", k));
        header_row.push(format!("lower confidence bound (k={})", k));
        header_row.push(format!("upper confidence bound (k={})", k));
    }

    let mut wtr = WriterBuilder::new()
        .delimiter(b',')
        .from_path(&args.outfile)?;
    wtr.write_record(&header_row)?;

    for dist in edit_distances.clone() {
        let mut row = vec![dist.to_string()];
        for k in &args.k_range {
            // mean
            row.push(format!("{}", means.get(*k, dist)));
            row.push(format!("{}", lower_bounds.get(*k, dist)));
            row.push(format!("{}", upper_bounds.get(*k, dist)));
        }
        wtr.write_record(&row)?;
    }

    println!(r#"Done, see output "{}""#, args.outfile);
    Ok(())
}