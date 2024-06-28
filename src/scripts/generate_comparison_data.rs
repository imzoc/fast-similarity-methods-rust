use std::fs::File;
use std::io::{Result, Write,BufReader,BufRead};
use std::collections::HashMap;

use csv::ReaderBuilder;
use serde::Deserialize;

use similarity_methods::utils::tensor;
#[allow(unused_imports)]
use similarity_methods::utils::sequence;

struct Metadata {
    sequence_length: usize,
    min_edit_distance: usize,
    max_edit_distance: usize,
}

/* Using serde to parse CSV data. */
#[derive(Debug, Deserialize)]
struct DatabaseRecord {
    base_sequence: String,
    modified_sequence: String,
    edit_distance: usize,
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

    fn insert(&mut self, k: &usize, edit_distance: &usize, value: f64) {
        let key = (k.clone(), edit_distance.clone());
        self.data.insert(key, value);
    }

    fn update(&mut self, k: &usize, edit_distance: &usize, value: f64) {
        let key = (*k, *edit_distance);
        if self.data.contains_key(&key) {
            *self.data.get_mut(&key).unwrap() += value;
        } else {
            self.insert(k, edit_distance, value);
        }
    }
    
    fn get(&self, k: &usize, edit_distance: &usize) -> &f64 {
        self.data.get(&(*k, *edit_distance)).unwrap()
    }

}
/*
 * This function parses metadata from the database file generated by generate_string_pairs.rs.
 * I'm having trouble with the parse() method, which says that:
 * the trait `From<ParseIntError>` is not implemented for `std::io::Error`. 
 */
fn read_metadata(sequence_database_file_name: &String) -> Result<Metadata> {
    let file = File::open(sequence_database_file_name)?;
    let mut reader = BufReader::new(file);

    let mut first_line = String::new();
    reader.read_line(&mut first_line)?;
    let sequence_length: usize = first_line.trim().strip_prefix("# base string length=")
        .unwrap().parse()?;

    let mut second_line = String::new();
    reader.read_line(&mut second_line)?;
    let min_edit_distance: usize = second_line.trim().strip_prefix("# minimum edit distance=")
        .unwrap().parse()?;

    let mut third_line = String::new();
    reader.read_line(&mut third_line)?;
    let max_edit_distance: usize = third_line.trim().strip_prefix("# maximum edit distance=")
        .unwrap().parse()?;

    let metadata = Metadata {
        sequence_length: sequence_length,
        min_edit_distance: min_edit_distance,
        max_edit_distance: max_edit_distance
    };
    Ok(metadata)
}

/* Please see this repo's README file, it contains the pseudocode for this script.
 */
fn main() -> Result<()> {
    // Filenames and hyper-parameters -- should the CLI handle this?
    let k_range: Vec<usize> = [2,3,4].to_vec();
    let sequence_database_file_name = "../../tests/inputs/sequences_1000.csv".to_string();
    let comparison_data_file_name = "../../tests/inputs/data.csv".to_string();

    // Parse database-related metadata -- Can I do this better?
    let metadata = read_metadata(&sequence_database_file_name)?;
    let _sequence_length = metadata.sequence_length;
    let max_edit_distance = metadata.max_edit_distance;
    let min_edit_distance = metadata.min_edit_distance;

    // Rolling magnitude and variability data.
    let mut edit_distance_sums = KAndEditDistanceHashMap::new();
    let mut edit_distance_squared_sums = KAndEditDistanceHashMap::new();
    let mut edit_distance_counts = KAndEditDistanceHashMap::new();

    let file = File::open(sequence_database_file_name)?;
    let mut rdr = ReaderBuilder::new().from_reader(file);
    for result in rdr.deserialize() {
        let record: DatabaseRecord = result?; // serde deserializes this for us!!!
        for k in k_range {
            let estimated_distance = tensor::l2norm(&record.base_sequence.chars().collect(), &record.modified_sequence.chars().collect(), k);
            edit_distance_sums.update(&k, &record.edit_distance, estimated_distance);
            edit_distance_squared_sums.update(&k, &record.edit_distance, estimated_distance * estimated_distance);
            edit_distance_counts.update(&k, &record.edit_distance, 1 as f64);
        }
    }

    // Compute mean and confidence intervals from rolling data.
    let mut lower_bounds = KAndEditDistanceHashMap::new();
    let mut upper_bounds = KAndEditDistanceHashMap::new();
    let mut means = KAndEditDistanceHashMap::new();
    for k in k_range {
        for edit_distance in min_edit_distance..max_edit_distance {
            let sum = edit_distance_sums.get(&k, &edit_distance);
            let squared_sum = edit_distance_squared_sums.get(&k, &edit_distance);
            let count = edit_distance_counts.get(&k, &edit_distance);

            let mean = sum / count;
            let variance = (squared_sum - mean * sum) / count;
            let mean_se = (variance / count).sqrt();

            means.update(&k, &edit_distance, mean);
            lower_bounds.update(&k, &edit_distance, mean - mean_se);
            upper_bounds.update(&k, &edit_distance, mean + mean_se);
        }
    }

    // Create the header row for comparison data file
    let mut header_row: Vec<String> = Vec::new();
    header_row.push("edit distance".to_string());
    for k in k_range {
        header_row.push(format!("mean (k={})", k));
        header_row.push(format!("lower confidence bound (k={})", k));
        header_row.push(format!("upper confidence bound (k={})", k));
    }
    let header_row = header_row.join(",");

    let mut file = File::create(comparison_data_file_name)?;
    writeln!(&mut file, "{}", header_row)?;

    for edit_distance in min_edit_distance..max_edit_distance {
        let mut row_values: Vec<String> = Vec::new();
        row_values.push(edit_distance.to_string());
        for k in k_range {
            row_values.push(format!("{}", means.get(&k, &edit_distance))); // mean
            row_values.push(format!("{}", lower_bounds.get(&k, &edit_distance)));
            row_values.push(format!("{}", upper_bounds.get(&k, &edit_distance)));
        }
        let row = row_values.join(",");
        writeln!(&mut file, "{}", row)?;
    }
    Ok(())
}