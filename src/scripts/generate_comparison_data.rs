use std::fs::File;
use std::io::{self, Result, Write,BufReader,BufRead};
use std::collections::HashMap;

use csv::ReaderBuilder;

use similarity_methods::utils::tensor;
use similarity_methods::utils::sequence;

fn read_metadata(sequence_database_file_name: &String) -> Result<(usize, (usize, usize)), Box<dyn Error>> {
    let file = File::open(sequence_database_file_name)?;
    let mut reader = BufReader::new(file);

    let mut first_line = String::new();
    reader.read_line(&mut first_line)?;
    let sequence_length: usize = first_line.trim().strip_prefix("# base string length=")
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Invalid metadata format"))?
        .parse().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let mut second_line = String::new();
    reader.read_line(&mut second_line)?;
    let min_edit_distance: usize = second_line.trim().strip_prefix("# minimum edit distance=")
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Invalid metadata format"))?
        .parse().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let mut third_line = String::new();
    reader.read_line(&mut third_line)?;
    let max_edit_distance: usize = third_line.trim().strip_prefix("# maximum edit distance=")
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Invalid metadata format"))?
        .parse().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok((sequence_length, (min_edit_distance, max_edit_distance)))
}

fn generate_comparison_data(sequence_database_file_name: &String) -> Result<(HashMap<usize, Vec<f64>>, HashMap<usize, Vec<f64>>, Vec<usize>), Box<dyn Error>> {
    // Initialize the data structures this function will return. 
    // edit_distance_sums and edit_distance_squared_sums are HashMaps that contain
    // the sums and squaned sums of estimated distances between strings for a given k and edit distance.
    let mut edit_distance_sums: HashMap<usize, Vec<f64>> = HashMap::new();
    let mut edit_distance_squared_sums: HashMap<usize, Vec<f64>> = HashMap::new();
    for k in K_RANGE {
        edit_distance_sums.insert(k, vec![0.0; max_edit_distance]);
        edit_distance_squared_sums.insert(k, vec![0.0; max_edit_distance]);
    }
    // edit_distance_counts contains the number of string pairs with a true edit distance given by the index.
    let mut edit_distance_counts: Vec<usize> = vec![0; max_edit_distance];

    let file = File::open(sequence_database_file_name)?;
    let mut rdr = ReaderBuilder::new().from_reader(file);

    for result in rdr.records() {
        let record = result?;
        let string1_chars: Vec<char> = record[0].chars().collect(); // (string of ~1000 letters)
        let string2_chars: Vec<char> = record[1].chars().collect(); // (string of ~1000 letters)
        let edit_distance: usize = record[2].parse().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        edit_distance_counts[edit_distance] += 1;
        
        // COMPUTATIONS
        for k in K_RANGE {
            let estimated_distance = tensor::l2norm(&string1_chars, &string2_chars, k);
            edit_distance_sums.get_mut(&k).unwrap()[edit_distance] += estimated_distance;
            edit_distance_squared_sums.get_mut(&k).unwrap()[edit_distance] += estimated_distance * estimated_distance;
        }
    }
    Ok((edit_distance_sums, edit_distance_squared_sums, edit_distance_counts))
}

fn main() -> Result<()> {
    const K_RANGE: [usize; 3] = [2,3,4];
    let sequence_database_file_name = "sequences_1000.csv".to_string();

    let (sequence_length, (min_edit_distance, max_edit_distance)) = match read_metadata(&sequence_database_file_name) {
        Ok(a, (b, c)) => (a, (b, c)),
        Err(e) => {
            println!("error reading sequence length: {}", err);
            process::exit(1);
        }
    }

    let (edit_distance_sums: HashMap<usize, Vec<f64>>, edit_distance_squared_sums: HashMap<usize, Vec<f64>>, edit_distance_counts: Vec<usize>) = match generate_comparison_data(&sequence_database_file_name, &K_RANGE) {
        Ok(a, b, c) => (a, b, c),
        Err(e) => {
            println!("error generating estimated distance data from the file: {}", err);
            process::exit(1);
        }
    }

    let mut confidence_intervals: HashMap<usize, Vec<(f64, f64)>> = HashMap::new();
    let mut means: HashMap<usize, Vec<f64>> = HashMap::new();
    for k in K_RANGE {
        confidence_intervals.insert(k, vec![(0.0, 0.0); max_edit_distance]);

        means.insert(k, vec![0.0; max_edit_distance]);
        for edit_distance in min_edit_distance..max_edit_distance {
            confidence_intervals.get_mut(&k).unwrap()[edit_distance] = {
                let sum = edit_distance_sums[&k][edit_distance];
                let squared_sum = edit_distance_squared_sums[&k][edit_distance];
                let count = edit_distance_count[edit_distance] as f64;

                let mean = sum / count;
                let variance = (squared_sum - mean * sum) / count;
                let mean_SE = (variance / count as f64).sqrt();
                means.get_mut(&k).unwrap()[edit_distance] = mean;
                (mean - mean_SE, mean + mean_SE)
            }
        }
    }

    let mut file = File::create("data.csv")?;

    let mut column_names: Vec<String> = Vec::new();
    column_names.push("edit distance".to_string());
    for k in K_RANGE {
        column_names.push(format!("mean (k={})", k));
        column_names.push(format!("lower bound (k={})", k));
        column_names.push(format!("upper bound (k={})", k));
    }

    let header_row = column_names.join(",");
    writeln!(&mut file, "{}", header_row)?;

    for edit_distance in min_edit_distance..max_edit_distance {
        let mut row_values: Vec<String> = Vec::new();
        row_values.push(edit_distance.to_string());
        for k in K_RANGE {
            row_values.push(means[&k][edit_distance].to_string()); // mean
            row_values.push(format!("{}", confidence_intervals[&k][edit_distance].0.to_string()));
            row_values.push(format!("{}", confidence_intervals[&k][edit_distance].1.to_string()));
        }
        let row = row_values.join(",");
        writeln!(&mut file, "{}", row)?;
    }
    Ok(())
}