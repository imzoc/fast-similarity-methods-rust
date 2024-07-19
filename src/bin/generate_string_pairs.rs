use std::fs::File;
use std::io::{Write, Result};

use alignment_free_methods::utils::sequence::generate_pair;

/*
 * This script generates a bunch of string pairs with varying pre-computed
 * Levenshtein distances, from 0 to 0.4 * string length.
 */
fn main() -> Result<()> {
    let n_sequences = 5;
    let sequence_length = 1000;
    let starting_edit_distance = 0;
    let max_edit_distance = (0.4 * sequence_length as f64) as usize;
    let edit_distance_range = starting_edit_distance..max_edit_distance;

    let mut file = File::create(format!("tests/inputs/sequences_{:?}.csv", sequence_length))?;
    writeln!(&mut file, "base_sequence,modified_sequence,edit_distance")?;

    for _ in 0..n_sequences {
        for edit_distance in edit_distance_range.clone() {
            let (base_sequence_chars, modified_sequence_chars) = generate_pair(sequence_length, edit_distance);

            let base_string: String = base_sequence_chars.iter().collect();
            let modified_string: String = modified_sequence_chars.iter().collect();
            writeln!(&mut file, "{:?},{:?},{}", base_string, modified_string, edit_distance)?;
        }
    }

    Ok(())
}