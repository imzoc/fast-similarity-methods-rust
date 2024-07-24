use std::fs::File;
use std::io::{self, Result, Write};
use bio::io::fasta;
use alignment_free_methods::utils::sequence::generate_pair;

/*
 * This script generates a bunch of string pairs with varying pre-computed
 * Levenshtein distances, from 0 to 0.4 * string length.
 */
fn main() -> Result<()> {
    let n_sequences = 5;
    let sequence_length = 1000;
    let starting_edit_distance = 0;
    let max_edit_distance = (0.6 * sequence_length as f64) as usize;
    let edit_distance_range = starting_edit_distance..max_edit_distance;

    let path = "tests/artificial.fasta";
    let file = File::create(path)?;
    let mut writer = fasta::Writer::new(file);

    let mut i = 0;
    for _ in 0..n_sequences {
        for edit_distance in edit_distance_range.clone() {
            let (base_sequence_chars, modified_sequence_chars) = generate_pair(sequence_length, edit_distance);

            let base_string: String = base_sequence_chars.iter().collect();
            let modified_string: String = modified_sequence_chars.iter().collect();
            i += 1;
            writer.write_record(&fasta::Record::with_attrs(
                &format!("FAKE_{}", i), 
                Some(&format!("Artificially generated seed. {} edits away from FAKE_{}", edit_distance, i + 1)),
                base_string.as_bytes())
            )?;
            i += 1;
            writer.write_record(&fasta::Record::with_attrs(
                &format!("FAKE_{}", i), 
                Some(&format!("Artificially generated seed. {} edits away from FAKE {}", edit_distance, i - 1)),
                modified_string.as_bytes())
            )?;
            print!("\rString pair #{:?} has been processed!", i);
            io::stdout().flush()?;
        }
    }

    Ok(())
}