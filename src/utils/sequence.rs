use rand::Rng;
use edit_distance::edit_distance;

const ALPHABET: [char; 4] = ['A', 'C', 'T', 'G'];

/*
 * This function takes an integer n and returns an owned, n-long
 * random string over the DNA alphabet.
 */
pub fn generate_sequence_chars(n: usize) -> Vec<char> {
    let mut rng = rand::thread_rng();

    // Generate n random indices within the range of DNA_ALPHABET
    let sequence_chars: Vec<char> = (0..n).map(|_| {
        ALPHABET[rng.gen_range(0..ALPHABET.len())]
    }).collect();
    sequence_chars
}

pub fn generate_pair(base_length: usize, distance: usize) -> (Vec<char>, Vec<char>) {
    let base_sequence_chars: Vec<char> = generate_sequence_chars(base_length);
    let mut modified_sequence_chars: Vec<char> = base_sequence_chars.clone();

    let base_string: String = base_sequence_chars.clone().into_iter().collect();
    let modified_string: String = modified_sequence_chars.clone().into_iter().collect();
    let mut distance_needed = distance as isize - edit_distance(&base_string, &modified_string) as isize;
    while distance_needed > 0 {
        for _ in 0..distance_needed {
            modify(&mut modified_sequence_chars);
        }

        let base_string: String = base_sequence_chars.clone().into_iter().collect();
        let modified_string: String = modified_sequence_chars.clone().into_iter().collect();
        distance_needed = distance as isize - edit_distance(&base_string, &modified_string) as isize;
    }

    (base_sequence_chars, modified_sequence_chars)
}

pub fn modify_string(s: &mut String) {
    let mut chars: Vec<char> = s.chars().collect();
    modify(&mut chars);
    *s = chars.iter().collect();
}

/* This function randomly modifies one index in a Vec<char>.
*/
pub fn modify(chars: &mut Vec<char>) {
    if chars.is_empty() {
        return;
    }

    let mut rng = rand::thread_rng();

    let index = rng.gen_range(0..chars.len());
    let new_character = ALPHABET[rng.gen_range(0..ALPHABET.len())];
    let which_operation = rng.gen_range(0..3);
    match which_operation {
        0 => { // Insert
            chars.insert(index, new_character);
        }
        1 => { // Delete
            chars.remove(index);
        }
        2 => { // Change the letter
            chars[index] = new_character
        }
        _ => unreachable!(), // This should never happen with range 0..3
    }
}