use anyhow::{anyhow, bail, Result};
use seahash::SeaHasher;
use std::hash::{Hash, Hasher};
use std::collections::HashMap;

#[derive(Debug)]
pub struct Tensor {
    pub data: Vec<usize>,
    pub shape: Vec<usize>,
}

#[allow(dead_code)]
impl Tensor {
    pub fn new_empty(k: usize) -> Tensor {
        let shape = vec![4; k];
        let size = shape.iter().product();
        let data = vec![0; size];

        Tensor { data, shape }
    }

    pub fn construct_with_step(sequence: &[char], k: usize, step: usize) -> Result<Self> {
        let mut tensor = Tensor::new_empty(k);
        let kmers = generate_kmers_with_step(sequence, k, step)?;
        tensor.populate(&kmers)?;
        Ok(tensor)
    }

    pub fn construct(sequence: &[char], k: usize) -> Result<Self> {
        Tensor::construct_with_step(sequence, k, 1)
    }

    pub fn shape(&self) -> &Vec<usize> {
        &self.shape
    }

    pub fn data(&mut self) -> &mut Vec<usize> {
        &mut self.data
    }

    pub fn data_ref(&self) -> &Vec<usize> {
        &self.data
    }

    pub fn to_vec(&self) -> Vec<usize> {
        self.data.clone()
    }

    pub fn get(&self, indices: &[usize]) -> Option<&usize> {
        let index = self.calculate_index(indices)?;
        self.data.get(index)
    }

    pub fn get_mut(&mut self, indices: &[usize]) -> Option<&mut usize> {
        let index = self.calculate_index(indices)?;
        self.data.get_mut(index)
    }

    fn calculate_index(&self, indices: &[usize]) -> Option<usize> {
        if indices.len() != self.shape.len() {
            return None;
        }
        let mut index = 0;
        let mut stride = 1;
        for (i, &dim_index) in indices.iter().rev().enumerate() {
            if dim_index >= self.shape[self.shape.len() - 1 - i] {
                return None;
            }
            index += dim_index * stride;
            stride *= self.shape[self.shape.len() - 1 - i];
        }
        Some(index)
    }

    pub fn populate(&mut self,  kmers: &Vec<&[char]>) -> Result<()> {
        let dimensions = self.shape.len();
        //const ALPHABET: [char; 4] = ['A', 'C', 'T', 'G']; // unused

        for kmer in kmers {
            if kmer.len() != dimensions {
                bail!(
                    "K ({}) does not match tensor dimensions ({})",
                    kmer.len(),
                    dimensions
                );
            }

            // Convert k-mer to multi-dimensional index
            let indices: Vec<usize> = kmer
                .iter()
                .map(|&c| match c {
                    'A' => 0,
                    'C' => 1,
                    'T' => 2,
                    'G' => 3,
                    _ => 666,
                })
                .collect();

            // Calculate the position in the tensor
            let index = self.calculate_index(&indices).ok_or_else(|| {
                anyhow!("Failed to calculate index for k-mer: {kmer:?}")
            })?;

            // Increment the value at the calculated position
            self.data[index] += 1;
        }

        Ok(())
    }
}

// This function simply takes a string and returns the set of k-mers.
pub fn generate_kmers(sequence: &[char], k: usize) -> Result<Vec<&[char]>> {
    generate_kmers_with_step(sequence, k, 1)
}

pub fn generate_kmers_with_step(
    sequence: &[char],
    k: usize,
    step: usize
) -> Result<Vec<&[char]>> {
    let mut kmers = vec![];

    if k > sequence.len() {
        bail!(
            "K ({}) is less than string length {}",
            k,
            sequence.len()
        );
    }

    for i in (0..=(sequence.len() - k)).step_by(step) {
        kmers.push(&sequence[i..i + k]);
    }
    Ok(kmers)
}

/* This function calculates the euclidean distance between kmer-vector representations of strings.
 * This is calculated as the l^2 norm of one vector subtracted from the other.
 */
pub fn kmer_euclidean_distance(
    base_seq: &[char],
    mod_sequence: &[char],
    k: usize,
    step: usize
) -> Result<f64> {
    let base_kmers = generate_kmers_with_step(base_seq, k, step)?;
    let mod_kmers = generate_kmers_with_step(mod_sequence, k, step)?;

    kmer_euclidean_distance_from_kmers(base_kmers, mod_kmers, k)
}

pub fn kmer_euclidean_distance_from_kmers(
    base_kmers: Vec<&[char]>,
    mod_kmers: Vec<&[char]>,
    k: usize
) -> Result<f64> {
    // At the moment, Tensor cannot handle unified kmers (it would require an extra dimension)
    // I can modify Tensor. However, I think I can do this computation without a tensor.
    // That is, I think I can construct the kmer occurrence vector without constructing the tensor first.
    let mut base_tensor = Tensor::new_empty(k);
    let mut mod_tensor = Tensor::new_empty(k);
    base_tensor.populate(&base_kmers)?;
    mod_tensor.populate(&mod_kmers)?;

    let base_kmer_vec = base_tensor.to_vec();
    let mod_kmer_vec = mod_tensor.to_vec();

    let mut distance = 0;
    for (x1, x2) in base_kmer_vec
        .iter()
        .zip(mod_kmer_vec.iter())
    {
        distance += (x1.max(x2) - x1.min(x2)).pow(2);
    }
    Ok((distance as f64).sqrt())
}

/* This function calculates the cosine similarity between two strings' kmeer vectors.
 * cosine similarity = (v1 * v2) / (||v1|| * ||v2||)
 */
pub fn kmer_cosine_similarity(
    base_seq: &[char],
    mod_seq: &[char],
    k: usize,
    step: usize
) -> Result<f64> {
    let base_seq_kmer_vector =
        Tensor::construct_with_step(&base_seq, k, step)?
            .to_vec()
            .into_iter()
            .map(|v| v as f64);
    let mod_seq_kmer_vector =
        Tensor::construct_with_step(&mod_seq, k, step)?
            .to_vec()
            .into_iter()
            .map(|v| v as f64);

    let mut base_seq_magnitude = 0.0;
    let mut mod_seq_magnitude = 0.0;
    let mut dot_product = 0.0;
    for (base_element, mod_element) in
        base_seq_kmer_vector.zip(mod_seq_kmer_vector)
    {
        base_seq_magnitude += base_element * base_element;
        mod_seq_magnitude += mod_element * mod_element;
        dot_product += base_element * mod_element;
    }
    if base_seq_magnitude == 0.0 || mod_seq_magnitude == 0.0 {
        return Ok(0.0);
    }
    Ok(dot_product / (base_seq_magnitude.sqrt() * mod_seq_magnitude.sqrt()))
}

/* This function calculates the euclidean distance between kmer-vector representations
* of minimizer representations of strings. */
pub fn minimizer_euclidean_distance(
    base_seq: &[char],
    mod_seq: &[char],
    k: usize,
    w: usize,
    step: usize
) -> Result<f64> {
    let base_minimizers = generate_minimizers(base_seq, k, w, step)?;
    let mod_minimizers = generate_minimizers(mod_seq, k, w, step)?;

    kmer_euclidean_distance_from_kmers(base_minimizers, mod_minimizers, k)
}

/* This function generates a set of minimizers from a string. It uses my_hash_function
 * with a set seed (69). It only uses ONE hash function.
 */
pub fn generate_minimizers(
    seq: &[char],
    k: usize,
    w: usize,
    step: usize
) -> Result<Vec<&[char]>> {
    if k > w {bail!("k ({}) > w ({})", k, w);}
    if w > seq.len() {bail!("w ({}) > sequence length ({})", w, seq.len());}

    let mut minimizers: Vec<&[char]> = Vec::new();
    for idx in (0..seq.len() - w).step_by(step) {
        let window = &seq[idx..idx + w];
        let minimizer = generate_single_minimizer(window, k)?;
        minimizers.push(minimizer);
    }
    Ok(minimizers)
}

pub fn generate_single_minimizer(window: &[char], k: usize) -> Result<&[char]> {
    let seed: u64 = 69;
    let kmers = generate_kmers(window, k)?;
    let minimizer = kmers
            .iter()
            .min_by_key(|item| my_hash_function(item, seed))
            .unwrap();
    Ok(minimizer)
}

/* This is a simple hash function that maps items (kmers) to u64 integers.
* It uses a seed, so this is technically a series of hash functions.
*/
fn my_hash_function<T: Hash> (item: &T, seed: u64) -> u64 {
    let mut hasher = SeaHasher::with_seeds(seed, seed, seed, seed);
    item.hash(&mut hasher);
    hasher.finish()
}

/* This struct manages ownership of string slices. Strobemers are constructed from strobes,
 * each of which is its own borrowed string slice. These slices need to be concatenated
 * into a bigger string slice, so this struct takes ownership so the concatenated string
 * slice stays alive while it's being used for calculations.
 */
struct Strobemer {
    chars: Vec<char>
}

impl Strobemer {
    pub fn new(strobes: Vec<&[char]>) -> Self {
        let strobemer_length: usize = strobes.iter()
            .map(|strobe| strobe.len()).sum();
        let mut strobemer = Vec::with_capacity(strobemer_length);
        for strobe in strobes {
            strobemer.extend_from_slice(strobe);
        }
        Strobemer {chars: strobemer}
    }

    pub fn data_ref(&self) -> &[char] {
        &self.chars
    }
}

/*  This function generates a string's set of strobemers. Strobemers are quirky,
so here's an intuitive breakdown of what strobemers are and how they're generated:

* Stromebers are the concatenation of strobes. There are *order* strobes in each
strobemer.
* Each strobemer is drawn from a span of length strobemer_span_length =
            strobe_length + (strobe_window_gap + strobe_window_length) * order,
representing the first strobe and all of the remaining strobe windows.
* In any sequence, there are
            (sequence.len() - strobemer_span_length) / step
strobemers. 
* The n'th strobemer is drawn from
            sequence[n * step..n * step + order * l].
* In each strobemer there are *order* strobes.
* In each strobemer, the first strobe's index is alwoys at the very start of the
strobemer's span. It is NOT determined by the argmin of a hash function.
* The n'th strobe (not including the first strobe, of course) is the l-mer minimizer
of a window within the strobemer's span. Specifically, that window is:
        window_end = strobe_length + (strobe_window_gap + strobe_window_length) * n
        window_end = window_start - strobe_window_length
        strobemer_span[window_start..window_end].
*/
pub fn strobemer_euclidean_distance(
    base_seq: &[char],
    mod_seq: &[char],
    order: usize,                   // the number of concatenated strobes in each strobemer.
    strobe_length: usize,           // the length of each strobe.
    strobe_window_gap: usize,       // the gap between strobe selection windows.
    strobe_window_length: usize,    // the size of each strobe's selection window.
    step: usize                     // the space between strobemer selection windows.
)-> Result<f64> {
    if strobe_window_length < strobe_length {
        bail!("Strobe window length ({}) is smaller than the length of the strobe ({})",
            strobe_window_length, strobe_length);
    }

    let base_strobemers = generate_strobemers(
        base_seq,
        order,
        strobe_length,
        strobe_window_gap,
        strobe_window_length,
        step
    )?;
    let base_strobemers_chars: Vec<&[char]> = base_strobemers
        .iter().map(|strobemer| strobemer.data_ref()).collect();
    let mod_strobemers = generate_strobemers(
        mod_seq,
        order,
        strobe_length,
        strobe_window_gap,
        strobe_window_length,
        step
    )?;
    let mod_strobemers_chars: Vec<&[char]> = mod_strobemers
        .iter().map(|strobemer| strobemer.data_ref()).collect();

    kmer_euclidean_distance_from_kmers(
        base_strobemers_chars,
        mod_strobemers_chars,
        order * strobe_length)
}

/* HELPER FUNCTION FOR strobemer_euclidean_distance()
 * This function generates a string's set of strobemers.
 */
fn generate_strobemers(
    seq: &[char],
    order: usize,
    strobe_length: usize,
    strobe_window_gap: usize,
    strobe_window_length: usize,
    step: usize
) -> Result<Vec<Strobemer>> {
    let mut strobemers: Vec<Strobemer> = Vec::new(); // custom struct to manage slice ownership.
    let strobemer_span = strobe_length +
        order * (strobe_window_length + strobe_window_gap);
    let last_strobemer_start_index = seq.len() - strobemer_span; // try + 1?

    for idx in (0..last_strobemer_start_index).step_by(step) {
        let strobemer_window = &seq[idx..idx + strobemer_span];
        let strobemer = generate_single_strobemer(
            strobemer_window,
            order,
            strobe_length,
            strobe_window_gap,
            strobe_window_length
        )?;
        strobemers.push(strobemer);
    }
    Ok(strobemers)
}

/* This function generates the first strobemer from a designated window. This means the first
 * strobe is always seq[0..=l]. The subsequent strobes are generated by generate_single_minimizer()
 * from the window designated by the strobemer paper.
 *
 * To generate a string's entire set of strobemers, use this function as a HELPER FUNCTION.
 * Pass it a reference to the substring you want to generate a strobemer from, with the first
 * strobemer being seq[0..=l].
 */
fn generate_single_strobemer(
    strobemer_window: &[char],
    order: usize,
    strobe_length: usize,
    strobe_window_gap: usize,
    strobe_window_length: usize
) -> Result<Strobemer> {
    let mut strobes: Vec<&[char]> = Vec::new();
    let first_strobe = &strobemer_window[0..strobe_length];
    strobes.push(first_strobe);
    for n in 1..order {
        let window_end = strobe_length + (strobe_window_gap + strobe_window_length) * n;
        let window_start = window_end - strobe_window_length;
        let strobe_window = &strobemer_window[window_start..window_end];
        let strobe = generate_single_minimizer(strobe_window, strobe_length)?;
        strobes.push(strobe);
    }
    Ok(Strobemer::new(strobes))
}

/*  */
pub fn unified_minimizer_euclidean_distance(
    base_seq: &[char],
    mod_seq: &[char],
    k: usize,
    w: usize
) -> Result<f64> {
    let base_minimizers = generate_unified_minimizers(&base_seq, w, k)?;
    let mod_minimizers = generate_unified_minimizers(&mod_seq, w, k)?;
    kmer_euclidean_distance_from_unified_kmers(base_minimizers, mod_minimizers, k, w)
}

/* This function generates a set of minimizers from a string. It uses my_hash_function
 * with a set seed (69). It only uses ONE hash function.
 */
pub fn generate_unified_minimizers(
    seq: &[char],
    w: usize,
    k: usize
) -> Result<HashMap<&[char], usize>> {
    let seed: u64 = 69;
    /* This is a simple hash function that maps items (kmers) to u64 integers.
    * It uses a seed, so this is technically a series of hash functions.
    */
    fn my_unified_hash_function<T: Hash> (item: &T, occurrence_number: usize, seed: u64) -> u64 {
        let mut hasher = SeaHasher::with_seeds(seed, seed, seed, seed);
        (item, occurrence_number).hash(&mut hasher);
        hasher.finish()
    }

    if k > w {bail!("k ({}) > w ({})", k, w);}
    if w > seq.len() {bail!("w ({}) > sequence length ({})", w, seq.len());}

    let mut unified_minimizers: HashMap<&[char], usize> = HashMap::new();
    for idx in 0..seq.len() - w {
        let window = &seq[idx..idx + w];
        let kmers: Vec<&[char]> = generate_kmers(window, k)?;
        let unified_kmers: Vec<(&[char], usize)> = kmers
            .iter()
            .cloned()
            .zip(kmers
                .iter()
                .map(|&kmer| *unified_minimizers.get(&kmer).unwrap_or(&0_usize)))
            .collect();
        let minimizer = unified_kmers
                .iter()
                .min_by_key(|(item, o_n)| my_unified_hash_function(item, *o_n, seed))
                .unwrap();
        match unified_minimizers.get_mut(minimizer.0) {
            Some(occurrence_number) => *occurrence_number += 1,
            _ => {unified_minimizers.insert(minimizer.0, minimizer.1);}
        }
    }
    Ok(unified_minimizers)
}

pub fn kmer_euclidean_distance_from_unified_kmers(
    base_minimizers: HashMap<&[char], usize>,
    mod_minimizers: HashMap<&[char], usize>,
    k: usize,
    w: usize
) -> Result<f64> {
    Ok(0.0)
}

/*  */
pub fn ordered_minimizer_euclidean_distance(
    base_seq: &[char],
    mod_seq: &[char],
    k: usize,
    w: usize
) -> Result<f64> {
    unimplemented!();
}

/* This is the final goal of my project! */
pub fn tensor_slide_sketch(base_seq: &[char],
    mod_seq: &[char],
    k: usize,
    w: usize
) -> Result<f64> {
    unimplemented!();
}

#[cfg(test)]
mod similarity_method_tests {
    use pretty_assertions::assert_eq;

    use crate::utils::tensor::{
        kmer_euclidean_distance,
        kmer_cosine_similarity,
        minimizer_euclidean_distance,
        strobemer_euclidean_distance
    };

    #[test]
    fn kmer_euclidean_distance_test_1() {
        let k = 4;
        let step = k / 2;

        let seq1: Vec<char> = "ACGTTGCA".chars().collect();
        let seq2 = seq1.clone();
        let expected = 0.0;

        let res = kmer_euclidean_distance(&seq1, &seq2, k, step);
        assert!(res.is_ok());
        assert_eq!(res.unwrap(), expected);
    }

    #[test]
    fn  kmer_cosine_similarity_test_1() {
        let k = 4;
        let step = k / 2;

        let seq1: Vec<char> = "ACGTTGCA".chars().collect();
        let seq2 = seq1.clone();
        let expected = 0.0;
        let res =  kmer_cosine_similarity(&seq1, &seq2, k, step);
        assert!(res.is_ok());
        assert_eq!(res.unwrap(), expected);
    }

    #[test]
    fn minimizer_euclidean_distance_test_1() {
        let k = 4;
        let w = 6;
        let step = w / 2;

        let seq1: Vec<char> = "ACGTTGCA".chars().collect();
        let seq2 = seq1.clone();
        let expected = 0.0;
        let res =  minimizer_euclidean_distance(&seq1, &seq2, k, w, step);
        assert!(res.is_ok());
        assert_eq!(res.unwrap(), expected);
    }
    #[test]
    fn strobemer_euclidean_distance_test_1() {
        let order = 2;
        let l = 3;
        let w_min = 6;
        let w_max  = 12;
        let step = (w_max - w_min) / 2;

        let seq1: Vec<char> = "ACGTTGCA".chars().collect();
        let seq2 = seq1.clone();
        let expected = 0.0;
        let res = strobemer_euclidean_distance(&seq1, &seq2, order, l, w_min, w_max, step);
        assert!(res.is_ok());
        assert_eq!(res.unwrap(), expected);
    }
}


#[cfg(test)]
mod tensor_struct_tests {
    use super::Tensor;
    use pretty_assertions::assert_eq;

    #[test]
    fn tensor_new_empty() {
        // Should 0 be allowed?
        let t = Tensor::new_empty(0);
        assert_eq!(t.shape.len(), 0);
        assert_eq!(t.shape, &[]);
        assert_eq!(t.data.len(), 1);
        assert_eq!(t.data, vec![0]);

        let t = Tensor::new_empty(2);
        assert_eq!(t.shape.len(), 2);
        assert_eq!(t.shape, [4; 2]);
        assert_eq!(t.data.len(), 16);
        assert_eq!(t.data, [0; 16]);
    }

    #[test]
    fn tensor_construct() {
        // Use a less trivial sequence?
        let res = Tensor::construct(&['A', 'C', 'G', 'T'], 2);
        assert!(res.is_ok());

        let t = res.unwrap();
        assert_eq!(t.data.len(), 16);
        assert_eq!(t.data, [0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,]);
        assert_eq!(t.shape.len(), 2);
        assert_eq!(t.shape, [4; 2]);
    }
}



#[cfg(test)]
mod standalone_function_tests {
    use super::generate_kmers;
    use pretty_assertions::assert_eq;

    #[test]
    fn test_generate_kmers_k1() {
        let seq: Vec<char> = "ACGT"
            .chars().collect();
        let kmers = generate_kmers(&seq, 1)
            .expect("Failed at k-mer generation");
        assert_eq!(kmers.len(), 4);
        for (i, c) in "ACGT".chars().enumerate() {
            assert_eq!(kmers[i], [c]);
        }

        let kmers = generate_kmers(&['A', 'C', 'G', 'T'], 2)
            .expect("Failed at k-mer generation");
        assert_eq!(kmers.len(), 3);
        assert_eq!(kmers[0], ['A', 'C']);
        assert_eq!(kmers[1], ['C', 'G']);
        assert_eq!(kmers[2], ['G', 'T']);
    }
}