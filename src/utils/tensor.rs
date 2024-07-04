use anyhow::{anyhow, bail, Result};
use seahash::SeaHasher;
use std::hash::{Hash, Hasher};

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

    pub fn construct(sequence: &[char], k: usize) -> Result<Self> {
        let mut tensor = Tensor::new_empty(k);
        let kmers = generate_kmers(sequence, k)?;
        tensor.populate(&kmers)?;
        Ok(tensor)
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

    pub fn populate(&mut self, kmers: &Vec<&[char]>) -> Result<()> {
        let dimensions = self.shape.len();
        //const ALPHABET: [char; 4] = ['A', 'C', 'T', 'G']; // unused

        for kmer in kmers {
            if kmer.len() != dimensions {
                bail!(
                    "K-mer length {} does not match tensor dimensions {}",
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

#[derive(Debug)]
pub struct SequenceBasedParameters {
    pub k: usize,
    pub w: usize,
    pub base_sequence: Vec<char>,
    pub modified_sequence: Vec<char>,
}

#[derive(Debug)]
pub struct KmerBasedParameters<'a> {
    pub k: usize,
    pub w: usize,
    pub base_kmers: Vec<&'a [char]>,
    pub modified_kmers: Vec<&'a [char]>,
}


// This function simply takes a string and returns the set of k-mers.
pub fn generate_kmers(sequence: &[char], k: usize) -> Result<Vec<&[char]>> {
    let mut kmers = vec![];

    if k > sequence.len() {
        bail!(
            "K ({}) is less than string length {}",
            k,
            sequence.len()
        );
    }

    for i in 0..=(sequence.len() - k) {
        kmers.push(&sequence[i..i + k]);
    }
    Ok(kmers)
}

/* This function calculates the euclidean distance between kmer-vector representations of strings.
 * This is calculated as the l^2 norm of one vector subtracted from the other.
 */
pub fn euclidean_distance(params: SequenceBasedParameters) -> Result<f64> {
    let base_kmers = generate_kmers(&params.base_sequence, params.k)?;
    let mod_kmers = generate_kmers(&params.modified_sequence, params.k)?;

    euclidean_distance_from_kmers(KmerBasedParameters {
        base_kmers: base_kmers, modified_kmers: mod_kmers,
        k: params.k, w: params.w
    })
}

/* This function does the same thing, but takes k-mer sets instead (this is important
 * for minimizer methods)
 */
pub fn euclidean_distance_from_kmers(params: KmerBasedParameters) -> Result<f64> {
    let mut base_tensor = Tensor::new_empty(params.k);
    let mut mod_tensor = Tensor::new_empty(params.k);
    base_tensor.populate(&params.base_kmers)?;
    mod_tensor.populate(&params.modified_kmers)?;

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
pub fn cosine_similarity(params: SequenceBasedParameters) -> Result<f64> {
    let base_sequence_kmer_vector =
        Tensor::construct(&params.base_sequence, params.k)?
            .to_vec()
            .into_iter()
            .map(|v| v as f64);
    let modified_sequence_kmer_vector =
        Tensor::construct(&params.modified_sequence, params.k)?
            .to_vec()
            .into_iter()
            .map(|v| v as f64);

    let mut base_sequence_magnitude = 0.0;
    let mut modified_sequence_magnitude = 0.0;
    let mut dot_product = 0.0;
    for (base_element, mod_element) in
        base_sequence_kmer_vector.zip(modified_sequence_kmer_vector)
    {
        base_sequence_magnitude += base_element * base_element;
        modified_sequence_magnitude += mod_element * mod_element;
        dot_product += base_element * mod_element;
    }
    if base_sequence_magnitude == 0.0 || modified_sequence_magnitude == 0.0 {
        return Ok(0.0);
    }
    Ok(dot_product / (base_sequence_magnitude.sqrt() * modified_sequence_magnitude.sqrt()))
}

/* This is a simple hash function that maps items (kmers) to u64 integers.
 * It uses a seed, so this is technically a series of hash functions.
 */
fn my_hash_function<T: Hash> (item: &T, seed: u64) -> u64 {
    let mut hasher = SeaHasher::with_seeds(seed, seed, seed, seed);
    item.hash(&mut hasher);
    hasher.finish()
}
/* This function generates a set of minimizers from a string. It uses my_hash_function
 * with a set seed (69). It only uses ONE hash function.
 */
pub fn generate_minimizers(seq: &Vec<char>, w: usize, k: usize) -> Result<Vec<&[char]>> {
    let seed: u64 = 69;

    if k > w {bail!("k ({}) > w ({})", k, w);}
    if w > seq.len() {bail!("w ({}) > sequence length ({})", w, seq.len());}

    let mut minimizers: Vec<&[char]> = Vec::new();
    for idx in 0..seq.len() - w {
        let base_window = &seq[idx..idx + w];
        let base_kmers = generate_kmers(base_window, k)?;
        let base_minimizer = base_kmers
                .iter()
                .min_by_key(|item| my_hash_function(item, seed))
                .unwrap();
        minimizers.push(base_minimizer);
    }
    Ok(minimizers)
}

/* This function calculates the euclidean distance between kmer-vector representations
* of minimizer representations of strings. */
pub fn minimizer_euclidean_distance(params: SequenceBasedParameters) -> Result<f64> {
    let base_minimizers = generate_minimizers(&params.base_sequence, params.w, params.k)?;
    let mod_minimizers = generate_minimizers(&params.modified_sequence, params.w, params.k)?;

    euclidean_distance_from_kmers(KmerBasedParameters {
        base_kmers: base_minimizers, modified_kmers: mod_minimizers,
        k: params.k, w: params.w
    })
}

/* I haven't implemented this yet and I also don't know how it works :) */
pub fn strobemer(_params: SequenceBasedParameters) -> Result<f64> {
    unimplemented!();
}

/* This is the final goal of my project! */
pub fn tensor_slide_sketch(_params: SequenceBasedParameters) -> Result<f64> {
    unimplemented!();
}

#[cfg(test)]
mod similarity_method_tests {
    use crate::utils::tensor::{cosine_similarity, minimizer_euclidean_distance, strobemer};

    use super::{SequenceBasedParameters, euclidean_distance};
    use pretty_assertions::assert_eq;

    // KMER VECTOR L2 NORM TESTS
    fn test_euclidean_distance(p: SequenceBasedParameters, expected: f64)  {
        let res = euclidean_distance(p);
        assert!(res.is_ok());
        assert_eq!(res.unwrap(), expected);
    }
    #[test]
    fn test_euclidean_distance_equal_k1() {
        // k=1, same number of each character, should result in euclidean_distance = 0.0
        let p = SequenceBasedParameters {
            k: 1,
            w: 1, // irrelevant for euclidean_distance
            base_sequence: "ACGT".chars().collect(),
            modified_sequence: "GCAT".chars().collect(),
        };
        test_euclidean_distance(p, 0.0);
    }
    #[test]
    fn test_euclidean_distance_equal_k2() {
        // k=2, same number of 2-mers, should result in euclidean_distance = 0.0
        let p = SequenceBasedParameters {
            k: 2,
            w: 1, // irrelevant for euclidean_distance
            base_sequence: "ACGTA".chars().collect(),
            modified_sequence: "CGTAC".chars().collect(),
        };
        test_euclidean_distance(p, 0.0);
    }
    #[test]
    fn test_euclidean_distance_unequal_k1() {
        // k=1, unequal character set
        let p = SequenceBasedParameters {
            k: 1,
            w: 1, // irrelevant for euclidean_distance
            base_sequence: "GTTTTT".chars().collect(),
            modified_sequence: "ATTTTT".chars().collect(),
        };
        test_euclidean_distance(p, 1.4142135623730951); // square root 2
    }
    #[test]
    fn test_euclidean_distance_unequal_k2() {
        // k=2, unequal 2-mer set, should result in euclidean_distance = sqrt(2)
        let p = SequenceBasedParameters {
            k: 1,
            w: 1, // irrelevant for euclidean_distance
            base_sequence: "GTTTTT".chars().collect(),
            modified_sequence: "ATTTTT".chars().collect(),
        };
        test_euclidean_distance(p, 1.4142135623730951);
    }

    // COSINE SIMILARITY TESTS
    fn test_cosine_similarity(p: SequenceBasedParameters, expected: f64) {
        let res = cosine_similarity(p);
        assert!(res.is_ok());
        assert_eq!(res.unwrap(), expected);
    }
    #[test]
    fn test_cosine_similarity_equal_k1() {
        // k=2, unequal 2-mer set, should result in euclidean_distance = sqrt(2)
        let p = SequenceBasedParameters {
            k: 1,
            w: 1, // irrelevant for euclidean_distance
            base_sequence: "ACTG".chars().collect(),
            modified_sequence: "CTAG".chars().collect(),
        };
        test_cosine_similarity(p, 1.0);
    }
    #[test]
    fn test_cosine_similarity_equal_k2() {
        // k=2, unequal 2-mer set, should result in euclidean_distance = sqrt(2)
        let p = SequenceBasedParameters {
            k: 2,
            w: 1, // irrelevant for euclidean_distance
            base_sequence: "ACTGA".chars().collect(),
            modified_sequence: "CTGAC".chars().collect(),
        };
        test_cosine_similarity(p, 1.0);
    }
    #[test]
    fn test_cosine_similarity_unequal_k2() {
        // k=2, unequal 2-mer set, should result in euclidean_distance = sqrt(2)
        let p = SequenceBasedParameters {
            k: 2,
            w: 1, // irrelevant for euclidean_distance
            base_sequence: "ACTGG".chars().collect(),
            modified_sequence: "CTGAC".chars().collect(),
        };
        test_cosine_similarity(p, 0.75);
    }
    #[test]
    fn test_cosine_similarity_unequal_k1() {
        // k=2, unequal 2-mer set, should result in euclidean_distance = sqrt(2)
        let p = SequenceBasedParameters {
            k: 1,
            w: 1, // irrelevant for euclidean_distance
            base_sequence: "GTT".chars().collect(),
            modified_sequence: "ATT".chars().collect(),
        };
        test_cosine_similarity(p, 0.7999999999999998);
    }

    // MINIMIZER L2 NORM TESTS
    fn test_minimizer_euclidean_distance(p: SequenceBasedParameters, expected: f64) {
        let res = minimizer_euclidean_distance(p);
        assert!(res.is_ok());
        assert_eq!(res.unwrap(), expected);
    }
    #[test]
    fn test_minimizer_euclidean_distance_equal_w1() {
        // Some less trivial example?
        let p = SequenceBasedParameters {
            k: 1,
            w: 1,
            base_sequence: "AAA".chars().collect(),
            modified_sequence: "AAA".chars().collect(),
        };
        test_minimizer_euclidean_distance(p, 0.0);
    }
    #[test]
    fn test_minimizer_euclidean_distance_equal_w3() {
        // Some less trivial example?
        let p = SequenceBasedParameters {
            k: 2,
            w: 3,
            base_sequence: "AAAAA".chars().collect(),
            modified_sequence: "AAAAA".chars().collect(),
        };
        test_minimizer_euclidean_distance(p, 0.0);
    }
    #[test]
    fn test_minimizer_euclidean_distance_unequal_k3_w_equals_l() {
        // Some less trivial example?
        let p = SequenceBasedParameters {
            k: 2,
            w: 2,
            base_sequence: "AAAAA".chars().collect(),
            modified_sequence: "AAACC".chars().collect(),
        };
        test_minimizer_euclidean_distance(p, 0.0);
    }
    #[test]
    fn test_minimizer_euclidean_distance_unequal_k2_w4() {
        let p = SequenceBasedParameters {
            k: 2,
            w: 4,
            base_sequence: "ACTC".chars().collect(),
            modified_sequence: "ACTG".chars().collect(),
        };
        test_minimizer_euclidean_distance(p, 0.0);
    }
    #[test]
    fn test_minimizer_euclidean_distance_unequal_k1_w4_long() {
        // I should study this one a bit more...
        let p = SequenceBasedParameters {
            k: 2,
            w: 4,
            base_sequence: "AAAATTTT".chars().collect(),
            modified_sequence: "AAAAGGGG".chars().collect(),
        };
        test_minimizer_euclidean_distance(p, 0.0);
    }

    fn test_strobemer(p: SequenceBasedParameters, expected: f64) {
        let res = strobemer(p);
        assert!(res.is_ok());
        assert_eq!(res.unwrap(), expected);
    }
    #[test]
    fn test_strobemer_equal_k1() {
        let p = SequenceBasedParameters {
            k: 1,
            w: 1,
            base_sequence: vec!['A', 'C', 'G', 'T'],
            modified_sequence: vec!['T', 'C', 'A', 'T'],
        };
        test_strobemer(p, 0.0);
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