use anyhow::{anyhow, bail, Result};
use seahash::SeaHasher;
use std::hash::{Hash, Hasher};

#[derive(Debug)]
pub struct Tensor {
    pub data: Vec<usize>,
    pub shape: Vec<usize>,
}

#[derive(Debug)]
pub struct Parameters {
    pub k: usize,
    pub w: usize,
    pub base_sequence: Vec<char>,
    pub modified_sequence: Vec<char>,
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
        let kmers = generate_kmers(sequence, k);
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

pub fn generate_kmers(sequence: &[char], k: usize) -> Vec<&[char]> {
    let mut kmers = vec![];

    if k > sequence.len() {
        return kmers;
    }

    for i in 0..=(sequence.len() - k) {
        kmers.push(&sequence[i..i + k]);
    }
    kmers
}

/* This similarity estimation method returns the l2norm between two strings'
 * kmer vectors. Recall, a string's kmer vector is a vector where every
 * dimension represents the number of times a certain kmer appears in that string.
 * 
 * Thus, a low l2 norm predicts that two strings are similar, and vice versa.
 */
pub fn l2norm(params: Parameters) -> Result<f64> {
    let base_sequence_kmer_vector =
        Tensor::construct(&params.base_sequence, params.k)?.to_vec();
    let modified_sequence_kmer_vector =
        Tensor::construct(&params.modified_sequence, params.k)?.to_vec();

    let mut distance = 0;
    for (x1, x2) in base_sequence_kmer_vector
        .iter()
        .zip(modified_sequence_kmer_vector.iter())
    {
        distance += (x1.max(x2) - x1.min(x2)).pow(2);
    }
    Ok((distance as f64).sqrt())
}

pub fn cosine_similarity(params: Parameters) -> Result<f64> {
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
    for (base_char, mod_char) in
        base_sequence_kmer_vector.zip(modified_sequence_kmer_vector)
    {
        //base_sequence_magnitude +=
        //    base_char.clone() as f64 * base_char.clone() as f64;
        //modified_sequence_magnitude +=
        //    mod_char.clone() as f64 * mod_char.clone() as f64;
        //dot_product += base_char.clone() as f64 * mod_char.clone() as f64;
        base_sequence_magnitude += base_char * base_char;
        modified_sequence_magnitude += mod_char * mod_char;
        dot_product += base_char * mod_char;
    }
    if base_sequence_magnitude == 0.0 || modified_sequence_magnitude == 0.0 {
        return Ok(0.0);
    }
    Ok(dot_product / (base_sequence_magnitude * modified_sequence_magnitude))
}


fn my_hash_function<T: Hash> (item: &T, seed: u64) -> u64 {
    let mut hasher = SeaHasher::with_seeds(seed, seed, seed, seed);
    item.hash(&mut hasher);
    hasher.finish()
}

pub fn minimizer_l2_norm(params: Parameters) -> Result<f64> {
    let mut base_minimizers: Vec<&[char]> = Vec::new();
    let mut mod_minimizers: Vec<&[char]> = Vec::new();

    let smallest_sequence_length = std::cmp::min(
        params.base_sequence.len(),
        params.modified_sequence.len(),
    );

    for idx in 0..smallest_sequence_length-params.w {
        let base_window = &params.base_sequence[idx..idx + params.w];
        let mod_window = &params.modified_sequence[idx..idx + params.w];

        let base_kmers = generate_kmers(base_window, params.k);
        let mod_kmers = generate_kmers(mod_window, params.k);

        let seed: u64 = 69;
        let base_minimizer = base_kmers
                .iter()
                .min_by_key(|item| my_hash_function(item, seed))
                .ok_or_else(||
                    anyhow!("Window={:?}; k={:?}; idx={:?}; smallest sequence length={:?}", base_window, &params.k, idx, smallest_sequence_length)
                );
        let mod_minimizer = mod_kmers
                .iter()
                .min_by_key(|item| my_hash_function(item, seed))
                .ok_or_else(||
                    anyhow!("Mod window: {:?}; mod {:?}-mers: {:?}", mod_window, &params.k, mod_kmers)
                );

        base_minimizers.push(base_minimizer?);
        mod_minimizers.push(mod_minimizer?);
    }

    let mut base_tensor = Tensor::new_empty(params.k);
    base_tensor.populate(&base_minimizers)?;
    let mut mod_tensor = Tensor::new_empty(params.k);
    mod_tensor.populate(&mod_minimizers)?;

    let mut distance = 0;
    for (x1, x2) in base_tensor
        .data_ref()
        .iter()
        .zip(mod_tensor.data_ref().iter())
    {
        distance += (x1.max(x2) - x1.min(x2)).pow(2);
    }

    Ok((distance as f64).sqrt())
}

pub fn strobemer(_params: Parameters) -> Result<f64> {
    unimplemented!();
}

#[cfg(test)]
mod similarity_method_tests {
    use crate::utils::tensor::{cosine_similarity, minimizer_l2_norm, strobemer};

    use super::{Parameters, l2norm};
    use pretty_assertions::assert_eq;

    fn test_l2norm(p: Parameters, expected: f64)  {
        let res = l2norm(p);
        assert!(res.is_ok());
        assert_eq!(res.unwrap(), expected);
    }

    // START l2_norm TESTS
    #[test]
    fn test_l2norm_equal_k1() {
        // k=1, same number of each character, should result in l2norm = 0.0
        let p = Parameters {
            k: 1,
            w: 1, // irrelevant for l2norm
            base_sequence: "ACGT".chars().collect(),
            modified_sequence: "GCAT".chars().collect(),
        };
        test_l2norm(p, 0.0);
    }

    #[test]
    fn test_l2norm_equal_k2() {
        // k=2, same number of 2-mers, should result in l2norm = 0.0
        let p = Parameters {
            k: 2,
            w: 1, // irrelevant for l2norm
            base_sequence: "ACGTA".chars().collect(),
            modified_sequence: "CGTAC".chars().collect(),
        };
        test_l2norm(p, 0.0);
    }

    #[test]
    fn test_l2norm_unequal_k1() {
        // k=1, unequal character set, should result in l2norm = sqrt(2)
        let p = Parameters {
            k: 1,
            w: 1, // irrelevant for l2norm
            base_sequence: "GTTTTT".chars().collect(),
            modified_sequence: "ATTTTT".chars().collect(),
        };
        test_l2norm(p, 1.4142135623730951);
    }

    #[test]
    fn test_l2norm_unequal_k2() {
        // k=2, unequal 2-mer set, should result in l2norm = sqrt(2)
        let p = Parameters {
            k: 1,
            w: 1, // irrelevant for l2norm
            base_sequence: "GTTTTT".chars().collect(),
            modified_sequence: "ATTTTT".chars().collect(),
        };
        test_l2norm(p, 1.4142135623730951);
    }

    fn test_cosine_similarity(p: Parameters, expected: f64) {
        let res = cosine_similarity(p);
        assert!(res.is_ok());
        assert_eq!(res.unwrap(), expected);
    }

    #[test]
    fn test_cosine_similarity_equal_k1() {
        // k=2, unequal 2-mer set, should result in l2norm = sqrt(2)
        let p = Parameters {
            k: 1,
            w: 1, // irrelevant for l2norm
            base_sequence: "ACCT".chars().collect(),
            modified_sequence: "ACCG".chars().collect(),
        };
        test_cosine_similarity(p, 1.4142135623730951);
    }

    fn test_minimizer_l2_norm(p: Parameters, expected: f64) {
        let res = minimizer_l2_norm(p);
        assert!(res.is_ok());
        assert_eq!(res.unwrap(), expected);
    }

    #[test]
    fn test_minimizer_l2_norm_equal_k1() {
        // Some less trivial example?
        let p = Parameters {
            k: 1,
            w: 1,
            base_sequence: vec!['A', 'C', 'G', 'T'],
            modified_sequence: vec!['T', 'C', 'A', 'T'],
        };
        test_minimizer_l2_norm(p, 1.4142135623730951);
    }


    fn test_strobemer(p: Parameters, expected: f64) {
        let res = strobemer(p);
        assert!(res.is_ok());
        assert_eq!(res.unwrap(), expected);
    }
    #[test]
    fn test_strobemer_trivial_k1() {
        let p = Parameters {
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
    fn test_generate_kmers() {
        let kmers = generate_kmers(&['A', 'C', 'G', 'T'], 1);
        assert_eq!(kmers.len(), 4);
        for (i, c) in "ACGT".chars().enumerate() {
            assert_eq!(kmers[i], [c]);
        }

        let kmers = generate_kmers(&['A', 'C', 'G', 'T'], 2);
        assert_eq!(kmers.len(), 3);
        assert_eq!(kmers[0], ['A', 'C']);
        assert_eq!(kmers[1], ['C', 'G']);
        assert_eq!(kmers[2], ['G', 'T']);
    }
}