use crate::utils::sequence;
use anyhow::{anyhow, bail, Result};
use std::hash::{Hash, Hasher};

use seahash::SeaHasher;

#[derive(Debug)]
pub struct Kmer {
    pub chars: Vec<char>,
    pub k: usize,
}

pub struct Tensor {
    pub data: Vec<usize>,
    pub shape: Vec<usize>,
}

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

    pub fn populate(&mut self, kmers: &Vec<Kmer>) -> Result<()> {
        let dimensions = self.shape.len();
        const ALPHABET: [char; 4] = ['A', 'C', 'T', 'G'];

        for kmer in kmers {
            if kmer.chars.len() != dimensions {
                bail!(
                    "K-mer length {} does not match tensor dimensions {}",
                    kmer.chars.len(),
                    dimensions
                );
            }

            // Convert k-mer to multi-dimensional index
            let indices: Vec<usize> = kmer
                .chars
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

pub fn generate_kmers(sequence: &[char], k: usize) -> Vec<Kmer> {
    let mut kmers = vec![];

    if k > sequence.len() {
        return kmers;
    }

    for i in 0..=(sequence.len() - k) {
        let kmer = Kmer {
            chars: sequence[i..i + k].to_vec(),
            k,
        };
        kmers.push(kmer);
    }
    kmers
}

pub fn l2norm(params: Parameters) -> Result<f64> {
    let base_sequence_kmer_vector =
        kmer_tensor(&params.base_sequence, params.k)?.to_vec();
    let modified_sequence_kmer_vector =
        kmer_tensor(&params.modified_sequence, params.k)?.to_vec();

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
        kmer_tensor(&params.base_sequence, params.k)?
            .to_vec()
            .into_iter()
            .map(|v| v as f64);
    let modified_sequence_kmer_vector =
        kmer_tensor(&params.modified_sequence, params.k)?
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

fn my_hash_function<T: Hash>(item: &T, seed: u64) -> usize {
    let mut hasher = SeaHasher::with_seeds(seed, seed, seed, seed);
    item.hash(&mut hasher);
    hasher.finish() as usize
}

pub fn minimizer_l2_norm(params: Parameters) -> Result<f64> {
    let mut base_minimizers: Vec<Kmer> = Vec::new();
    let mut mod_minimizers: Vec<Kmer> = Vec::new();

    let smallest_sequence_length = std::cmp::min(
        params.base_sequence.len(),
        params.modified_sequence.len(),
    );

    for idx in 0..smallest_sequence_length {
        let base_window = &params.base_sequence[idx..params.w].to_vec();
        let mod_window = &params.modified_sequence[idx..params.w].to_vec();

        let base_kmers = sequence::generate_kmers(base_window, params.k);
        let mod_kmers = sequence::generate_kmers(mod_window, params.k);

        let seed: u64 = 69;
        let base_minimizer = Kmer {
            chars: base_kmers
                .iter()
                .min_by_key(|item| my_hash_function(item, seed))
                .unwrap()
                .clone(),
            k: params.k,
        };
        let mod_minimizer = Kmer {
            chars: mod_kmers
                .iter()
                .min_by_key(|item| my_hash_function(item, seed))
                .unwrap()
                .clone(),
            k: params.k,
        };

        base_minimizers.push(base_minimizer);
        mod_minimizers.push(mod_minimizer);
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

pub fn kmer_tensor(sequence: &Vec<char>, k: usize) -> Result<Tensor> {
    Tensor::construct(sequence, k)
}
