use anyhow::{anyhow, bail, Result};

use crate::utils::representation_methods;

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
        let kmers = representation_methods::generate_kmers_with_step(sequence, k, step)?;
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