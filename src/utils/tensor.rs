use crate::utils::sequence;
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
    pub fn new(shape: Vec<usize>, default_value: usize) -> Self {
        let size = shape.iter().product();
        let data = vec![default_value; size];
        Tensor { data, shape }
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


    pub fn populate(&mut self, kmers: &Vec<Vec<char>>) -> Result<(), String> {
        let dimensions = self.shape.len();
        const ALPHABET: [char; 4] = ['A', 'C', 'T', 'G'];

        for kmer in kmers {
            if kmer.len() != dimensions {
                return Err(format!("K-mer length {} does not match tensor dimensions {}", kmer.len(), dimensions));
            }

            // Convert k-mer to multi-dimensional index
            let indices: Vec<usize> = kmer.iter().map(|&c| match c {
                'A' => 0,
                'C' => 1,
                'T' => 2,
                'G' => 3,
                _ => 666,
            }).collect();

            // Calculate the position in the tensor
            let index = self.calculate_index(&indices)
                .ok_or_else(|| format!("Failed to calculate index for k-mer: {:?}", kmer))?;

            // Increment the value at the calculated position
            self.data[index] += 1;
        }

        Ok(())
    }
}

pub fn l2norm(params: Parameters) -> f64 {
    let chars1_tensor = kmer_tensor(&params.base_sequence, params.k);
    let chars2_tensor = kmer_tensor(&params.modified_sequence, params.k);

    let chars1_vector = chars1_tensor.data_ref();
    let chars2_vector = chars2_tensor.data_ref();

    let mut distance = 0;
    for (x1, x2) in chars1_vector.iter().zip(chars2_vector.iter()) {
        distance += (x1.max(x2) - x1.min(x2)).pow(2);
    }
    let distance_f64 = distance as f64;
    distance_f64.sqrt()
}

pub fn kmer_tensor(sequence: &Vec<char>, k: usize) -> Tensor {
    let kmers = sequence::generate_kmers(&sequence, k);

    let mut tensor = Tensor::new(vec![4; k], 0);
    tensor.populate(&kmers).expect("Failed to populate tensor");
    tensor
}