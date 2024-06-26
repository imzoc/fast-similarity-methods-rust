use utils::tensor::Tensor;
use utils::kmer::generate_kmers;


fn main() {
    let sequence = String::from("ACTGACTG");
    let k = 2;
    let kmers = generate_kmers(&sequence, k);
    println!("{:?}", kmers);

    let mut tensor = Tensor::new(vec![4, 4], 0);
    println!("Tensor shape: {:?}", tensor.shape());
    match tensor.populate(&kmers) {
        Ok(_) => println!("Tensor populated successfully."),
        Err(e) => println!("Error populating tensor: {}", e),
    }
    println!("K-mer frequency vector: {:?}", tensor.data_ref());
}