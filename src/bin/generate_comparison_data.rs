use std::time::Instant;
use std::io::{self, Write};
use bio::io::fasta;
use anyhow::Result;
use clap::Parser;
use rusqlite::{params, Connection};


#[derive(Debug, Parser)]
#[command(about, author, version)]
/// Similarity Methods
struct Args {
    #[arg(short='r', long)]
    references_file: String,
    #[arg(short='q', long)]
    query_file: String,
    #[arg(short='m', long)]
    similarity_method: String,
    #[arg(long, default_value_t = 1)]
    step: usize,
    
    order: usize,
    #[arg(short='l', long)]
    strobe_length: usize,
    #[arg(long="w-gap")]
    strobe_window_gap: usize,
    #[arg(long="w_len")]
    strobe_window_length: usize,
}

#[allow(unused_imports)]
use alignment_free_methods::utils::sequence;
use alignment_free_methods::utils::representation_methods;



// --------------------------------------------------
fn main() {
    if let Err(e) = run(Args::parse()) {
        eprintln!("{e}");
        std::process::exit(1);
    }
}

// --------------------------------------------------
// See this repo's README file for pseudocode
fn run(args: Args) -> Result<()> {
    let start = Instant::now();

    let conn = Connection::open("seeds.db")?;
    conn.execute(
        "CREATE TABLE IF NOT EXISTS seeds (
            id INTEGER PRIMARY KEY,
            query_name TEXT NOT NULL,
            reference_name TEXT NOT NULL,
            seed_name TEXT NOT NULL,
            score TEXT NOT NULL
        )",
        [],
    )?;
    let mut reference_reader =
        fasta::Reader::from_file(&args.references_file)?;
    let mut query_reader = 
        fasta::Reader::from_file(&args.query_file)?;

    let mut i = 0;
    for query_record in query_reader.records() {
        let query_record = query_record?;
        let query_seq = query_record.seq();

        for reference_record in reference_reader.records() {
            i += 1;
            let reference_record = reference_record?;
            let reference_seq = reference_record.seq();
            
            let estimated_distance: f64 = representation_methods::strobemer_similarity::<&str>(
                query_seq,
                reference_seq,
                &args.similarity_method,
                args.order.clone(),
                args.strobe_length.clone(),
                args.strobe_window_gap.clone(),
                args.strobe_window_length.clone(),
                args.step.clone()
            )?;

            // enter score into database
            let seed_name = format!("({},{},{},{},{})-strobemers",
                &args.order,
                &args.strobe_length,
                &args.strobe_window_gap,
                &args.strobe_window_length,
                &args.step
            );
            let query_name = query_record.id();
            let reference_name = reference_record.id();
            conn.execute(
                "INSERT INTO seeds (query_name, reference_name, seed_name, score) VALUES (?1, ?2, ?3, ?4)",
                params![query_name, reference_name, seed_name, estimated_distance],
            )?;
            print!("\rString pair #{:?} has been processed!", i);
            io::stdout().flush().unwrap();
        }

    }

    /*
    let duration = start.elapsed();
    let mut file = File::create(format!("{}{}", &args.outfile, ".benchmark.log"))?;
    write!(file, "{:?}", duration)?;

    */

    Ok(())
}