use std::time::Instant;
use std::io::{self, Write};
use std::path::Path;
use bio::io::fasta;
use anyhow::Result;
use clap::Parser;
use rusqlite::{params, Connection};

use alignment_free_methods::utils::representation_methods;

#[derive(Debug, Parser)]
#[command(about, author, version)]
struct Args {
    #[arg(short='r', long)]
    references_file: String,
    #[arg(short='q', long)]
    query_file: String,
    #[arg(short='d')]
    database: String,
    #[arg(short='m', long)]
    similarity_method: String,
    #[arg(short='s', default_value_t = 1)]
    step: usize,
    
    #[arg(short='o')]
    order: usize,
    #[arg(short='l', long)]
    strobe_length: usize,
    #[arg(long="w-gap")]
    strobe_window_gap: usize,
    #[arg(long="w-len")]
    strobe_window_length: usize,
}

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
    let project_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let db_file = Path::new(&project_dir)
        .join("tests/outputs")
        .join(&args.database);
    let conn = Connection::open(db_file)?;
    conn.execute(
        "CREATE TABLE IF NOT EXISTS seeds (
            id INTEGER PRIMARY KEY,
            query_name TEXT NOT NULL,
            reference_name TEXT NOT NULL,
            seed_name TEXT NOT NULL,
            score TEXT NOT NULL,
            time TEXT NOT NULL
        )",
        [],
    )?;
    let query_file = Path::new(&project_dir)
        .join("tests/inputs")
        .join(&args.query_file);
    let query_reader = fasta::Reader::from_file(query_file)?;

    let mut i = 0;
    for query_record in query_reader.records() {
        let query_record = query_record?;
        let query_seq = query_record.seq();

        let reference_file = Path::new(&project_dir)
            .join("tests/inputs")
            .join(&args.references_file);
        let reference_reader = fasta::Reader::from_file(reference_file)?;
        for reference_record in reference_reader.records() {
            i += 1;
            let reference_record = reference_record?;
            let reference_seq = reference_record.seq();
            
            let start = Instant::now();
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
            let duration = start.elapsed().subsec_millis();

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
            println!("{}", estimated_distance);
            conn.execute(
                "INSERT INTO seeds (query_name, reference_name, seed_name, score, time) VALUES (?1, ?2, ?3, ?4, ?5)",
                params![query_name, reference_name, seed_name, estimated_distance, duration],
            )?;
            print!("\rString pair #{:?} has been processed!", i);
            io::stdout().flush().unwrap();
        }

    }
    Ok(())
}