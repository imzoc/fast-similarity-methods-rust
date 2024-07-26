use std::time::Instant;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use bio::io::fasta;
use anyhow::Result;
use clap::Parser;
use rusqlite::{self, params, Connection};

use alignment_free_methods;

#[derive(Debug, Parser)]
#[command(about, author, version)]
struct CommonArgs {
    #[arg(required = true, value_name = "REF FILE", help = "References FASTA file")]
    references_file: String,
    #[arg(required = true, value_name = "QUERY FILE", help = "Query FAST(A/Q) file")]
    query_file: String,
    #[arg(short='m', value_name = "METHOD()", help = "E.g. jaccard_similarity")]
    similarity_method: String,
    #[arg(short='s', default_value_t = 1, value_name = "INT PARAM")]
    step: usize
}

#[derive(Debug, Parser)]
struct StrobemerArgs {
    #[command(flatten)]
    common: CommonArgs,

    #[arg(short='p', value_name = "STRING")]
    protocol: String,
    #[arg(short='o', value_name = "INT")]
    order: usize,
    #[arg(short='l', value_name = "INT", long)]
    strobe_length: usize,
    #[arg(long="w-gap", value_name = "INT")]
    strobe_window_gap: usize,
    #[arg(long="w-len", value_name = "INT")]
    strobe_window_length: usize,
}

// --------------------------------------------------
fn main() {
    if let Err(e) = run(StrobemerArgs::parse()) {
        eprintln!("{e}");
        std::process::exit(1);
    }
}

fn initialize_comparison_db(filename: PathBuf) -> Result<Connection> {
    let conn = Connection::open(filename)?;
    conn.execute(
        "CREATE TABLE IF NOT EXISTS comparisons (
            query_name TEXT NOT NULL,
            reference_name TEXT NOT NULL,
            seed_name TEXT NOT NULL,
            score TEXT NOT NULL,
            time TEXT NOT NULL
        )",
        [],
    )?;
    Ok(conn)
}

fn initialize_seeds_db(filename: PathBuf) -> Result<Connection> {
    let conn = Connection::open(filename)?;
    conn.execute(
        "CREATE TABLE IF NOT EXISTS seeds (
            id INTEGER PRIMARY KEY,
            seq_name TEXT NOT NULL,
            seed_name TEXT NOT NULL,
            seed TEXT NOT NULL,
            time TEXT NOT NULL
        )",
        [],
    )?;
    Ok(conn)
}
fn fetch_seeds(conn: &Connection, seq_name: &str, seed_name: &str) -> Result<Vec<Vec<char>>> {
    let mut seeds = Vec::new();
    let mut stmt = conn.prepare(
        "SELECT seed FROM seeds WHERE seq_name = ?1 AND seed_name = ?2")?;
    let mut rows = stmt.query(params![seq_name, seed_name])?;

    while let Some(row) = rows.next()? {
        let value: String = row.get(0)?; // Assuming the column type is TEXT
        seeds.push(value.chars().collect());
    }
    Ok(seeds)
}

fn store_seeds(conn: &Connection, seeds: &Vec<Vec<char>>, seq_name: &str, seed_name: &str) -> Result<()> {
    for seed in seeds {
        conn.execute(
            "INSERT INTO seeds (seq_name, seed_name, seed, time) VALUES (?1, ?2, ?3, ?4)",
            params![seq_name, seed_name, String::from_iter(seed), 0],
        )?;
    }
    
    Ok(())
}

// --------------------------------------------------
// See this repo's README file for pseudocode
fn run(args: StrobemerArgs) -> Result<()> {
    let seed_name = format!("({},{},{},{},{})-strobemers",
        &args.order,
        &args.strobe_length,
        &args.strobe_window_gap,
        &args.strobe_window_length,
        &args.common.step
    );
    let project_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let comparison_db_conn = initialize_comparison_db(
        Path::new(&project_dir).join("tests/outputs/comparisons.db")
    )?;
    let seeds_db_conn = initialize_seeds_db(
        Path::new(&project_dir).join("tests/outputs/seeds.db")
    )?;

    let query_reader = fasta::Reader::from_file(
        Path::new(&project_dir)
            .join("tests/inputs")
            .join(&args.common.query_file))?;

    let mut seeds_generated = 0;
    let mut seeds_retrieved = 0;
    let mut i = 0;
    for query_record in query_reader.records() {
        let query_record = query_record?;
        let query_seeds: Vec<Vec<char>> = match seeds_db_conn.prepare(
            "SELECT 1 FROM seeds WHERE seq_name = ?1 AND seed_name = ?2")?
            .exists(params![query_record.id(), seed_name])? {
            true => {
                seeds_retrieved += 1;
                fetch_seeds(&seeds_db_conn, query_record.id(), &seed_name)?
            },
            false => {
                let seeds = alignment_free_methods::generate_strobemers(
                    query_record.seq(),
                    args.order.clone(),
                    args.strobe_length.clone(),
                    args.strobe_window_gap.clone(),
                    args.strobe_window_length.clone(),
                    args.common.step.clone(),
                    None
                )?;
                store_seeds(&seeds_db_conn, &seeds, query_record.id(), &seed_name)?;
                seeds_generated += 1;
                seeds
            }
        };
        let reference_reader = fasta::Reader::from_file(
            Path::new(&project_dir)
                .join("tests/inputs")
                .join(&args.common.references_file)
        )?;
        for reference_record in reference_reader.records() {
            i += 1;
            let reference_record = reference_record?;
            //let reference_seeds: Vec<Vec<char>> = match seeds_db_conn.prepare(
            //    "SELECT seed FROM seeds WHERE seq_name = ?1 AND seed_name = ?2")?
            //    .exists(params![reference_record.id(), seed_name])? {
            let reference_seeds: Vec<Vec<char>> = match false {
                true => {
                    seeds_retrieved += 1;
                    fetch_seeds(&seeds_db_conn, reference_record.id(), &seed_name)?
                },
                false => {
                    let seeds = alignment_free_methods::generate_strobemers(
                        reference_record.seq(),
                        args.order.clone(),
                        args.strobe_length.clone(),
                        args.strobe_window_gap.clone(),
                        args.strobe_window_length.clone(),
                        args.common.step.clone(),
                        None
                    )?;
                    // store_seeds(&seeds_db_conn, &seeds, reference_record.id(), &seed_name)?;
                    seeds_generated += 1;
                    seeds
                }
            };

            print!("\rseeds generated:{:?} | seeds retrieved:{:?}", seeds_generated, seeds_retrieved);
            io::stdout().flush().unwrap();
            
            let start = Instant::now();
            let estimated_distance: f64 = alignment_free_methods::jaccard_similarity(
                &reference_seeds,
                &query_seeds,
            )?;
            let duration = start.elapsed().subsec_millis();
            comparison_db_conn.execute(
                "INSERT OR REPLACE INTO comparisons (query_name, reference_name, seed_name, score, time) VALUES (?1, ?2, ?3, ?4, ?5)",
                params![query_record.id(), reference_record.id(), seed_name, estimated_distance, duration],
            )?;

        }
    }
    Ok(())
}