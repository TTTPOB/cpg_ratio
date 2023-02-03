use bio::io::fasta::{self, Record};
use clap::{Parser, Subcommand, ValueEnum};
use std::io::{self, prelude::*, BufWriter};

#[derive(Parser)]
struct App {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    Cg(CgArgs),
    #[clap(name = "cpg")]
    CpG(CpGArgs),
}

#[derive(Parser, Debug, Clone, Copy)]
struct CgArgs {
    #[clap(short, long)]
    window_size: u64,
}

#[derive(Parser, Debug, Clone, Copy)]
struct CpGArgs {
    #[clap(short, long)]
    window_size: u64,
}

struct BdgRecord {
    seqname: String,
    start: u64,
    end: u64,
    score: f64,
}
impl BdgRecord {
    fn new(seqname: String, start: u64, end: u64, score: f64) -> Self {
        Self {
            seqname,
            start,
            end,
            score,
        }
    }
    fn writeln(&self, writer: &mut impl Write) -> io::Result<()> {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}",
            self.seqname, self.start, self.end, self.score
        )
    }
    fn write(&self, writer: &mut impl Write) -> io::Result<()> {
        write!(
            writer,
            "{}\t{}\t{}\t{}",
            self.seqname, self.start, self.end, self.score
        )
    }
}

fn count_bigram(seq: &[u8], bigram: &[u8; 2]) -> u64 {
    seq.windows(2).filter(|w| w == bigram).count() as u64
}
fn count_base(seq: &[u8], bases: &[u8]) -> u64 {
    seq.iter().filter(|b| bases.contains(b)).count() as u64
}

struct CpGCounter {
    record: Record,
    window_size: u64,
    current_idx: usize,
}

impl Iterator for CpGCounter {
    type Item = BdgRecord;
    fn next(&mut self) -> Option<Self::Item> {
        if self.current_idx + self.window_size as usize > self.record.seq().len() {
            return None;
        }
        let start = self.current_idx;
        let end = self.current_idx + self.window_size as usize;
        let seq = &self.record.seq()[start..end];
        let cpg_count = count_bigram(seq, b"CG");
        let total_bigram_count = seq.len() - 1;
        let score = (cpg_count as f64) / (total_bigram_count as f64);
        self.current_idx += self.window_size as usize;
        Some(BdgRecord::new(
            self.record.id().to_string(),
            start as u64,
            end as u64,
            score,
        ))
    }
}

struct CGCounter {
    record: Record,
    window_size: u64,
    current_idx: usize,
}
impl Iterator for CGCounter {
    type Item = BdgRecord;
    fn next(&mut self) -> Option<Self::Item> {
        if self.current_idx + self.window_size as usize > self.record.seq().len() {
            return None;
        }
        let start = self.current_idx;
        let end = self.current_idx + self.window_size as usize;
        let seq = &self.record.seq()[start..end];
        let cg_count = count_base(seq, b"CG");
        let total_base_count = seq.len();
        let score = (cg_count as f64) / (total_base_count as f64);
        self.current_idx += self.window_size as usize;
        Some(BdgRecord::new(
            self.record.id().to_string(),
            start as u64,
            end as u64,
            score,
        ))
    }
}


fn main() {
    let app = App::parse();
    let input = io::stdin().lock();
    let mut output = BufWriter::new(io::stdout());
    let fasta = fasta::Reader::new(input);
    let mut records = fasta.records();
    while let Some(Ok(record)) = records.next() {
        match app.command {
            Commands::Cg(args) => {
                let mut counter = CGCounter {
                    record,
                    window_size: args.window_size,
                    current_idx: 0,
                };
                while let Some(record) = counter.next() {
                    record.writeln(&mut output).unwrap();
                }
            },
            Commands::CpG(args) => {
                let mut counter = CpGCounter {
                    record,
                    window_size: args.window_size,
                    current_idx: 0,
                };
                while let Some(record) = counter.next() {
                    record.writeln(&mut output).unwrap();
                }
            }
        }
    }
}
