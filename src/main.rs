use bio::io::fasta::{self, Record};
use clap::{Parser, Subcommand};
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
#[derive(Debug, Clone, Copy)]
enum ScoreTypes {
    CpG,
    CG,
}
trait Scoring {
    fn score(&self, seq: &[u8]) -> f64;
}
impl Scoring for ScoreTypes {
    fn score(&self, seq: &[u8]) -> f64 {
        match self {
            ScoreTypes::CpG => {
                let cg_count = count_bigram(seq, b"CG");
                let total_bigram_count = seq.len() as u64 - 1;
                (cg_count as f64) / (total_bigram_count as f64)
            }
            ScoreTypes::CG => {
                let cg_count = count_base(seq, b"CG");
                let total_base_count = seq.len() as u64;
                (cg_count as f64) / (total_base_count as f64)
            }
        }
    }
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

struct CGCounter<T: Scoring> {
    record: Record,
    window_size: u64,
    current_idx: usize,
    score_type: T,
}
impl<T: Scoring> CGCounter<T> {
    fn start(&self) -> usize {
        self.current_idx
    }
    fn end(&self) -> usize {
        std::cmp::min(
            self.current_idx + self.window_size as usize,
            self.record.seq().len(),
        )
    }
    fn step(&mut self) {
        self.current_idx += self.window_size as usize;
    }
    fn seq(&self) -> &[u8] {
        &self.record.seq()[self.start()..self.end()]
    }
    fn score(&self) -> f64 {
        self.score_type.score(self.seq())
    }
}

impl<T: Scoring> Iterator for CGCounter<T> {
    type Item = BdgRecord;
    fn next(&mut self) -> Option<Self::Item> {
        if self.start() >= self.record.seq().len() {
            return None;
        }
        let bdg = BdgRecord::new(
            self.record.id().to_string(),
            self.start() as u64,
            self.end() as u64,
            self.score(),
        );
        self.step();
        Some(bdg)
    }
}

fn main() {
    let app = App::parse();
    let input = io::stdin().lock();
    let mut output = BufWriter::new(io::stdout());
    let fasta = fasta::Reader::new(input);
    let mut records = fasta.records();
    while let Some(Ok(record)) = records.next() {
        let bdg_records = match &app.command {
            &Commands::CpG(args) => CGCounter {
                record,
                window_size: args.window_size,
                current_idx: 0,
                score_type: ScoreTypes::CpG,
            },
            &Commands::Cg(args) => CGCounter {
                record,
                window_size: args.window_size,
                current_idx: 0,
                score_type: ScoreTypes::CG,
            },
        };
        for bdg in bdg_records {
            bdg.writeln(&mut output).unwrap();
        }
    }
}
