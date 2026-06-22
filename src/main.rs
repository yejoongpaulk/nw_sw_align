use clap::{Parser, Subcommand};
use rand::prelude::*;
use std::time::Instant;

// Import your updated library methods
use nw_sw_align::{
    reconstruct_alignment_with_coords, 
    scoring_matrix_linear_gap, 
    scoring_matrix_linear_gap_naive, 
    scoring_matrix_linear_gap_optimized_extreme,
};

/// High-Performance Sequence Alignment Tool (Rust Engine)
#[derive(Parser)]
#[command(name = "align-engine", author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Run linear global alignment on two explicit string sequences
    Global {
        #[arg(short, long, help = "First sequence string (e.g. ATCG)")]
        seq1: String,
        #[arg(short, long, help = "Second sequence string (e.g. ATGC)")]
        seq2: String,
        #[arg(long, default_value_t = 1, help = "Match score")]
        match_score: i32,
        #[arg(long, default_value_t = -1, help = "Mismatch penalty")]
        mismatch_penalty: i32,
        #[arg(long, default_value_t = -2, help = "Gap penalty")]
        gap_penalty: i32,
    },
    /// Run linear local alignment on two explicit string sequences
    Local {
        #[arg(short, long, help = "First sequence string (e.g. ATCG)")]
        seq1: String,
        #[arg(short, long, help = "Second sequence string (e.g. ATGC)")]
        seq2: String,
        #[arg(long, default_value_t = 1, help = "Match score")]
        match_score: i32,
        #[arg(long, default_value_t = -1, help = "Mismatch penalty")]
        mismatch_penalty: i32,
        #[arg(long, default_value_t = -2, help = "Gap penalty")]
        gap_penalty: i32,
    },
    /// Execute multi-tier performance comparisons using randomized sequence datasets
    Compare {
        #[arg(short, long, default_value_t = 50000, help = "Length of the sequence parameters")]
        size: usize,
        #[arg(short, long, action = clap::ArgAction::SetTrue, help = "Execute local alignment benchmarks (Smith-Waterman) instead of global")]
        local: bool,
    },
}

fn main() {
    let cli = Cli::parse();

    match &cli.command {
        Commands::Global { seq1, seq2, match_score, mismatch_penalty, gap_penalty } => {
            run_alignment(seq1, seq2, *match_score, *mismatch_penalty, *gap_penalty, false);
        }
        Commands::Local { seq1, seq2, match_score, mismatch_penalty, gap_penalty } => {
            run_alignment(seq1, seq2, *match_score, *mismatch_penalty, *gap_penalty, true);
        }
        Commands::Compare { size, local } => {
            run_performance_comparison(*size, *local);
        }
    }
}

/// Standardized alignment executor for global/local modes with execution profiling
fn run_alignment(seq1: &str, seq2: &str, match_score: i32, mismatch_penalty: i32, gap_penalty: i32, is_local: bool) {
    let seq1_bytes = seq1.as_bytes();
    let seq2_bytes = seq2.as_bytes();
    let m = seq2_bytes.len() + 1;

    println!("Running {} Alignment...", if is_local { "Local" } else { "Global" });
    
    // Start timing the execution block
    let start_time = Instant::now();

    // 1. Calculate matrix using extreme raw pointer layout
    let (score, starting_pos, matrix) = scoring_matrix_linear_gap_optimized_extreme(
        seq1_bytes, seq2_bytes, match_score, mismatch_penalty, gap_penalty, is_local
    );
    
    // 2. Extract aligned text and coordinate boundaries via traceback
    let (align1, align2, s1_start, _s1_end, s2_start, _s2_end) = 
        reconstruct_alignment_with_coords(seq1_bytes, seq2_bytes, &matrix, m, starting_pos);
    
    // Stop timing before terminal printing begins so output latency doesn't skew benchmarks
    let execution_duration = start_time.elapsed();

    // Pipe the verified start coordinates straight into the chunk visualizer
    print_alignment_line_by_line(&align1, &align2, s1_start, s2_start, 60);
    
    println!("Final Alignment Score: {}", score);
    println!("Execution Duration: {:?}", execution_duration);
}

/// Automated random data tracking and multi-tier profiling harness
fn run_performance_comparison(size: usize, is_local: bool) {
    let letters = ["A", "C", "T", "G"];
    let mut rng = rand::rng();

    println!("Mode: {} Alignment Profiling", if is_local { "LOCAL" } else { "GLOBAL" });
    println!("Generating two randomized sequence strands of size {}...", size);
    let seq1: String = (0..size).map(|_| *letters.choose(&mut rng).unwrap()).collect();
    let seq2: String = (0..size).map(|_| *letters.choose(&mut rng).unwrap()).collect();

    let seq1_bytes = seq1.as_bytes();
    let seq2_bytes = seq2.as_bytes();
    let m = seq2_bytes.len() + 1;

    let match_score = 1;
    let mismatch_penalty = -1;
    let gap_penalty = -2;

    println!("\nExecuting comparisons (Matrix Dimensions: {} x {})...", size + 1, size + 1);

    // Tier 1: Extreme Optimization Pass
    let start = Instant::now();
    let (score_ext, starting_pos_ext, matrix_ext) = scoring_matrix_linear_gap_optimized_extreme(
        seq1_bytes, seq2_bytes, match_score, mismatch_penalty, gap_penalty, is_local
    );
    let _ = reconstruct_alignment_with_coords(seq1_bytes, seq2_bytes, &matrix_ext, m, starting_pos_ext);
    let duration_ext = start.elapsed();
    println!("\nExtreme Optimized Approach:\n  Score: {}\n  Duration: {:?}", score_ext, duration_ext);

    // Tier 2: Standard Optimization Pass
    let start = Instant::now();
    let (score_opt, starting_pos_opt, matrix_opt) = scoring_matrix_linear_gap(
        seq1_bytes, seq2_bytes, match_score, mismatch_penalty, gap_penalty, is_local
    );
    let _ = reconstruct_alignment_with_coords(seq1_bytes, seq2_bytes, &matrix_opt, m, starting_pos_opt);
    let duration_opt = start.elapsed();
    println!("\nStandard Optimized Approach:\n  Score: {}\n  Duration: {:?}", score_opt, duration_opt);

    // Tier 3: Naive Pass (Only run if matrix is reasonably sized to prevent terminal lockups)
    if size <= 30000 {
        let start = Instant::now();
        let (score_naive, starting_pos_naive, matrix_naive) = scoring_matrix_linear_gap_naive(
            seq1_bytes, seq2_bytes, match_score, mismatch_penalty, gap_penalty, is_local
        );
        let _ = reconstruct_alignment_with_coords(seq1_bytes, seq2_bytes, &matrix_naive, m, starting_pos_naive);
        let duration_naive = start.elapsed();
        println!("\nNaive Approach:\n  Score: {}\n  Duration: {:?}", score_naive, duration_naive);
    } else {
        println!("\nNaive Approach:\n  Skipped (Size > 30,000 threshold to prevent reallocation lockup).");
    }
}

/// Dynamically generates relation markers and prints the alignment track line-by-line,
/// utilizing explicit 1-based start positions provided straight from your traceback function.
pub fn print_alignment_line_by_line(
    seq1_aligned: &str, 
    seq2_aligned: &str, 
    mut s1_idx: usize, 
    mut s2_idx: usize, 
    chunk_size: usize
) {
    let seq1_chars: Vec<char> = seq1_aligned.chars().collect();
    let seq2_chars: Vec<char> = seq2_aligned.chars().collect();
    let total_len = seq1_chars.len();
    let mut start = 0;

    println!("\n=================== BLAST-STYLE VISUAL ALIGNMENT ===================");
    println!("");
    while start < total_len {
        let end = std::cmp::min(start + chunk_size, total_len);
        let chunk1: String = seq1_chars[start..end].iter().collect();
        let chunk2: String = seq2_chars[start..end].iter().collect();
        
        let chunk_m: String = seq1_chars[start..end]
            .iter()
            .zip(seq2_chars[start..end].iter())
            .map(|(c1, c2)| {
                if *c1 == '-' || *c2 == '-' { ' ' }
                else if c1 == c2 { '|' }
                else { '.' }
            })
            .collect();
            
        // Count characters in this specific line chunk window (excluding gaps)
        let s1_chars_in_chunk = chunk1.chars().filter(|&c| c != '-').count();
        let s2_chars_in_chunk = chunk2.chars().filter(|&c| c != '-').count();
        
        // Calculate the row end values safely
        let s1_line_end = s1_idx + s1_chars_in_chunk - if s1_chars_in_chunk > 0 { 1 } else { 0 };
        let s2_line_end = s2_idx + s2_chars_in_chunk - if s2_chars_in_chunk > 0 { 1 } else { 0 };
        
        // Print tracks bounded by clean visual coordinates
        println!("Seq1:  {:5} {} {:<5}", s1_idx, chunk1, s1_line_end);
        println!("Match:       {}", chunk_m);
        println!("Seq2:  {:5} {} {:<5}", s2_idx, chunk2, s2_line_end);
        println!();
        
        // Step active index parameters forward for the next line chunk window iteration
        s1_idx += s1_chars_in_chunk;
        s2_idx += s2_chars_in_chunk;
        start += chunk_size;
    }
    println!("====================================================================\n");
}
