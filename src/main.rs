use nw_sw_align::{scoring_matrix_linear_gap, scoring_matrix_linear_gap_naive, reconstruct_alignment};
use std::time::Instant;
use rand::prelude::*;


fn main() {
    // set of letters
    let letters = ["A", "C", "T", "G"];

    // random generator
    let mut rng = rand::rng();

    // length of sequence for seq1 and seq2 respectively
    let sequence_length_1 = 50000;
    let sequence_length_2 = 50000;


    // sample case
    let seq1: String =  (0..sequence_length_1).map(|_| *letters.choose(&mut rng).unwrap()).collect();
    let seq2: String = (0..sequence_length_2).map(|_| *letters.choose(&mut rng).unwrap()).collect();
    let seq1_bytes = &seq1.as_bytes();
    let seq2_bytes = &seq2.as_bytes();

    let match_score = 1;
    let mismatch_penalty = -1;
    let gap_penalty = -2;
    let is_local = false;

    // measure elapsed time start for optimized implementation
    let start = Instant::now();

    // get score, starting position, and matrix
    let (score, starting_pos, matrix) = scoring_matrix_linear_gap(seq1_bytes, seq2_bytes, match_score, mismatch_penalty, gap_penalty, is_local);
    let _alignments = reconstruct_alignment(seq1_bytes, seq2_bytes, matrix, starting_pos);

    // get elapsed time final
    let duration = start.elapsed();

    // print out alignments for optimized function
    println!();
    println!("Optimized Approach:");
    println!("Score: {}", score);
    // println!("Alignment:");
    // println!("{}", alignments[0]);
    // println!("{}", alignments[1]);
    println!("Duration: {:?}", duration);
    
    // measure elapsed time start for naive implementation
    let start = Instant::now();
    
    // get score, starting position, and matrix
    let (score, starting_pos, matrix) = scoring_matrix_linear_gap_naive(seq1_bytes, seq2_bytes, match_score, mismatch_penalty, gap_penalty, is_local);
    let _alignments = reconstruct_alignment(seq1_bytes, seq2_bytes, matrix, starting_pos);
    
    // get elapsed time final
    let duration = start.elapsed();
    
    // print out alignments for optimized function
    println!();
    println!("Naive Approach:");
    println!("Score: {}", score);
    // println!("Alignment:");
    // println!("{}", alignments[0]);
    // println!("{}", alignments[1]);
    println!("Duration: {:?}", duration);
}
