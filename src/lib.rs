use std::alloc::{alloc, dealloc, Layout};


const NUM_DIRECTIONS: usize = 3;


/// Comparison function for linear gap alignment. It
/// expects a 3-element array, such that the first element
/// is the "left" score, the second is the "diagonal" score, and
/// the third is the "top" score.
/// 
/// Returns a 3-element boolean array, with each element
/// being True if its respective direction (i.e. first element
/// corresponds to "left", second element corresponds to "diagonal", etc.) is
/// the direction that is "maximum."
fn score_comparison(list: &[i32; NUM_DIRECTIONS], is_local: bool) -> (i32, [bool; NUM_DIRECTIONS]) {
    // first find maximum
    let mut max;

    max = if list[0] > list[1] {list[0]} else {list[1]};
    max = if max > list[2] {max} else {list[2]};

    // if local alignment is specified, compare maximum against "0"
    if is_local {
        max = if max > 0 {max} else {0};
    }

    // now, determine which directions the alignments will go
    let mut directions = [false; NUM_DIRECTIONS];

    // if local alignment and score is 0, then just return
    // a cell telling the program to stop at this cell
    if max == 0 && is_local {
        return (max, directions)
    }

    for i in 0..list.len() {
        if list[i] == max {directions[i] = true};
    }

    (max, directions)
}


/// Build scoring matrix for either global or local alignment
/// with a linear gap scoring system.
fn scoring_matrix_linear_gap(seq1: &[u8], seq2: &[u8], max_function: fn(&[i32; NUM_DIRECTIONS], bool) -> (i32, [bool; NUM_DIRECTIONS]),
match_score: i32, mismatch_penalty: i32, gap_penalty: i32, is_local: bool) -> (i32, (usize, usize), Vec<Vec<(i32, [bool; NUM_DIRECTIONS])>>) {
    // prepare sizes
    let n = seq1.len() + 1;
    let m = seq2.len() + 1;

    // prepare max, i_max, j_max
    let mut max = 0;
    let mut i_max = 0;
    let mut j_max = 0;

    // prev and curr layouts and arrays
    let prev_layout = Layout::array::<i32>(m).unwrap();

    // allocate arrays
    let prev_ptr = unsafe {alloc(prev_layout)} as *mut i32;

    // create vector of vectors
    let mut vectors = Vec::with_capacity(n);

    // initialize first row of penalties
    unsafe {
        // create first vector row
        let mut vec_to_add = Vec::with_capacity(m);

        // create very first cell (0)
        vec_to_add.push((0, [false, false, false]));
        *prev_ptr.add(0) = 0;

        // set the memory pointer values to the first row of penalties
        for j in 1..m {
            if is_local {
                *prev_ptr.add(j) = 0;
                vec_to_add.push((0, [false, false, false]));
            } else {
                // calculate the current score for the first row
                let value = j as i32 * gap_penalty;

                // set the value to the ptr
                *prev_ptr.add(j) = value;
                vec_to_add.push((value, [true, false, false]));
            }
        }

        vectors.push(vec_to_add);
    }

    // create mutable array of left, diag, and top
    // scores, respectively
    let mut scores = [0, 0, 0];


    // temp variables for current max
    // and the prev score[0] value
    let mut curr;
    let mut prev;
    
    // obtain the optimal alignment score
    // calculate row-by-row dynamic programming matrix
    for i in 1..n {
        // "left" value
        scores[0] = if is_local {0} else {i as i32 * gap_penalty};
        
        // prepare a vector to add to vectors
        let mut vec_to_add = Vec::with_capacity(m);
        
        // prepare gap penalty at current row, first column
        let obj_to_add = if is_local {(scores[0], [false, false, false])} else {(scores[0], [false, false, true])};
        vec_to_add.push(obj_to_add);

        prev = scores[0];
        scores[0] += gap_penalty;

        // calculate scores
        for j in 1..m {
            unsafe {
                // get diagonal score
                scores[1] = *prev_ptr.add(j - 1) + (if seq1[i - 1] == seq2[j - 1] {match_score} else {mismatch_penalty});
    
                // get top score
                scores[2] = *prev_ptr.add(j) + gap_penalty;
            }

            // get max score and directions
            let output_tup = max_function(&scores, is_local);
            curr = output_tup.0;

            // push current score and directions to the current row
            vec_to_add.push(output_tup);

            // write down the "prev" score
            unsafe {
                *prev_ptr.add(j - 1) = prev;
            }

            // set new "prev" score
            prev = curr;

            // set max
            if curr > max {
                max = curr;
                i_max = i;
                j_max = j;
            }

            // update the "left" score to the current score
            scores[0] = curr + gap_penalty;
        }

        // add new vector to overall vector of vectors
        vectors.push(vec_to_add);

        // update the last cell with the "left" score in the "previous row"
        unsafe {
            *prev_ptr.add(m - 1) = prev;
        }
    }

    // deallocate all unsafe memory
    unsafe {
        dealloc(prev_ptr as *mut u8, prev_layout);
    }


    if is_local {
        return (max, (i_max, j_max), vectors);
    }

    (vectors[n - 1][m - 1].0, (n - 1, m - 1), vectors)
}


/// Reconstruct a global or local alignment, given
/// the full scoring matrix.
fn reconstruct_alignment(seq1: &str, seq2: &str, matrix: Vec<Vec<(i32, [bool; NUM_DIRECTIONS])>>, starting_pos: (usize, usize), is_local: bool) -> Vec<[String; 2]> {
    // get starting position
    let (mut i, mut j) = starting_pos;

    let seq1 = seq1.as_bytes();
    let seq2 = seq2.as_bytes();

    // create vector of alignments
    let mut alignments = vec![[String::new(), String::new()]];

    // local alignment
    if is_local {
        // construct alignment until a "0" is found
        loop {
            let (score, directions) = matrix[i][j];
            if score != 0 {
                // count the number of possible paths
                let mut count = 0;

                // prefer the score that says to go "align"
                let go_left = directions[0];
                let go_diag = directions[1];
                let go_up = directions[2];

                // go left
                if go_left {
                    alignments[0][0].push('-');
                    alignments[0][1].push(seq2[j - 1] as char);
                    j -= 1;
                    continue;
                }

                if go_diag {
                    alignments[0][0].push(seq1[i - 1] as char);
                    alignments[0][1].push(seq2[j - 1] as char);
                    i -= 1;
                    j -= 1;
                    continue;
                }

                if go_up {
                    alignments[0][0].push(seq1[i - 1] as char);
                    alignments[0][1].push('-');
                    i -= 1;
                    continue;
                }
            } else {
                break;
            }
        }
    } else {
        // construct alignment until a "0" is found
        while i > 0 || j > 0 {
            let (score, directions) = matrix[i][j];
            
            // count the number of possible paths
            let mut count = 0;

            // prefer the score that says to go "align"
            let go_left = directions[0];
            let go_diag = directions[1];
            let go_up = directions[2];

            // go left
            if go_left {
                alignments[0][0].push('-');
                alignments[0][1].push(seq2[j - 1] as char);
                j -= 1;
                continue;
            }

            if go_diag {
                alignments[0][0].push(seq1[i - 1] as char);
                alignments[0][1].push(seq2[j - 1] as char);
                i -= 1;
                j -= 1;
                continue;
            }

            if go_up {
                alignments[0][0].push(seq1[i - 1] as char);
                alignments[0][1].push('-');
                i -= 1;
                continue;
            }
        }
    }

    alignments[0][0] = alignments[0][0].chars().rev().collect();
    alignments[0][1] = alignments[0][1].chars().rev().collect();

    alignments
}