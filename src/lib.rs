use std::alloc::{alloc, dealloc, Layout};


// constants for hard-coded direction values.
const NO_DIRECTION: i32 = -1;
const LEFT_DIRECTION: i32 = 0;
const DIAGONAL_DIRECTION: i32 = 1;
const UP_DIRECTION: i32 = 2;


/// Comparison function for linear gap alignment. It
/// expects three scores: the "left," "diagonal," and "up" scores.
/// 
/// It returns the maximum score (i32) and the direction (i32) representing a
/// hard-coded direction (see constants).
fn score_comparison(left: i32, diagonal: i32, up: i32, is_local: bool) -> (i32, i32) {
    // maximum and direction values
    let mut max;
    let mut direction;

    // compare the "left" score vs. the "diagonal" score, then the "up" score
    if left > diagonal {
        max = left;
        direction = LEFT_DIRECTION;
    } else {
        max = diagonal;
        direction = DIAGONAL_DIRECTION;
    }

    if up > max {
        max = up;
        direction = UP_DIRECTION;
    }

    // if local alignment is specified, compare maximum against "0." If it turns out that
    // the "0" score wins, then that means, in local alignment, it should return "no direction."
    if is_local {
        max = if max > 0 {max} else {0};
        return (max, NO_DIRECTION);
    }

    (max, direction)
}


/// Naive version of "scoring_matrix_linear_gap" that uses only safe functions and
/// constructs in Rust.
pub fn scoring_matrix_linear_gap_naive(seq1: &[u8], seq2: &[u8], match_score: i32, mismatch_penalty: i32, gap_penalty: i32, is_local: bool) -> (i32, (usize, usize), Vec<Vec<(i32, i32)>>) {
    // prepare sizes
    let n = seq1.len() + 1;
    let m = seq2.len() + 1;

    // prepare max, i_max, j_max for local alignment
    let mut max = 0;
    let mut i_max = 0;
    let mut j_max = 0;

    // create matrix
    let mut matrix: Vec<Vec<(i32, i32)>> = Vec::new();

    // initialize first row
    let mut first_row = Vec::new();

    // push first cell
    first_row.push((0, NO_DIRECTION));

    for j in 1..m {
        if is_local {
            first_row.push((0, NO_DIRECTION));
        } else {
            // calculate the current score for the first row
            let value = j as i32 * gap_penalty;

            // set the value to the ptr
            first_row.push((value, LEFT_DIRECTION));
        }
    }

    matrix.push(first_row);

    // for each row, first push the "first column" cell, then
    // push each subsequent cell
    for i in 1..n {
        // create row
        let mut row = Vec::new();

        row.push((i as i32 * gap_penalty, UP_DIRECTION));

        // update max for local alignment
        if is_local {
            if row[0].0 > max {
                max = row[0].0;
                i_max = i;
                j_max = 0;
            }
        }

        // add each column
        for j in 1..m {
            row.push(
                score_comparison(
                row[j - 1].0 + gap_penalty,
                matrix[i - 1][j - 1].0 + if seq1[i - 1] == seq2[j - 1] {match_score} else {mismatch_penalty},
                matrix[i - 1][j].0 + gap_penalty,
                is_local
                )
            );

            // update max for local alignment
            if is_local {
                if row[j].0 > max {
                    max = row[j].0;
                    i_max = i;
                    j_max = j;
                }
            }
        }

        matrix.push(row);
    }

    if is_local {
        return (max, (i_max, j_max), matrix);
    }

    (matrix[n - 1][m - 1].0, (n - 1, m - 1), matrix)
}


/// Build scoring matrix for either global or local alignment
/// with a linear gap scoring system.
pub fn scoring_matrix_linear_gap(seq1: &[u8], seq2: &[u8], match_score: i32, mismatch_penalty: i32, gap_penalty: i32, is_local: bool) -> (i32, (usize, usize), Vec<Vec<(i32, i32)>>) {
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
        vec_to_add.push((0, NO_DIRECTION));
        *prev_ptr.add(0) = 0;

        // set the memory pointer values to the first row of penalties
        for j in 1..m {
            if is_local {
                *prev_ptr.add(j) = 0;
                vec_to_add.push((0, NO_DIRECTION));
            } else {
                // calculate the current score for the first row
                let value = j as i32 * gap_penalty;

                // set the value to the ptr
                *prev_ptr.add(j) = value;
                vec_to_add.push((value, LEFT_DIRECTION));
            }
        }

        vectors.push(vec_to_add);
    }


    // temp variables holding left, diagonal, and top
    // scores, respectively
    let mut temp_left;
    let mut temp_diagonal;
    let mut temp_up;

    // temp variables for current max
    // and the prev score[0] value
    let mut curr;
    let mut prev;
    
    // obtain the optimal alignment score
    // calculate row-by-row dynamic programming matrix
    for i in 1..n {
        // "left-most" score (the first cell in the row)
        temp_left = if is_local {0} else {i as i32 * gap_penalty};
        
        // update max for local alignment
        if is_local {
            if temp_left > max {
                max = temp_left;
                i_max = i;
                j_max = 0;
            }
        }

        // prepare a vector to add to vectors
        let mut vec_to_add = Vec::with_capacity(m);
        
        // add gap penalty and direction; if local, then indicate no direction
        let obj_to_add = if is_local {(temp_left, NO_DIRECTION)} else {(temp_left, UP_DIRECTION)};
        vec_to_add.push(obj_to_add);

        // prepare to compare against diagonal and up scores
        prev = temp_left;
        temp_left += gap_penalty;


        // calculate scores
        for j in 1..m {
            unsafe {
                // get diagonal score
                temp_diagonal = *prev_ptr.add(j - 1) + (if seq1.get_unchecked(i - 1) == seq2.get_unchecked(j - 1) {match_score} else {mismatch_penalty});
    
                // get up score
                temp_up = *prev_ptr.add(j) + gap_penalty;
            }

            // get max score and directions
            let output_tup = score_comparison(temp_left, temp_diagonal, temp_up, is_local);
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
            if is_local {
                if curr > max {
                    max = curr;
                    i_max = i;
                    j_max = j;
                }
            }

            // update the "left" score to the current score
            temp_left = curr + gap_penalty;
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
/// the full scoring matrix and a starting position.
pub fn reconstruct_alignment(seq1: &[u8], seq2: &[u8], matrix: Vec<Vec<(i32, i32)>>, starting_pos: (usize, usize)) -> Vec<String> {
    // get starting position
    let (mut i, mut j) = starting_pos;

    // create vector of alignments
    let mut alignments = vec![String::new(), String::new()];

    // loop through the dp matrix until a cell with "NO_DIRECTION" is found
    while matrix[i][j].1 != NO_DIRECTION {
        // get the direction of the current cell
        let (_, direction) = matrix[i][j];

        // go left, diagonal, or up
        if direction == LEFT_DIRECTION {
            // left direction - add a gap to the first sequence,
            // and iterate the current cell back one column
            alignments[0].push('-');
            alignments[1].push(seq2[j - 1] as char);
            
            j -= 1;
        } else if direction == DIAGONAL_DIRECTION {
            // diagonal direction - match (or mismatch) the two
            // characters, and iterate current cell back one column
            // and one row
            alignments[0].push(seq1[i - 1] as char);
            alignments[1].push(seq2[j - 1] as char);
            
            i -= 1;
            j -= 1;
        } else if direction == UP_DIRECTION {
            // up direction - add a gap to the second sequence,
            // and iterate the current cell back one row
            alignments[0].push(seq1[i - 1] as char);
            alignments[1].push('-');
            
            i -= 1;
        } else {
            break;
        }
    }

    // reconstruct and reverse the string alignments
    alignments[0] = alignments[0].chars().rev().collect();
    alignments[1] = alignments[1].chars().rev().collect();

    // return alignments
    alignments
}