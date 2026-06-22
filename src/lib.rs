use std::alloc::{alloc, dealloc, Layout};

// constants for hard-coded direction values.
const NO_DIRECTION: i32 = -1;
const LEFT_DIRECTION: i32 = 0;
const DIAGONAL_DIRECTION: i32 = 1;
const UP_DIRECTION: i32 = 2;

/// Comparison function for linear gap alignment.
fn score_comparison(left: i32, diagonal: i32, up: i32, is_local: bool) -> (i32, i32) {
    let mut max;
    let mut direction;

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

    if is_local {
        max = if max > 0 { max } else { 0 };
        return (max, NO_DIRECTION);
    }

    (max, direction)
}

/// Naive version that uses a single safe flat 1D vector instead of nested vectors.
pub fn scoring_matrix_linear_gap_naive(
    seq1: &[u8],
    seq2: &[u8],
    match_score: i32,
    mismatch_penalty: i32,
    gap_penalty: i32,
    is_local: bool,
) -> (i32, (usize, usize), Vec<(i32, i32)>) {
    let n = seq1.len() + 1;
    let m = seq2.len() + 1;

    let mut max = 0;
    let mut i_max = 0;
    let mut j_max = 0;

    // Allocate single flat contiguous vector up front
    let mut matrix = Vec::with_capacity(n * m);

    // Initialize first row
    matrix.push((0, NO_DIRECTION));
    for j in 1..m {
        if is_local {
            matrix.push((0, NO_DIRECTION));
        } else {
            matrix.push((j as i32 * gap_penalty, LEFT_DIRECTION));
        }
    }

    // Populate matrix row by row using flat index mapping: (i * m) + j
    for i in 1..n {
        let first_col_val = i as i32 * gap_penalty;
        matrix.push((first_col_val, UP_DIRECTION));

        if is_local && first_col_val > max {
            max = first_col_val;
            i_max = i;
            j_max = 0;
        }

        for j in 1..m {
            let left_score = matrix[(i * m) + (j - 1)].0 + gap_penalty;
            let diag_score = matrix[((i - 1) * m) + (j - 1)].0
                + if seq1[i - 1] == seq2[j - 1] { match_score } else { mismatch_penalty };
            let up_score = matrix[((i - 1) * m) + j].0 + gap_penalty;

            let output_tup = score_comparison(left_score, diag_score, up_score, is_local);
            matrix.push(output_tup);

            if is_local && output_tup.0 > max {
                max = output_tup.0;
                i_max = i;
                j_max = j;
            }
        }
    }

    if is_local {
        return (max, (i_max, j_max), matrix);
    }
    (matrix[(n * m) - 1].0, (n - 1, m - 1), matrix)
}

/// Optimized version utilizing row-caching pointers alongside a pre-allocated flat output vector.
pub fn scoring_matrix_linear_gap(
    seq1: &[u8],
    seq2: &[u8],
    match_score: i32,
    mismatch_penalty: i32,
    gap_penalty: i32,
    is_local: bool,
) -> (i32, (usize, usize), Vec<(i32, i32)>) {
    let n = seq1.len() + 1;
    let m = seq2.len() + 1;

    let mut max = 0;
    let mut i_max = 0;
    let mut j_max = 0;

    let prev_layout = Layout::array::<i32>(m).unwrap();
    let prev_ptr = unsafe { alloc(prev_layout) } as *mut i32;

    // Allocate safe flat vector to send directly to Python/FFI layers
    let mut matrix = Vec::with_capacity(n * m);

    unsafe {
        matrix.push((0, NO_DIRECTION));
        *prev_ptr.add(0) = 0;

        for j in 1..m {
            if is_local {
                *prev_ptr.add(j) = 0;
                matrix.push((0, NO_DIRECTION));
            } else {
                let value = j as i32 * gap_penalty;
                *prev_ptr.add(j) = value;
                matrix.push((value, LEFT_DIRECTION));
            }
        }
    }

    let mut temp_left;
    let mut temp_diagonal;
    let mut temp_up;
    let mut curr;
    let mut prev;

    for i in 1..n {
        temp_left = if is_local { 0 } else { i as i32 * gap_penalty };
        if is_local && temp_left > max {
            max = temp_left;
            i_max = i;
            j_max = 0;
        }

        let obj_to_add = if is_local { (temp_left, NO_DIRECTION) } else { (temp_left, UP_DIRECTION) };
        matrix.push(obj_to_add);

        prev = temp_left;
        temp_left += gap_penalty;

        for j in 1..m {
            unsafe {
                temp_diagonal = *prev_ptr.add(j - 1)
                    + (if seq1.get_unchecked(i - 1) == seq2.get_unchecked(j - 1) {
                        match_score
                    } else {
                        mismatch_penalty
                    });
                temp_up = *prev_ptr.add(j) + gap_penalty;
            }

            let output_tup = score_comparison(temp_left, temp_diagonal, temp_up, is_local);
            curr = output_tup.0;
            matrix.push(output_tup);

            unsafe {
                *prev_ptr.add(j - 1) = prev;
            }
            prev = curr;

            if is_local && curr > max {
                max = curr;
                i_max = i;
                j_max = j;
            }
            temp_left = curr + gap_penalty;
        }

        unsafe {
            *prev_ptr.add(m - 1) = prev;
        }
    }

    unsafe {
        dealloc(prev_ptr as *mut u8, prev_layout);
    }

    if is_local {
        return (max, (i_max, j_max), matrix);
    }
    (matrix[(n * m) - 1].0, (n - 1, m - 1), matrix)
}


/// The absolute extreme optimization: operates entirely on raw pointers for both
/// the row-cache and the massive output matrix, converting to a flat Vec only at the end.
pub fn scoring_matrix_linear_gap_optimized_extreme(
    seq1: &[u8],
    seq2: &[u8],
    match_score: i32,
    mismatch_penalty: i32,
    gap_penalty: i32,
    is_local: bool,
) -> (i32, (usize, usize), Vec<(i32, i32)>) {
    let n = seq1.len() + 1;
    let m = seq2.len() + 1;
    let total_elements = n * m;

    let mut max = 0;
    let mut i_max = 0;
    let mut j_max = 0;

    // 1. Allocate the small 1D row-caching buffer
    let prev_layout = Layout::array::<i32>(m).unwrap();
    let prev_ptr = unsafe { alloc(prev_layout) } as *mut i32;
    if prev_ptr.is_null() {
        std::alloc::handle_alloc_error(prev_layout);
    }

    // 2. Allocate the massive output matrix directly on the heap as a raw block
    let matrix_layout = Layout::array::<(i32, i32)>(total_elements).unwrap();
    let matrix_ptr = unsafe { alloc(matrix_layout) } as *mut (i32, i32);
    if matrix_ptr.is_null() {
        // Clean up row cache if matrix allocation fails
        unsafe { dealloc(prev_ptr as *mut u8, prev_layout); }
        std::alloc::handle_alloc_error(matrix_layout);
    }

    // Initialize first row using raw pointers
    unsafe {
        *matrix_ptr.add(0) = (0, NO_DIRECTION);
        *prev_ptr.add(0) = 0;

        for j in 1..m {
            if is_local {
                *prev_ptr.add(j) = 0;
                *matrix_ptr.add(j) = (0, NO_DIRECTION);
            } else {
                let value = j as i32 * gap_penalty;
                *prev_ptr.add(j) = value;
                *matrix_ptr.add(j) = (value, LEFT_DIRECTION);
            }
        }
    }

    let mut temp_left;
    let mut temp_diagonal;
    let mut temp_up;
    let mut curr;
    let mut prev;

    // Main row-by-row dynamic programming execution block
    for i in 1..n {
        let row_offset = i * m;
        temp_left = if is_local { 0 } else { i as i32 * gap_penalty };
        
        if is_local && temp_left > max {
            max = temp_left;
            i_max = i;
            j_max = 0;
        }

        let obj_to_add = if is_local { (temp_left, NO_DIRECTION) } else { (temp_left, UP_DIRECTION) };
        unsafe {
            *matrix_ptr.add(row_offset) = obj_to_add;
        }

        prev = temp_left;
        temp_left += gap_penalty;

        for j in 1..m {
            unsafe {
                // Fetch scores sequentially via raw pointers
                temp_diagonal = *prev_ptr.add(j - 1)
                    + (if seq1.get_unchecked(i - 1) == seq2.get_unchecked(j - 1) {
                        match_score
                    } else {
                        mismatch_penalty
                    });
                temp_up = *prev_ptr.add(j) + gap_penalty;
            }

            let output_tup = score_comparison(temp_left, temp_diagonal, temp_up, is_local);
            curr = output_tup.0;
            
            unsafe {
                // Write directly to the massive output matrix block
                *matrix_ptr.add(row_offset + j) = output_tup;
                *prev_ptr.add(j - 1) = prev;
            }
            prev = curr;

            if is_local && curr > max {
                max = curr;
                i_max = i;
                j_max = j;
            }
            temp_left = curr + gap_penalty;
        }

        unsafe {
            *prev_ptr.add(m - 1) = prev;
        }
    }

    // Deallocate the small row-caching layout buffer
    unsafe {
        dealloc(prev_ptr as *mut u8, prev_layout);
    }

    // 3. Reconstruct the raw matrix pointer into a standard safe flat Vec for zero-copy FFI
    let matrix_vec = unsafe {
        Vec::from_raw_parts(matrix_ptr, total_elements, total_elements)
    };

    if is_local {
        return (max, (i_max, j_max), matrix_vec);
    }
    (matrix_vec[total_elements - 1].0, (n - 1, m - 1), matrix_vec)
}


/// Reconstructs a global or local alignment from the flat 1D matrix layout,
/// returning the aligned sequences along with the exact start and end coordinates.
/// 
/// Returns: (seq1_aligned, seq2_aligned, seq1_start, seq1_end, seq2_start, seq2_end)
pub fn reconstruct_alignment_with_coords(
    seq1: &[u8],
    seq2: &[u8],
    matrix: &[(i32, i32)],
    cols: usize, // This must be the value of 'm' (seq2.len() + 1)
    starting_pos: (usize, usize),
) -> (String, String, usize, usize, usize, usize) {
    let (mut i, mut j) = starting_pos;
    
    // BLAST End coordinates correspond to where our traceback loops BEGIN
    let seq1_end = i;
    let seq2_end = j;

    let mut alignments = vec![String::new(), String::new()];

    // Trace back through the flat 1D matrix layout
    while matrix[(i * cols) + j].1 != NO_DIRECTION {
        let direction = matrix[(i * cols) + j].1;

        if direction == LEFT_DIRECTION {
            alignments[0].push('-');
            alignments[1].push(seq2[j - 1] as char);
            j -= 1;
        } else if direction == DIAGONAL_DIRECTION {
            alignments[0].push(seq1[i - 1] as char);
            alignments[1].push(seq2[j - 1] as char);
            i -= 1;
            j -= 1;
        } else if direction == UP_DIRECTION {
            alignments[0].push(seq1[i - 1] as char);
            alignments[1].push('-');
            i -= 1;
        } else {
            break;
        }
    }

    // BLAST Start coordinates correspond to where our traceback loops FINISH (adding 1 for 1-based index alignment match strings)
    let seq1_start = i + 1;
    let seq2_start = j + 1;

    alignments[0] = alignments[0].chars().rev().collect();
    alignments[1] = alignments[1].chars().rev().collect();

    (alignments[0].clone(), alignments[1].clone(), seq1_start, seq1_end, seq2_start, seq2_end)
}
