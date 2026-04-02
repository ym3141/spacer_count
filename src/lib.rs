use pyo3::prelude::*;

/// A Python module implemented in Rust.
#[pymodule]
mod pw_align {
    use pyo3::prelude::*;
    use std::cmp::{max, min};
    use std::io::{self, Write};

    #[pyfunction]
    fn correct_seq_mt(ref_seqs: Vec<String>, query_seqs: Vec<String>, max_flex: isize, threads: usize) -> Vec<String> {

        print!("    Aligning with {} threads: ", threads);
        io::stdout().flush().unwrap();

        if threads == 1 {
            let mut results = Vec::new();
            let progress_interval = max(1, query_seqs.len() / 50); // Print progress every 2% of the sequences
            for (idx, seq) in query_seqs.iter().enumerate() {
                if idx % progress_interval == 0 {
                    print!("\u{25A0}");
                    io::stdout().flush().unwrap();
                }
                let result = correct_seq(&ref_seqs, seq, max_flex);
                results.push(result);
            }
            println!(); // Move to the next line after progress bar
            return results;

        } else if threads > 1 {
            use std::thread;
            use std::sync::Arc;

            let chunk_size = (query_seqs.len() + threads - 1) / threads; // Calculate chunk size for each thread
            let mut shared_seq_chunks = Vec::new();
            for i in 0..threads {
                let start = i * chunk_size;
                let end = min(start + chunk_size, query_seqs.len());
                if start < end {
                    shared_seq_chunks.push(Arc::new(query_seqs[start..end].to_vec()));
                }
            }

            let shared_ref_seqs = Arc::new(ref_seqs);
            let mut handles = Vec::new();
            let progress_interval = chunk_size / 50 * threads; // Print progress after each chunk is processed

            for i in 0..threads {
                let ref_seqs_clone = Arc::clone(&shared_ref_seqs);
                let seq_chunk_clone = Arc::clone(&shared_seq_chunks[i]);
                let handle = thread::spawn(move || {
                    let mut chunk_results = Vec::new();
                    for (idx, seq) in seq_chunk_clone.iter().enumerate() {

                        if (idx + i) % progress_interval == 0 {
                            print!("\u{25A0}");
                            io::stdout().flush().unwrap();
                        }

                        let result = correct_seq(&ref_seqs_clone, seq, max_flex);
                        chunk_results.push(result);
                    }
                    chunk_results
                });
                handles.push(handle);
            }

            let mut results_vec = Vec::new();
            for handle in handles {
                let result = handle.join().unwrap();
                results_vec.push(result);
            }
            
            let mut all_results = Vec::new();
            for chunk in results_vec {
                all_results.extend(chunk);
            }
            println!(); // Move to the next line after progress bar
            return all_results;

        } else {
            return Vec::new(); // return an empty vector if threads is 0 or negative, which should not happen
        }
    }

    // #[pyfunction]
        fn correct_seq(ref_seqs: &Vec<String>, query_seq: &str, max_flex: isize) -> String {

        let query_len = query_seq.len() as isize;
        let query_count_A = query_seq.chars().filter(|&c| c == 'A').count() as isize;

        for ref_seq in ref_seqs {
            // Check if the length difference exceeds max_flex
            if (query_len - ref_seq.len() as isize).abs() > max_flex {
                continue; // Skip if the length difference exceeds max_flex
            }

            // quick check by counting A to make sure it does not exceed max_flex
            let ref_count_A = ref_seq.chars().filter(|&c| c == 'A').count() as isize;
            if (query_count_A - ref_count_A).abs() > max_flex {
                continue; // Skip if the count of 'A' differs by more than max_flex
            }

            if nw_align_bool(&ref_seq, query_seq, max_flex) {
                return ref_seq.clone(); // Return the first reference sequence that passes the checks
                }

        }
        "".to_string() // Return a empty string if no reference sequence matches
    }


    fn nw_align_bool(seq1: &str, seq2: &str, max_flex: isize) -> bool {
        const MATCH_SCORE: isize = 1;
        const MISMATCH_PENALTY: isize = 0;
        const GAP_PENALTY: isize = 0;

        let len1: isize = seq1.len() as isize;
        let len2: isize = seq2.len() as isize;

        // println!("Comparing sequences: {} and {}", seq1, seq2);

        if max_flex > min(len1, len2) {
            panic!("max_flex cannot be greater than the length of the sequences.");
        }

        if len1 < len2 {
            let (seq1, seq2) = (seq2, seq1); // Ensure seq1 is the shorter sequence
            let (len1, len2) = (len2, len1);
        }

        // Initialize the scoring matrix, seq1 is longer on top and seq2 is shorter on the left
        let mut matrix = vec![vec![0; (len1 + 1) as usize]; (len2 + 1) as usize];

        let seq1_vec: Vec<char> = seq1.chars().collect();
        let seq2_vec: Vec<char> = seq2.chars().collect();


        for i in 1..len1 + 1 {
            for j in 1..len2 + 1 {
                let idx = i as usize;
                let jdx = j as usize;

                let matched: bool = seq1_vec[idx - 1] == seq2_vec[jdx - 1];

                let score_diagonal = matrix[jdx - 1][idx - 1] + if matched {
                    MATCH_SCORE
                } else {
                    MISMATCH_PENALTY
                };

                let score_up = matrix[jdx - 1][idx] + GAP_PENALTY;
                let score_left = matrix[jdx][idx - 1] + GAP_PENALTY;

                matrix[jdx][idx] = max(score_diagonal, max(score_up, score_left));
            }
        }

        // for row in &matrix {
        //     println!("{:?}", row);
        // }

        if matrix[len2 as usize][len1 as usize] >= len1 - max_flex {
            return true; // Alignment score meets the threshold
        } else {
            return false; // Alignment score does not meet the threshold
        }

    }
}
