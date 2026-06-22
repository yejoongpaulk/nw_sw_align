# Needleman-Wunsch, Smith-Waterman Alignment

Optimized functions for calculating the dynamic programming matrices for alignments.


| Sequence Length | Total DP Matrix Cells | Naive Approach (Safe Rust Only) | Standard Optimized (Raw Pointers) | Extreme Optimized (Raw Pointers, Raw Vector) | Performance Speedup |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **10,000 × 10,000** | 100 Million | 453.01 ms | 192.40 ms | **184.87 ms** | **2.45x Faster** |
| **20,000 × 20,000** | **400 Million** | 1.77 s | 790.72 ms | **666.29 ms** | **2.67x Faster** |