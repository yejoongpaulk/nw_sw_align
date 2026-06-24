# Needleman-Wunsch, Smith-Waterman Alignment Program

Optimized functions for calculating an optimal alignment between two genetic sequences.

## Usage & CLI Reference

This project includes a built-in command-line interface powered by `clap`. You can pass parameters using explicit flags or accept the calibrated biological default scores.

### Global Commands
To check configuration options or the release version from any terminal state:
```bash
# Display the interactive help manual
cargo run --release -- --help

# Check the engine version
cargo run --release -- --version
```

### 1. Global Sequence Alignment (`global`)
Finds the optimal alignment across the entire length of two explicit sequence strings using a high-performance linear gap implementation of the Needleman-Wunsch algorithm.

#### Options & Flags
*   `--seq1 <STRING>` — **(Required)** First sequence string (e.g. `ATCG`)
*   `--seq2 <STRING>` — **(Required)** Second sequence string (e.g. `ATGC`)
*   `--match-score <INT>` — Score added for a character match *[Default: 1]*
*   `--mismatch-penalty <INT>` — Penalty subtracted for a mismatch *[Default: -1]*
*   `--gap-penalty <INT>` — Penalty subtracted for introducing a gap *[Default: -2]*

#### Example Execution
```bash
cargo run --release -- global --seq1 ATCGGTA --seq2 ATGCGTA --match-score 2 --gap-penalty -3
```

---

### 2. Local Sequence Alignment (`local`)
Isolates and extracts the highest-scoring local sequence patches between two distinct string strands using an optimized Smith-Waterman variant.

#### Options & Flags
*   Shares identical string inputs (`--seq1`, `--seq2`) and scoring flag modifiers (`--match-score`, `--mismatch-penalty`, `--gap-penalty`) with the global alignment engine.

#### Example Execution
```bash
cargo run --release -- local --seq1 CCCCCATCGGGGG --seq2 AAAAATCGTTTTT
```

---

### 3. Multi-Tier Performance Profiling (`compare`)
Generates two completely randomized sequence strands containing biological characters (`A`, `C`, `T`, `G`) to evaluate raw engine throughput. It tracks execution times side-by-side across three structural layers: **Extreme Raw Pointer Optimization**, **Standard Optimization**, and the **Naive Array Pass**.

#### Options & Flags
*   `--size <SIZE>` — Total length of the generated random test strings *[Default: 50000]*
*   `--local` — Switches the benchmarking target from Global to Local alignment tracking

> **Memory Threshold Note:** If `--size` is configured above `30000`, the engine automatically bypasses the Naive Pass calculation to prevent memory reallocation lockups.

#### Example Executions
```bash
# Benchmark global alignment engines on a 15,000 character matrix
cargo run --release -- compare --size 15000

# Benchmark local alignment engines on the maximum default matrix size (50k)
cargo run --release -- compare --local
```

## Performance Comparison

| Sequence Length | Total DP Matrix Cells | Naive Approach (Safe Rust Only) | Standard Optimized (Raw Pointers) | Extreme Optimized (Raw Pointers, Raw Allocations in Intermediate Calculations) | Performance Speedup, Naive to Extreme |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **10,000 × 10,000** | 100 Million | 453.01 ms | 192.40 ms | **184.87 ms** | **2.45x Faster** |
| **20,000 × 20,000** | **400 Million** | 1.77 s | 790.72 ms | **666.29 ms** | **2.67x Faster** |
