# Local Sequence Alignment Tool (Affine Gap Smith–Waterman)

## Overview
This application performs **local alignment of two nucleotide sequences** using the **Smith–Waterman algorithm with affine gap penalties**.  
It reads sequences from FASTA files and scoring parameters from a CSV file, then outputs the optimal local alignment and score.

This implementation supports:
- custom substitution matrices
- configurable gap penalties
- biologically accurate affine gap scoring
- file-based input/output

---

## Features
- Correct Smith–Waterman local alignment
- Affine gap penalty model (k−1 convention)
- Custom nucleotide scoring matrix
- CSV-based scoring configuration
- FASTA input support
- Alignment output saved to file
- Fully deterministic traceback

---

## Input Files

### 1. FASTA Files
Two FASTA files must exist in the same directory:

