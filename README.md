# DNA Sequence Analyzer

A Python tool for analyzing DNA sequences, including:

- Validating DNA sequences
- Calculating nucleotide frequencies and GC content
- Transcribing DNA to RNA
- Finding reverse complements
- Translating DNA to proteins
- Calculating codon usage
- Finding open reading frames (ORFs)
- Visualizing codon usage with plots

## Features

- Input DNA sequence manually or upload FASTA files
- Uses BioPython for reliable sequence handling
- Generates codon usage bar plots with matplotlib
- Finds ORFs across all six reading frames
- Handles errors gracefully and supports interactive analysis

## Installation

1. Clone this repository:

```bash
git clone https://github.com/MaheshPrasanthGovind/dna-sequence-analyzer.git
cd dna-sequence-analyzer

pip install biopython matplotlib
python dna_analyzer.py

from dna_analyzer import analyze_sequence

seq = "ATGCGTATCGTAGCTAGCTAGCTAA"
analyze_sequence(seq, seq_id="Example Sequence")

--- ANALYSIS FOR: Example Sequence ---
Sequence preview: ATGCGTATCGTAGCTAGCTAGCTAA... (Total length: 24 bp)
Valid DNA (A, T, C, G only): True
Nucleotide Frequencies: {'A': 6, 'T': 6, 'C': 6, 'G': 6}
GC Content: 50.0 %
Transcription (RNA preview): AUGCGUAUCGUAGCUAGCUAGCUAA ...
Reverse Complement (preview): TTAGCTAGCTAGCTACGATACGCAT ...
Protein Translation (first 20 codons): MRYVS*
Codon Usage (showing first 10 entries): {'ATG': 1, 'CGT': 1, ...}
Found 1 potential start codons (ATG).
Found 1 ORFs longer than 100 bp (considering all 6 frames).
Preview of found ORFs (first 3):
  ORF 1 (length: 24 bp): MRYVS

