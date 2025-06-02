# --- Setup for Google Colab ---
!pip install biopython matplotlib

# --- Imports ---
from Bio.Seq import Seq
from Bio import SeqIO
import re
import matplotlib.pyplot as plt
from collections import defaultdict
from google.colab import files
# --- CORE DNA ANALYSIS FUNCTIONS ---

def validate_dna(seq):
    """
    Checks if a given sequence contains only valid DNA nucleotides (A, T, C, G).
    Args:
        seq (str): The DNA sequence string.
    Returns:
        bool: True if valid DNA, False otherwise.
    """
    return all(base in 'ATCG' for base in seq.upper())

def nucleotide_frequency(seq):
    """
    Calculates the frequency of each nucleotide (A, T, C, G) in a DNA sequence.
    Args:
        seq (str): The DNA sequence string.
    Returns:
        dict: A dictionary with nucleotides as keys and their counts as values.
    """
    return {base: seq.upper().count(base) for base in 'ATCG'}

def gc_content(seq):
    """
    Calculates the GC (Guanine-Cytosine) content percentage of a DNA sequence.
    Args:
        seq (str): The DNA sequence string.
    Returns:
        float: The GC content as a percentage, rounded to two decimal places.
              Returns 0.0 if the sequence is empty to avoid ZeroDivisionError.
    """
    seq = seq.upper()
    if not seq: # Handle empty sequence to prevent division by zero
        return 0.0
    gc = seq.count('G') + seq.count('C')
    return round((gc / len(seq)) * 100, 2)

def transcribe(seq):
    """
    Transcribes a DNA sequence into an RNA sequence (replaces 'T' with 'U').
    Args:
        seq (str): The DNA sequence string.
    Returns:
        str: The RNA sequence string.
    """
    return seq.upper().replace("T", "U")

def reverse_complement(seq):
    """
    Generates the reverse complement of a DNA sequence.
    Args:
        seq (str): The DNA sequence string.
    Returns:
        str: The reverse complement sequence string.
    """
    return str(Seq(seq).reverse_complement())

def translate_sequence(seq):
    """
    Translates a DNA sequence into an amino acid sequence (protein).
    Translation stops at the first stop codon.
    Args:
        seq (str): The DNA sequence string.
    Returns:
        str: The translated protein sequence string.
    """
    try:
        return str(Seq(seq).translate(to_stop=True))
    except Exception as e:
        return f"Error during translation: {e}"

def codon_usage(seq):
    """
    Calculates the frequency of each codon (3-nucleotide triplet) in a DNA sequence.
    Args:
        seq (str): The DNA sequence string.
    Returns:
        dict: A dictionary with codons as keys and their counts as values.
    """
    codon_table = defaultdict(int) # Use defaultdict to easily increment counts for new codons
    # Iterate through the sequence in steps of 3 to get codons
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3].upper()
        # Basic validation for codon bases
        if len(codon) == 3 and all(base in 'ATCG' for base in codon):
            codon_table[codon] += 1
    return dict(codon_table) # Convert back to a regular dict for final output

def plot_codon_usage(codon_table, seq_id="Sequence"):
    """
    Generates and displays a bar plot of codon usage frequencies.
    Args:
        codon_table (dict): A dictionary of codon frequencies.
        seq_id (str): Identifier for the sequence, used in the plot title.
    """
    if not codon_table: # Handle empty codon table
        print(f"No codons to plot for {seq_id}.")
        return

    codons = list(codon_table.keys())
    counts = list(codon_table.values())

    # Sort codons for consistent plotting order
    sorted_codons = sorted(codons)
    sorted_counts = [codon_table[c] for c in sorted_codons]

    plt.figure(figsize=(15, 7))
    plt.bar(sorted_codons, sorted_counts, color='skyblue')
    plt.xlabel("Codons")
    plt.ylabel("Frequency")
    plt.title(f"Codon Usage for {seq_id}")
    plt.xticks(rotation=90, fontsize=8)
    plt.tight_layout()
    plt.show()

def find_motif(seq, motif):
    """
    Finds all occurrences of a given DNA motif within a sequence and returns their start positions.
    Args:
        seq (str): The DNA sequence string to search within.
        motif (str): The DNA motif to search for.
    Returns:
        list: A list of integer starting positions where the motif was found.
    """
    pattern = re.compile(motif, re.IGNORECASE)
    return [match.start() for match in pattern.finditer(seq)]

def find_orfs(seq, min_length_bp=100):
    """
    Finds potential Open Reading Frames (ORFs) within a DNA sequence across 6 frames.
    This version looks for regions starting with ATG and ending with a stop codon (*).
    Args:
        seq (str): The DNA sequence string.
        min_length_bp (int): Minimum length of the ORF in base pairs.
    Returns:
        list: A list of translated protein sequences for the found ORFs.
    """
    orfs = []
    # Define stop codons for DNA
    stop_codons_dna = ["TAA", "TAG", "TGA"]

    # Check all three forward reading frames
    for frame in range(3):
        current_seq_frame = Seq(seq[frame:].upper())
        trans = str(current_seq_frame.translate(to_stop=False)) # Translate without stopping for full protein

        # Iterate through the translated protein to find ORFs
        current_orf_start_aa = -1 # Start index in amino acids
        for i, aa in enumerate(trans):
            if aa == 'M' and current_orf_start_aa == -1: # Found a Methionine (start codon)
                current_orf_start_aa = i
            elif aa == '*' and current_orf_start_aa != -1: # Found a stop codon after a start
                orf_protein = trans[current_orf_start_aa:i]
                # Convert protein length back to base pairs for min_length_bp check
                if len(orf_protein) * 3 >= min_length_bp:
                    orfs.append(orf_protein)
                current_orf_start_aa = -1 # Reset for next ORF

    # Check all three reverse reading frames (using reverse complement)
    rev_comp_seq = reverse_complement(seq)
    for frame in range(3):
        current_seq_frame = Seq(rev_comp_seq[frame:].upper())
        trans = str(current_seq_frame.translate(to_stop=False))

        current_orf_start_aa = -1
        for i, aa in enumerate(trans):
            if aa == 'M' and current_orf_start_aa == -1:
                current_orf_start_aa = i
            elif aa == '*' and current_orf_start_aa != -1:
                orf_protein = trans[current_orf_start_aa:i]
                if len(orf_protein) * 3 >= min_length_bp:
                    orfs.append(orf_protein)
                current_orf_start_aa = -1

    return orfs

# --- Sequence Analysis Display ---

def analyze_sequence(seq_str, seq_id="Manual Input"):
    """
    Performs and prints various analyses on a single DNA sequence.
    Args:
        seq_str (str): The DNA sequence string to analyze.
        seq_id (str): An identifier for the sequence (e.g., header from FASTA).
    """
    print(f"\n--- ANALYSIS FOR: {seq_id} ---")
    print(f"Sequence preview: {seq_str[:100]}... (Total length: {len(seq_str)} bp)") # Print more of the sequence

    if not validate_dna(seq_str):
        print("Warning: Sequence contains non-DNA characters. Results may be inaccurate.")

    print("Valid DNA (A, T, C, G only):", validate_dna(seq_str))
    print("Nucleotide Frequencies:", nucleotide_frequency(seq_str))
    print("GC Content: ", gc_content(seq_str), "%")
    print("Transcription (RNA preview):", transcribe(seq_str[:100]), "...")
    print("Reverse Complement (preview):", reverse_complement(seq_str[:100]), "...")
    print("Protein Translation (first 20 codons):", translate_sequence(seq_str[:60])) # Translate first 60 bp (20 codons)

    codons = codon_usage(seq_str)
    # Print only the first 10 entries of codon usage for brevity
    print("Codon Usage (showing first 10 entries):", dict(list(codons.items())[:10]))
    plot_codon_usage(codons, seq_id=seq_id) # Show the codon usage plot

    start_codons_atg = find_motif(seq_str, "ATG") # Find all 'ATG' start codons
    print(f"Found {len(start_codons_atg)} potential start codons (ATG).")

    orfs = find_orfs(seq_str) # Find ORFs
    print(f"Found {len(orfs)} ORFs longer than 100 bp (considering all 6 frames).")
    if orfs:
        # Show a preview of some found ORFs
        print("Preview of found ORFs (first 3):")
        for i, orf in enumerate(orfs[:3]):
            print(f"  ORF {i+1} (length: {len(orf)*3} bp): {orf[:50]}...") # Show protein preview


# --- MAIN EXECUTION ---

def main():
    """
    Main function to run the DNA Sequence Analyzer.
    Prompts the user to enter DNA manually or upload a FASTA file.
    """
    print("=== DNA SEQUENCE ANALYZER (Colab Ready) ===")

    # Loop until valid input is received
    while True:
        choice = input("Enter DNA manually or upload from FASTA? (manual/fasta): ").strip().lower()

        if choice == "manual":
            user_seq = input("Enter DNA sequence: ").strip()
            if validate_dna(user_seq):
                analyze_sequence(user_seq, seq_id="Manual Input")
            else:
                print("Invalid DNA sequence. Please enter only A, T, C, G characters.")
            break # Exit loop after processing manual input

        elif choice == "fasta":
            print("\n📁 Please upload your FASTA file.")
            # This line opens the Colab file picker
            uploaded = files.upload()

            if not uploaded:
                print("No file uploaded. Please try again.")
                continue # Ask for choice again

            # Get the filename (assuming only one file is uploaded)
            fasta_filename = list(uploaded.keys())[0]
            print(f"File '{fasta_filename}' uploaded successfully.")

            try:
                # Use BioPython's SeqIO to parse the FASTA file
                # This is more robust than a custom parser
                records = list(SeqIO.parse(fasta_filename, "fasta"))

                if not records:
                    print("No sequences found in the FASTA file. Please check file format.")
                    # Optionally, you might want to remove the uploaded file here
                    # os.remove(fasta_filename) # Requires 'import os'
                    break

                for record in records:
                    # Convert Seq object to string for analysis functions
                    analyze_sequence(str(record.seq), seq_id=record.id)

            except FileNotFoundError: # Although 'files.upload()' handles this mostly
                print(f"Error: File '{fasta_filename}' not found. This shouldn't happen after upload.")
            except Exception as e: # Catch other potential errors during file processing
                print(f"An error occurred while processing the FASTA file: {e}")

            # Remove the uploaded file to clean up the Colab environment
            # This requires 'import os' at the top if you want to enable it
            # import os
            # os.remove(fasta_filename)
            break # Exit loop after processing FASTA file

        else:
            print("Invalid input. Please type 'manual' or 'fasta'.")

# This ensures main() runs only when the script is executed directly
if __name__ == "__main__":
    main()
