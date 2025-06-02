import streamlit as st
from Bio.Seq import Seq
from Bio import SeqIO
import matplotlib.pyplot as plt
import re
from collections import defaultdict
from io import StringIO

# --- DNA Analysis Functions ---
def validate_dna(seq):
    return all(base in 'ATCG' for base in seq.upper())

def nucleotide_frequency(seq):
    return {base: seq.upper().count(base) for base in 'ATCG'}

def gc_content(seq):
    seq = seq.upper()
    if not seq:
        return 0.0
    gc = seq.count('G') + seq.count('C')
    return round((gc / len(seq)) * 100, 2)

def transcribe(seq):
    return seq.upper().replace("T", "U")

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def translate_sequence(seq):
    try:
        return str(Seq(seq).translate(to_stop=True))
    except Exception as e:
        return f"Translation Error: {e}"

def codon_usage(seq):
    codon_table = defaultdict(int)
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3].upper()
        if len(codon) == 3 and all(base in 'ATCG' for base in codon):
            codon_table[codon] += 1
    return dict(codon_table)

def plot_codon_usage(codon_table, seq_id="Sequence"):
    if not codon_table:
        st.warning(f"No codons found in {seq_id}.")
        return

    sorted_codons = sorted(codon_table.keys())
    sorted_counts = [codon_table[c] for c in sorted_codons]

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.bar(sorted_codons, sorted_counts, color='skyblue')
    ax.set_title(f"Codon Usage - {seq_id}")
    ax.set_xlabel("Codons")
    ax.set_ylabel("Frequency")
    plt.xticks(rotation=90)
    st.pyplot(fig)

def find_motif(seq, motif):
    pattern = re.compile(motif, re.IGNORECASE)
    return [match.start() for match in pattern.finditer(seq)]

def find_orfs(seq, min_length_bp=100):
    orfs = []
    for strand_seq in [seq, reverse_complement(seq)]:
        for frame in range(3):
            current_seq = Seq(strand_seq[frame:].upper())
            trans = str(current_seq.translate(to_stop=False))
            start = -1
            for i, aa in enumerate(trans):
                if aa == 'M' and start == -1:
                    start = i
                elif aa == '*' and start != -1:
                    orf = trans[start:i]
                    if len(orf) * 3 >= min_length_bp:
                        orfs.append(orf)
                    start = -1
    return orfs

# --- Streamlit Interface ---
st.title("🧬 DNA Sequence Analyzer")
st.subheader("By Mahesh Prasanth Govind")
st.markdown("Upload a **FASTA** file or paste a DNA sequence manually to begin analysis.")

input_mode = st.radio("Input method:", ["Manual", "Upload FASTA"])

if input_mode == "Manual":
    seq_input = st.text_area("Enter your DNA sequence here (A, T, C, G):", height=200)
    seq_id = "Manual Input"

    if seq_input:
        if not validate_dna(seq_input):
            st.warning("The input contains non-DNA characters. Results may be unreliable.")
        
        st.subheader("Basic Analysis")
        st.write("**Length:**", len(seq_input))
        st.write("**Nucleotide Frequency:**", nucleotide_frequency(seq_input))
        st.write("**GC Content:**", gc_content(seq_input), "%")
        st.write("**Transcription:**", transcribe(seq_input[:100]), "...")
        st.write("**Reverse Complement:**", reverse_complement(seq_input[:100]), "...")
        st.write("**Translation (first 20 codons):**", translate_sequence(seq_input[:60]))

        st.subheader("Codon Usage")
        codons = codon_usage(seq_input)
        st.json(dict(list(codons.items())[:10]))  # Show first 10 codons
        plot_codon_usage(codons, seq_id)

        st.subheader("Start Codons (ATG)")
        atgs = find_motif(seq_input, "ATG")
        st.write(f"Found {len(atgs)} ATG motifs.")

        st.subheader("Open Reading Frames (ORFs)")
        orfs = find_orfs(seq_input)
        st.write(f"Found {len(orfs)} ORFs longer than 100 bp.")
        for i, orf in enumerate(orfs[:3]):
            st.code(f"ORF {i+1} (length: {len(orf)*3} bp): {orf[:50]}...", language="text")

elif input_mode == "Upload FASTA":
    fasta_file = st.file_uploader("Upload a FASTA file", type=["fasta", "fa", "txt"])
    if fasta_file is not None:
        try:
            contents = fasta_file.read().decode("utf-8")
            records = list(SeqIO.parse(StringIO(contents), "fasta"))

            if not records:
                st.error("No sequences found in the uploaded FASTA file.")
            else:
                for record in records:
                    seq_str = str(record.seq)
                    seq_id = record.id
                    st.markdown(f"## 🔍 Analysis for: `{seq_id}`")
                    st.write("**Length:**", len(seq_str))
                    st.write("**Nucleotide Frequency:**", nucleotide_frequency(seq_str))
                    st.write("**GC Content:**", gc_content(seq_str), "%")
                    st.write("**Transcription:**", transcribe(seq_str[:100]), "...")
                    st.write("**Reverse Complement:**", reverse_complement(seq_str[:100]), "...")
                    st.write("**Translation (first 20 codons):**", translate_sequence(seq_str[:60]))

                    st.subheader("Codon Usage")
                    codons = codon_usage(seq_str)
                    st.json(dict(list(codons.items())[:10]))
                    plot_codon_usage(codons, seq_id)

                    st.subheader("Start Codons (ATG)")
                    atgs = find_motif(seq_str, "ATG")
                    st.write(f"Found {len(atgs)} ATG motifs.")

                    st.subheader("Open Reading Frames (ORFs)")
                    orfs = find_orfs(seq_str)
                    st.write(f"Found {len(orfs)} ORFs longer than 100 bp.")
                    for i, orf in enumerate(orfs[:3]):
                        st.code(f"ORF {i+1} (length: {len(orf)*3} bp): {orf[:50]}...", language="text")
        except Exception as e:
            st.error(f"Error processing FASTA file: {e}")
