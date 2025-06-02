import streamlit as st
from Bio.Seq import Seq

from collections import defaultdict
import matplotlib.pyplot as plt
import re
import io

# --- Core Functions ---
def validate_dna(seq):
    return all(base in 'ATCG' for base in seq.upper())

def nucleotide_frequency(seq):
    return {base: seq.upper().count(base) for base in 'ATCG'}

def gc_content(seq):
    seq = seq.upper()
    return round((seq.count('G') + seq.count('C')) / len(seq) * 100, 2) if seq else 0.0

def transcribe(seq):
    return seq.upper().replace("T", "U")

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def translate_sequence(seq):
    try:
        return str(Seq(seq).translate(to_stop=True))
    except Exception as e:
        return f"Error during translation: {e}"

def codon_usage(seq):
    codon_table = defaultdict(int)
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3].upper()
        if len(codon) == 3 and all(base in 'ATCG' for base in codon):
            codon_table[codon] += 1
    return dict(codon_table)

def plot_codon_usage(codon_table):
    if not codon_table:
        st.warning("No codons to plot.")
        return
    sorted_codons = sorted(codon_table.keys())
    sorted_counts = [codon_table[c] for c in sorted_codons]
    plt.figure(figsize=(12, 5))
    plt.bar(sorted_codons, sorted_counts, color='skyblue')
    plt.xticks(rotation=90, fontsize=7)
    plt.title("Codon Usage")
    st.pyplot(plt)

def find_motif(seq, motif):
    pattern = re.compile(motif, re.IGNORECASE)
    return [match.start() for match in pattern.finditer(seq)]

def find_orfs(seq, min_length_bp=100):
    orfs = []
    for frame in range(3):
        current_seq_frame = Seq(seq[frame:].upper())
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

# --- Streamlit Layout ---
st.set_page_config(page_title="DNA Analyzer", layout="wide")
st.title("🧬 DNA Sequence Analyzer")

tab1, tab2 = st.tabs(["✍️ Manual Input", "📁 Upload FASTA File"])

with tab1:
    seq_input = st.text_area("Enter DNA Sequence (A, T, C, G only):", height=200)
    motif = st.text_input("Optional: Enter motif to search (e.g., ATG):")
    if st.button("Analyze Manual Sequence"):
        if not validate_dna(seq_input):
            st.error("Invalid DNA sequence. Only A, T, C, G are allowed.")
        else:
            st.success("✅ DNA Sequence is valid.")
            st.write(f"**Length:** {len(seq_input)} bp")
            st.write("**Nucleotide Frequencies:**", nucleotide_frequency(seq_input))
            st.write(f"**GC Content:** {gc_content(seq_input)}%")
            st.write("**RNA Transcription (first 100 bp):**", transcribe(seq_input[:100]))
            st.write("**Reverse Complement (first 100 bp):**", reverse_complement(seq_input[:100]))
            st.write("**Protein Translation (first 20 codons):**", translate_sequence(seq_input[:60]))

            codons = codon_usage(seq_input)
            plot_codon_usage(codons)

            if motif:
                positions = find_motif(seq_input, motif)
                st.write(f"Found **{len(positions)}** occurrences of motif `{motif}` at positions:", positions)

            orfs = find_orfs(seq_input)
            st.write(f"**ORFs longer than 100 bp found:** {len(orfs)}")
            for i, orf in enumerate(orfs[:3]):
                st.code(f"ORF {i+1} ({len(orf)*3} bp): {orf[:50]}...")

with tab2:
    uploaded_file = st.file_uploader("Upload a FASTA file", type=["fasta", "fa"])
    if uploaded_file:
        try:
            records = list(SeqIO.parse(io.StringIO(uploaded_file.getvalue().decode("utf-8")), "fasta"))
            if not records:
                st.error("No valid sequences found in file.")
            else:
                for record in records:
                    st.subheader(f"📄 {record.id}")
                    seq = str(record.seq)
                    st.write(f"**Length:** {len(seq)} bp")
                    st.write("**Nucleotide Frequencies:**", nucleotide_frequency(seq))
                    st.write(f"**GC Content:** {gc_content(seq)}%")
                    st.write("**RNA Transcription (first 100 bp):**", transcribe(seq[:100]))
                    st.write("**Reverse Complement (first 100 bp):**", reverse_complement(seq[:100]))
                    st.write("**Protein Translation (first 20 codons):**", translate_sequence(seq[:60]))
                    codons = codon_usage(seq)
                    plot_codon_usage(codons)
                    orfs = find_orfs(seq)
                    st.write(f"**ORFs longer than 100 bp found:** {len(orfs)}")
                    for i, orf in enumerate(orfs[:3]):
                        st.code(f"ORF {i+1} ({len(orf)*3} bp): {orf[:50]}...")
        except Exception as e:
            st.error(f"Error reading FASTA file: {e}")
