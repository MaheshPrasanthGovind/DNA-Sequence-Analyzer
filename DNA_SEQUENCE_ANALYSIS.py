import streamlit as st
from Bio.Seq import Seq
from Bio import SeqIO
import matplotlib.pyplot as plt
import re
from collections import defaultdict
from io import StringIO
import json
import csv
import base64

# --- DNA Analysis Functions ---
def validate_dna(seq):
    # Returns True if seq only has A,T,C,G and ambiguous bases (N,Y,R,etc)
    valid_bases = set('ATCGNRYWSMKHBVD')  # IUPAC codes for ambiguous bases included
    return all(base in valid_bases for base in seq.upper())

def sanitize_dna(seq):
    valid_bases = set('ATCGNRYWSMKHBVD')
    seq = seq.upper()
    removed = 0
    sanitized = []
    for base in seq:
        if base in valid_bases:
            sanitized.append(base)
        else:
            removed += 1
    return ''.join(sanitized), removed

def nucleotide_frequency(seq):
    bases = 'ATCG'
    freq = {base: seq.count(base) for base in bases}
    # Count ambiguous bases separately
    ambiguous = len(seq) - sum(freq.values())
    freq['Ambiguous'] = ambiguous
    return freq

def gc_content(seq):
    if not seq:
        return 0.0
    gc = seq.count('G') + seq.count('C')
    return round((gc / len(seq)) * 100, 2)

def transcribe(seq):
    return seq.replace("T", "U")

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
        codon = seq[i:i+3]
        if len(codon) == 3:
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
    """Find ORFs in all 6 frames. Returns list of dicts with keys: orf_seq, strand, frame, start_nt, end_nt"""
    orfs = []
    seq = seq.upper()
    rev_comp_seq = reverse_complement(seq)
    seq_len = len(seq)
    start_codon = 'ATG'
    stop_codons = {'TAA', 'TAG', 'TGA'}

    # Search forward strand frames +1, +2, +3
    for frame in range(3):
        trans = str(Seq(seq[frame:]).translate(to_stop=False))
        aa_start = None
        for i, aa in enumerate(trans):
            if aa == 'M' and aa_start is None:
                aa_start = i
            elif aa == '*' and aa_start is not None:
                orf_length = i - aa_start
                if orf_length * 3 >= min_length_bp:
                    start_nt = frame + aa_start * 3
                    end_nt = frame + i * 3 + 3  # +3 because stop codon included
                    orf_seq = seq[start_nt:end_nt]
                    orfs.append({
                        'orf_seq': orf_seq,
                        'strand': '+',
                        'frame': frame + 1,
                        'start_nt': start_nt + 1,  # 1-based indexing
                        'end_nt': end_nt,
                        'aa_seq': trans[aa_start:i]
                    })
                aa_start = None
        # Handle ORF without stop at end
        if aa_start is not None:
            end_nt = seq_len
            orf_seq = seq[frame + aa_start * 3:end_nt]
            if len(orf_seq) >= min_length_bp:
                orfs.append({
                    'orf_seq': orf_seq,
                    'strand': '+',
                    'frame': frame + 1,
                    'start_nt': frame + aa_start * 3 + 1,
                    'end_nt': end_nt,
                    'aa_seq': trans[aa_start:]
                })

    # Search reverse strand frames -1, -2, -3
    for frame in range(3):
        trans = str(Seq(rev_comp_seq[frame:]).translate(to_stop=False))
        aa_start = None
        for i, aa in enumerate(trans):
            if aa == 'M' and aa_start is None:
                aa_start = i
            elif aa == '*' and aa_start is not None:
                orf_length = i - aa_start
                if orf_length * 3 >= min_length_bp:
                    # Calculate coordinates relative to forward strand
                    end_nt = seq_len - (frame + aa_start * 3)
                    start_nt = seq_len - (frame + i * 3 + 3) + 1
                    orf_seq = rev_comp_seq[frame + aa_start * 3:frame + i * 3 + 3]
                    orfs.append({
                        'orf_seq': orf_seq,
                        'strand': '-',
                        'frame': -(frame + 1),
                        'start_nt': start_nt,
                        'end_nt': end_nt,
                        'aa_seq': trans[aa_start:i]
                    })
                aa_start = None
        # Handle ORF without stop at end
        if aa_start is not None:
            end_nt = 1
            start_nt = seq_len - (frame + aa_start * 3)
            orf_seq = rev_comp_seq[frame + aa_start * 3:]
            if len(orf_seq) >= min_length_bp:
                orfs.append({
                    'orf_seq': orf_seq,
                    'strand': '-',
                    'frame': -(frame + 1),
                    'start_nt': start_nt,
                    'end_nt': end_nt,
                    'aa_seq': trans[aa_start:]
                })
    return orfs

def translate_six_frames(seq):
    """Translate all 6 frames and return dict with frame name keys and aa seqs"""
    frames = {}
    seq = seq.upper()
    rev_comp = reverse_complement(seq)
    for i in range(3):
        frames[f"+{i+1}"] = str(Seq(seq[i:]).translate(to_stop=False))
        frames[f"-{i+1}"] = str(Seq(rev_comp[i:]).translate(to_stop=False))
    return frames

# --- Export helpers ---
def get_download_link(content, filename, mime):
    b64 = base64.b64encode(content.encode()).decode()
    return f'<a href="data:{mime};base64,{b64}" download="{filename}">Download {filename}</a>'

def orfs_to_fasta(orfs):
    fasta = []
    for i, orf in enumerate(orfs, 1):
        fasta.append(f">ORF_{i} strand={orf['strand']} frame={orf['frame']} start={orf['start_nt']} end={orf['end_nt']}\n{orf['orf_seq']}")
    return '\n'.join(fasta)

def aa_to_fasta(orfs):
    fasta = []
    for i, orf in enumerate(orfs, 1):
        fasta.append(f">ORF_{i} strand={orf['strand']} frame={orf['frame']} start={orf['start_nt']} end={orf['end_nt']}\n{orf['aa_seq']}")
    return '\n'.join(fasta)

# --- Streamlit Interface ---
st.set_page_config(page_title="DNA Sequence Analyzer", layout="wide")

st.markdown("""
    <style>
    html, body, [class*="css"] {
        background-color: #000000 !important;
        color: #f0f0f0 !important;
    }
    textarea, input {
        background-color: #333333 !important;
        color: #f0f0f0 !important;
    }
    .stCodeBlock, .stMarkdown {
        background-color: #2d2d2d !important;
        color: #f0f0f0 !important;
    }
    </style>
""", unsafe_allow_html=True)

st.title("🧬 DNA Sequence Analyzer")
st.subheader("By Mahesh Prasanth Govind")

st.markdown("Upload a **FASTA** file or paste a DNA sequence manually to begin analysis.")

input_mode = st.radio("Input method:", ["Manual", "Upload FASTA"])

show_6frames = st.checkbox("Show full translation in all 6 frames")

min_orf_len = st.number_input("Minimum ORF length (bp)", min_value=30, max_value=10000, value=100, step=10)

if input_mode == "Manual":
    raw_seq = st.text_area("Enter your DNA sequence here (A, T, C, G, ambiguous bases allowed):", height=200)
    seq_id = "Manual Input"
    seq_input, removed_count = sanitize_dna(raw_seq)

    if raw_seq:
        if removed_count > 0:
            st.warning(f"Removed {removed_count} invalid bases from input.")

        if not validate_dna(raw_seq):
            st.warning("Warning: Input contained ambiguous or invalid bases.")

        st.subheader("Basic Analysis")
        st.write("**Length:**", len(seq_input))
        st.write("**Nucleotide Frequency:**", nucleotide_frequency(seq_input))
        st.write("**GC Content:**", gc_content(seq_input), "%")
        st.write("**Transcription:**", transcribe(seq_input[:100]), "...")
        st.write("**Reverse Complement:**", reverse_complement(seq_input[:100]), "...")
        st.write("**Translation (first 20 codons):**", translate_sequence(seq_input[:60]))

        if show_6frames:
            st.subheader("Full 6-frame Translation")
            frames = translate_six_frames(seq_input)
            for frame_name, aa_seq in frames.items():
                st.markdown(f"**Frame {frame_name}**")
                st.text(aa_seq[:300] + " ...")

        st.subheader("Codon Usage")
        if len(seq_input) >= 6:
            codons = codon_usage(seq_input)
            st.json(dict(list(codons.items())[:10]))  # Show first 10
            plot_codon_usage(codons, seq_id)
        else:
            st.warning("Sequence too short for codon usage analysis.")

        st.subheader("Start Codons (ATG)")
        atgs = find_motif(seq_input, "ATG")
        st.write(f"Found {len(atgs)} ATG motifs.")

        st.subheader("Open Reading Frames (ORFs)")
        orfs = find_orfs(seq_input, min_length_bp=min_orf_len)
        st.write(f"Found {len(orfs)} ORFs longer than {min_orf_len} bp.")

        for i, orf in enumerate(orfs[:5]):
            st.code(f"ORF {i+1} (strand: {orf['strand']} frame: {orf['frame']} start: {orf['start_nt']} end: {orf['end_nt']} length: {len(orf['orf_seq'])} bp)\nAA seq: {orf['aa_seq'][:50]}...", language="text")

        # Export ORFs and AA seqs as FASTA
        if orfs:
            fasta_orfs = orfs_to_fasta(orfs)
            fasta_aa = aa_to_fasta(orfs)
            st.markdown(get_download_link(fasta_orfs, f"{seq_id}_orfs.fasta", "text/plain"), unsafe_allow_html=True)
            st.markdown(get_download_link(fasta_aa, f"{seq_id}_orfs_aa.fasta", "text/plain"), unsafe_allow_html=True)

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
                    raw_seq = str(record.seq)
                    seq_id = record.id
                    seq_str, removed_count = sanitize_dna(raw_seq)

                    st.markdown(f"## 🔍 Analysis for: `{seq_id}`")

                    if removed_count > 0:
                        st.warning(f"Removed {removed_count} invalid bases from input.")

                    if not validate_dna(raw_seq):
                        st.warning("Warning: Input contained ambiguous or invalid bases.")

                    st.write("**Length:**", len(seq_str))
                    st.write("**Nucleotide Frequency:**", nucleotide_frequency(seq_str))
                    st.write("**GC Content:**", gc_content(seq_str), "%")
                    st.write("**Transcription:**", transcribe(seq_str[:100]), "...")
                    st.write("**Reverse Complement:**", reverse_complement(seq_str[:100]), "...")
                    st.write("**Translation (first 20 codons):**", translate_sequence(seq_str[:60]))

                    if show_6frames:
                        st.subheader("Full 6-frame Translation")
                        frames = translate_six_frames(seq_str)
                        for frame_name, aa_seq in frames.items():
                            st.markdown(f"**Frame {frame_name}**")
                            st.text(aa_seq[:300] + " ...")

                    st.subheader("Codon Usage")
                    if len(seq_str) >= 6:
                        codons = codon_usage(seq_str)
                        st.json(dict(list(codons.items())[:10]))
                        plot_codon_usage(codons, seq_id)
                    else:
                        st.warning("Sequence too short for codon usage analysis.")

                    st.subheader("Start Codons (ATG)")
                    atgs = find_motif(seq_str, "ATG")
                    st.write(f"Found {len(atgs)} ATG motifs.")

                    st.subheader("Open Reading Frames (ORFs)")
                    orfs = find_orfs(seq_str, min_length_bp=min_orf_len)
                    st.write(f"Found {len(orfs)} ORFs longer than {min_orf_len} bp.")
                    for i, orf in enumerate(orfs[:5]):
                        st.code(f"ORF {i+1} (strand: {orf['strand']} frame: {orf['frame']} start: {orf['start_nt']} end: {orf['end_nt']} length: {len(orf['orf_seq'])} bp)\nAA seq: {orf['aa_seq'][:50]}...", language="text")

                    if orfs:
                        fasta_orfs = orfs_to_fasta(orfs)
                        fasta_aa = aa_to_fasta(orfs)
                        st.markdown(get_download_link(fasta_orfs, f"{seq_id}_orfs.fasta", "text/plain"), unsafe_allow_html=True)
                        st.markdown(get_download_link(fasta_aa, f"{seq_id}_orfs_aa.fasta", "text/plain"), unsafe_allow_html=True)

        except Exception as e:
            st.error(f"Error processing FASTA file: {e}")
