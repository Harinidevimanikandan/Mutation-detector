"""
AI-assisted Mutation Detection and Interpretation in Candida auris
--------------------------------------------------------------------
Upload a REFERENCE (wildtype) gene sequence and a QUERY (sample) sequence
in FASTA format. The app aligns them, finds DNA-level mutations, translates
to protein-level changes, checks them against a known-mutations table you
maintain, and gives a plain-English interpretation.

Run with:
    streamlit run app.py
"""

import streamlit as st
from Bio import SeqIO, Align
from Bio.Seq import Seq
import io

st.set_page_config(page_title="Candida auris Mutation Analyzer", layout="wide")

# ---------------------------------------------------------------------------
# EDIT THIS TABLE with real, literature-verified mutations for your organism.
# Format: protein position -> {wildtype AA, mutant AA, gene, significance}
# Fill this in from your actual source paper / NCBI / CDC data before your
# presentation. Leave empty entries as-is if you don't have verified data.
# ---------------------------------------------------------------------------
KNOWN_MUTATIONS = {
    # Example structure (replace with verified data):
    # 132: {"gene": "ERG11", "wt_aa": "Y", "mut_aa": "F", "significance": "Associated with azole resistance"},
}

# ---------------------------------------------------------------------------
def read_fasta(uploaded_file):
    """Parse an uploaded FASTA file and return (id, sequence string)."""
    text = uploaded_file.read().decode("utf-8")
    record = next(SeqIO.parse(io.StringIO(text), "fasta"))
    return record.id, str(record.seq).upper()


def align_sequences(ref_seq, query_seq):
    """Global pairwise alignment, returns aligned ref and query strings."""
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.match_score = 2
    aligner.mismatch_score = -1
    alignment = aligner.align(ref_seq, query_seq)[0]
    aligned_ref, aligned_query = alignment[0], alignment[1]
    return aligned_ref, aligned_query


def find_dna_mutations(aligned_ref, aligned_query):
    """Return list of (position_in_ref, ref_base, query_base)."""
    mutations = []
    ref_pos = 0
    for r_base, q_base in zip(aligned_ref, aligned_query):
        if r_base != "-":
            ref_pos += 1
        if r_base != q_base:
            mutations.append((ref_pos, r_base, q_base))
    return mutations


def translate_codon_changes(ref_seq, query_seq, dna_mutations):
    """
    For each DNA mutation, figure out which codon it falls in (assuming
    ref_seq starts at codon position 1, no frameshift) and report the
    amino acid change.
    """
    protein_changes = []
    seen_codons = set()
    for pos, r_base, q_base in dna_mutations:
        if r_base == "-" or q_base == "-":
            protein_changes.append({
                "dna_pos": pos, "type": "indel",
                "detail": f"Insertion/deletion near DNA position {pos} "
                           f"(ref='{r_base}', query='{q_base}') — may cause a frameshift."
            })
            continue

        codon_num = (pos - 1) // 3 + 1
        if codon_num in seen_codons:
            continue
        seen_codons.add(codon_num)

        start = (codon_num - 1) * 3
        ref_codon = ref_seq[start:start + 3]
        query_codon = query_seq[start:start + 3] if start + 3 <= len(query_seq) else None

        if len(ref_codon) == 3 and query_codon and len(query_codon) == 3:
            try:
                ref_aa = str(Seq(ref_codon).translate())
                query_aa = str(Seq(query_codon).translate())
            except Exception:
                continue
            if ref_aa != query_aa:
                protein_changes.append({
                    "dna_pos": pos, "type": "missense" if query_aa != "*" else "nonsense",
                    "codon_num": codon_num,
                    "ref_aa": ref_aa, "query_aa": query_aa,
                    "detail": f"Codon {codon_num}: {ref_aa} → {query_aa}"
                })
    return protein_changes


def interpret_mutation(change):
    """Rule-based plain-English interpretation. No AI API key required."""
    if change["type"] == "indel":
        return ("⚠️ This is an insertion/deletion. If it's not a multiple of 3 "
                "bases, it likely causes a frameshift, which usually has a major "
                "functional impact on the protein.")

    codon_num = change.get("codon_num")
    known = KNOWN_MUTATIONS.get(codon_num)
    if known:
        return (f"🔴 Matches a known mutation at codon {codon_num} "
                f"({known['gene']} {known['wt_aa']}{codon_num}{known['mut_aa']}): "
                f"{known['significance']}")

    if change["type"] == "nonsense":
        return ("⚠️ This introduces a premature stop codon (nonsense mutation), "
                "which typically truncates and disables the protein.")

    return ("ℹ️ This is a missense mutation not found in your known-mutations "
            "table. It changes the amino acid but its functional/clinical "
            "significance is unconfirmed — cross-check literature before "
            "reporting it as significant.")


# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
st.title("🧬 AI-assisted Mutation Detection — Candida auris")
st.caption("Upload a reference (wildtype) and a query (sample) FASTA sequence "
           "for the same gene to detect and interpret mutations.")

col1, col2 = st.columns(2)
with col1:
    ref_file = st.file_uploader("Reference (wildtype) FASTA", type=["fasta", "fa", "txt"])
with col2:
    query_file = st.file_uploader("Query (sample) FASTA", type=["fasta", "fa", "txt"])

if ref_file and query_file:
    ref_id, ref_seq = read_fasta(ref_file)
    query_id, query_seq = read_fasta(query_file)

    st.success(f"Loaded reference *{ref_id}* ({len(ref_seq)} bp) and "
               f"query *{query_id}* ({len(query_seq)} bp).")

    with st.spinner("Aligning sequences..."):
        aligned_ref, aligned_query = align_sequences(ref_seq, query_seq)
        dna_mutations = find_dna_mutations(aligned_ref, aligned_query)
        protein_changes = translate_codon_changes(ref_seq, query_seq, dna_mutations)

    st.subheader(f"🔎 Found {len(dna_mutations)} DNA-level difference(s)")
    if dna_mutations:
        st.dataframe(
            [{"Position": p, "Reference": r, "Query": q} for p, r, q in dna_mutations],
            use_container_width=True
        )
    else:
        st.info("No differences found — query matches reference exactly.")

    if protein_changes:
        st.subheader("🧪 Protein-level changes & interpretation")
        for change in protein_changes:
            st.markdown(f"*{change['detail']}*")
            st.write(interpret_mutation(change))
            st.divider()

    with st.expander("Show raw alignment"):
        st.code(f"REF:   {str(aligned_ref)}\nQUERY: {str(aligned_query)}")

else:
    st.info("⬆️ Upload both a reference and a query FASTA file to begin.")

st.divider()
st.caption(
    "Note: the known-mutations table in this app is a placeholder. Before "
    "presenting results, edit the KNOWN_MUTATIONS dictionary in app.py "
    "with verified positions/significance from your source literature "
    "(e.g., CDC or peer-reviewed papers on Candida auris resistance)."
)
