"""
Program: chips_app.py
Author: Alexandre Prefontaine
Affiliation: University of Manitoba, Manitoba Centre for Proteomics and Systems Biology
Date of Last Update: 2024-06-20

Streamlit web application for CHIPs (Chromatographic Hydrophobicity Index for Peptides).

Run with:
    streamlit run chips_app.py
"""

import streamlit as st
import pandas as pd
from pathlib import Path
from CHIPS_core import calculate_HI, VALID_AA, KNOWN_HEADERS


def is_valid_sequence(value: str) -> bool:
    """
    Returns True if the value looks like a peptide sequence.
    Rejects known header strings and sequences containing invalid characters.
    """
    cleaned = value.upper().strip()
    if cleaned.lower() in KNOWN_HEADERS:
        return False
    return all(c in VALID_AA for c in cleaned)


# --- Page configuration ---
st.set_page_config(page_title="CHIPs", layout="wide")

# --- Title and description ---
st.title("CHIPs — Chromatographic Hydrophobicity Index for Peptides")
st.write(
    "Predict peptide hydrophobicity under reversed-phase HPLC conditions "
    "using an additive model with logarithmic length correction."
)

# Working Equation
st.latex(r"HI = (1 - \alpha \ln n)\left(\sum_{i=1}^{n} c_i + \beta\right)")

# Citation
st.caption(
    "Prefontaine, Krokhin. (2024) *CHIPs: Chromatographic Hydrophobicity Index for Peptides.* "
    "Journal of Proteome Research. [Link to publication](https://doi.org/XXXXXXX)"
)

st.divider()

# --- Phase selector ---
phase = st.radio("Select ion pairing agent", options=["FA", "TFA"], horizontal=True)

st.divider()

# --- Input method selector ---
input_mode = st.radio(
    "Input method",
    options=["Manual entry", "Upload CSV"],
    horizontal=True
)

uploaded_file = None

if input_mode == "Manual entry":
    raw_input = st.text_area(
        "Enter peptide sequences (one per line)",
        height=200,
        placeholder="AANLVR\nFVTDIDELGK\nVTFLGLQHWVPELAR"
    )
    sequences = [
        s.strip().upper() for s in raw_input.splitlines()
        if s.strip() and is_valid_sequence(s)
    ]

else:
    uploaded_file = st.file_uploader(
        "Upload a CSV file (single column of sequences)",
        type=["csv"]
    )
    sequences = []

    if uploaded_file is not None:
        df_in = pd.read_csv(uploaded_file, header=None)
        raw_sequences = df_in.iloc[:, 0].astype(str).str.strip().str.upper().tolist()
        sequences = [s for s in raw_sequences if is_valid_sequence(s)]

st.divider()

# --- Run button and results ---
if st.button("Predict HI", type="primary"):

    if not sequences:
        st.warning("No valid sequences found. Please enter or upload peptide sequences.")
    else:
        results = []
        skipped = []

        for seq in sequences:
            try:
                hi = calculate_HI(seq, phase)
                results.append({"Peptide": seq, "HI": hi})
            except Exception as e:
                skipped.append(seq)

        df_out = pd.DataFrame(results)

        st.success(f"Predicted HI for {len(results)} peptides using {phase} coefficients.")

        if skipped:
            st.warning(f"Skipped {len(skipped)} invalid sequences: {', '.join(skipped)}")

        st.dataframe(df_out, use_container_width=True, hide_index=True)

        # --- Download button ---
        if uploaded_file is not None:
            output_name = Path(uploaded_file.name).stem + "_CHIPs_HI.csv"
        else:
            output_name = "CHIPs_predicted_HI.csv"

        csv_bytes = df_out.to_csv(index=False).encode("utf-8")

        st.download_button(
            label="Download results as CSV",
            data=csv_bytes,
            file_name=output_name,
            mime="text/csv"
        )

# --- Footer ---
st.divider()

col1, col2 = st.columns([1, 3])

with col1:
    st.image("assets/UM-logo-horizontal-CMYK.jpg", width=120)

with col2:
    st.markdown(
        """
        **Alexandre Préfontaine** — Krokhin Laboratory  
        Manitoba Centre for Proteomics and Systems Biology  
        [University of Manitoba](https://umanitoba.ca) · Department of Internal Medicine  
        """
    )