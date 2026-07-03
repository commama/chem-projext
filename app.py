import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from st_py3dmol import show_struct

# --- Page Configuration ---
st.set_page_config(page_title="PharmDesign Pro | Unit II", layout="wide")

# Custom CSS to hide "Made with Streamlit" and look like professional software
st.markdown("""
    <style>
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    .stApp {background-color: #f8f9fa;}
    .main-header {color: #1e3a8a; font-size: 2.5rem; font-weight: bold; margin-bottom: 20px;}
    </style>
    """, unsafe_allow_html=True)

# --- Sidebar Navigation ---
st.sidebar.title("🔬 Unit II: In Silico Design")
menu = st.sidebar.radio("Navigation", 
    ["Lipinski Checker", "Molecular Docking Demo", "ADMET Dashboard", "Case Study Hub"])

# --- Helper Functions ---
def calculate_lipinski(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        
        # Rule of 5 logic
        res = {
            "MW (<500 Da)": [mw, mw <= 500],
            "LogP (<5)": [logp, logp <= 5],
            "H-Donors (<5)": [hbd, hbd <= 5],
            "H-Acceptors (<10)": [hba, hba <= 10]
        }
        return res, mol
    return None, None

# --- Logic for Modules ---

if menu == "Lipinski Checker":
    st.markdown('<div class="main-header">Lipinski\'s Rule of Five Checker</div>', unsafe_allow_html=True)
    smiles_input = st.text_input("Enter SMILES String (e.g., Aspirin: CC(=O)OC1=CC=CC=C1C(=O)O)", "CC(=O)OC1=CC=CC=C1C(=O)O")
    
    data, mol = calculate_lipinski(smiles_input)
    if data:
        cols = st.columns(4)
        pass_count = 0
        for i, (prop, val) in enumerate(data.items()):
            with cols[i]:
                status = "✅ Pass" if val[1] else "❌ Fail"
                st.metric(prop, f"{val[0]:.2f}", status)
                if val[1]: pass_count += 1
        
        if pass_count >= 3:
            st.success("Verdict: Drug-Likely (Compliant)")
        else:
            st.error("Verdict: Poor Drug-Likeness")

elif menu == "Molecular Docking Demo":
    st.markdown('<div class="main-header">Molecular Docking Simulator</div>', unsafe_allow_header=True)
    st.info("Visualizing Ligand-Protein Binding Affinity at the Active Site")
    
    # Simplified visual demo using Py3Dmol
    # Displaying a sample protein-ligand pocket
    view = show_struct(pdb='1vsn', height=500) # Sample PDB
    st.write("Target: HIV-1 Protease | Ligand: Indinavir")

elif menu == "ADMET Dashboard":
    st.markdown('<div class="main-header">ADMET Property Dashboard</div>', unsafe_allow_html=True)
    
    smiles = st.text_input("Ligand SMILES", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C") # Caffeine
    
    # Mock prediction logic (Calculated via RDKit descriptors)
    mol = Chem.MolFromSmiles(smiles)
    tpsa = Descriptors.TPSA(mol)
    rotb = Lipinski.NumRotatableBonds(mol)
    
    fig = go.Figure(data=go.Scatterpolar(
      r=[tpsa, rotb, 50, 40, 60], # Mocked values
      theta=['TPSA','Rotatable Bonds','Solubility','BBB Permeability','GI Absorption'],
      fill='toself',
      line_color='#1e3a8a'
    ))
    fig.update_layout(polar=dict(radialaxis=dict(visible=True, range=[0, 150])), showlegend=False)
    
    st.plotly_chart(fig)
    st.caption("Radar chart showing pharmacokinetic profile balance.")

elif menu == "Case Study Hub":
    st.markdown('<div class="main-header">AI in Drug Discovery: Case Studies</div>', unsafe_allow_html=True)
    
    studies = {
        "DeepMind AlphaFold": "Predicting protein structures with atomic accuracy, accelerating target identification.",
        "Insilico Medicine": "First AI-discovered drug for Idiopathic Pulmonary Fibrosis to enter clinical trials.",
        "Exscientia": "Centaur Chemist platform used to design highly selective immuno-oncology ligands."
    }
    
    for title, desc in studies.items():
        with st.expander(title):
            st.write(desc)
            st.button(f"Read Full Paper: {title[:5]}")