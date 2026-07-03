import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, AllChem
import pandas as pd
import plotly.graph_objects as go
import st_py3dmol

# --- Page Configuration ---
st.set_page_config(page_title="PharmDesign Pro", layout="wide")

st.markdown("""
    <style>
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    .stApp {background-color: #fcfcfc;}
    .main-header {color: #1e3a8a; font-size: 2.2rem; font-weight: bold; margin-bottom: 20px; border-bottom: 2px solid #1e3a8a;}
    </style>
    """, unsafe_allow_html=True)

# --- Sidebar ---
st.sidebar.title("🔬 In Silico Drug Design")
menu = st.sidebar.radio("Analysis Modules", 
    ["Lipinski Checker", "Molecular Docking Demo", "ADMET Dashboard", "Case Study Hub"])

# --- Helper Functions ---
def get_lipinski_data(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        res = {
            "Molecular Weight": [mw, mw <= 500, 500],
            "LogP": [logp, logp <= 5, 5],
            "H-Donors": [hbd, hbd <= 5, 5],
            "H-Acceptors": [hba, hba <= 10, 10]
        }
        return res, mol
    return None, None

# --- Modules ---
if menu == "Lipinski Checker":
    st.markdown('<div class="main-header">Lipinski\'s Rule of Five Checker</div>', unsafe_allow_html=True)
    smiles_input = st.text_input("Enter SMILES String", "CC(=O)OC1=CC=CC=C1C(=O)O") # Aspirin
    
    data, mol = get_lipinski_data(smiles_input)
    if data:
        cols = st.columns(4)
        passes = 0
        for i, (prop, val) in enumerate(data.items()):
            with cols[i]:
                is_ok = val[1]
                st.metric(prop, f"{val[0]:.2f}", f"Limit: {val[2]}", delta_color="normal" if is_ok else "inverse")
                if is_ok: passes += 1
        
        if passes >= 3:
            st.success("Verdict: Molecule is likely orally bioactive (Drug-Like).")
        else:
            st.error("Verdict: Molecule fails drug-likeness criteria.")

elif menu == "Molecular Docking Demo":
    st.markdown('<div class="main-header">Molecular Docking Simulator</div>', unsafe_allow_html=True)
    st.write("Visualizing the interaction between Ligand and Protein Binding Pocket.")
    
    # 3D Visualization using st-py3dmol
    xyzview = st_py3dmol.make_viewer(height=500, width=800)
    # Using a PDB ID for the visual demo (HIV Protease)
    xyzview.addModel("pdb:1HSG")
    xyzview.setStyle({'cartoon': {'color': 'spectrum'}})
    xyzview.addStyle({'resn': 'MK1'}, {'stick': {'colorscheme': 'magentaCarbon'}}) # Highlight ligand
    xyzview.zoomTo()
    st_py3dmol.show(xyzview)
    st.info("The Magenta structure represents the ligand docked in the enzymatic pocket.")

elif menu == "ADMET Dashboard":
    st.markdown('<div class="main-header">ADMET Property Dashboard</div>', unsafe_allow_html=True)
    smiles = st.text_input("Ligand SMILES", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C") # Caffeine
    
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        tpsa = Descriptors.TPSA(mol)
        rotb = Lipinski.NumRotatableBonds(mol)
        
        # Radar Chart
        fig = go.Figure(data=go.Scatterpolar(
          r=[tpsa, rotb, 40, 70, 50],
          theta=['TPSA','Rotatable Bonds','Solubility','Absorption','Metabolism'],
          fill='toself'
        ))
        fig.update_layout(polar=dict(radialaxis=dict(visible=True, range=[0, 150])), showlegend=False)
        st.plotly_chart(fig)

elif menu == "Case Study Hub":
    st.markdown('<div class="main-header">AI Case Studies</div>', unsafe_allow_html=True)
    st.markdown("""
    - **AlphaFold (DeepMind):** Solved the 50-year-old 'protein folding problem', enabling rapid target identification.
    - **INS018_055 (Insilico Medicine):** The first AI-designed drug to reach Phase II clinical trials for IPF.
    - **Exscientia:** Developed the first AI-designed molecule to enter human clinical trials for OCD.
    """)
