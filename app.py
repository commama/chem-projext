import streamlit as st
import pandas as pd
import plotly.graph_objects as go
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
import st_py3dmol

# --- UI Configuration ---
st.set_page_config(page_title="In Silico Drug Design Lab", layout="wide")

# Professional styling
st.markdown("""
    <style>
    .main {background-color: #f5f7f9;}
    .stMetric {background-color: #ffffff; padding: 15px; border-radius: 10px; border: 1px solid #e0e0e0;}
    h1 {color: #003366;}
    </style>
    """, unsafe_allow_html=True)

# --- Sidebar ---
st.sidebar.title("🔬 Unit II: Lab Modules")
app_mode = st.sidebar.selectbox("Select Module", 
    ["Lipinski's Rule Checker", "Molecular Docking Visualizer", "ADMET Properties", "Case Study Hub"])

# --- Module 1: Lipinski ---
if app_mode == "Lipinski's Rule Checker":
    st.header("Lipinski's Rule of Five Checker")
    st.write("Evaluate the drug-likeness of a molecule based on its physicochemical properties.")
    
    smiles = st.text_input("Enter SMILES string", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C") # Caffeine default
    
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        
        cols = st.columns(4)
        cols[0].metric("Mol Weight", f"{mw:.1f}", "Limit < 500")
        cols[1].metric("LogP", f"{logp:.2f}", "Limit < 5")
        cols[2].metric("H-Donors", hbd, "Limit < 5")
        cols[3].metric("H-Acceptors", hba, "Limit < 10")
        
        # Logic check
        if mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10:
            st.success("✅ Verdict: Molecule is likely to be orally active.")
        else:
            st.warning("⚠️ Verdict: Molecule may have poor oral bioavailability.")
    else:
        st.error("Invalid SMILES string. Please enter a valid chemical structure.")

# --- Module 2: Docking ---
elif app_mode == "Molecular Docking Visualizer":
    st.header("Molecular Docking Simulator (3D)")
    st.info("Demonstrating Ligand-Protein interaction (Target: HIV-1 Protease)")
    
    # Create the 3D viewer
    view = st_py3dmol.make_viewer(height=500, width=800)
    # Using PDB 1HSG (HIV Protease with Indinavir)
    view.addModel("pdb:1hsg")
    view.setStyle({'cartoon': {'color': 'spectrum'}})
    view.addStyle({'resn': 'MK1'}, {'stick': {'colorscheme': 'magentaCarbon'}})
    view.zoomTo()
    st_py3dmol.show(view)
    st.caption("Visualization: Protein (Cartoon) and Ligand (Sticks) at the active site.")

# --- Module 3: ADMET ---
elif app_mode == "ADMET Properties":
    st.header("ADMET Property Dashboard")
    
    # Mock data for visualization based on a standard pharmacokinetic profile
    categories = ['Absorption', 'Distribution', 'Metabolism', 'Excretion', 'Toxicity']
    values = [85, 60, 45, 70, 20] # Representative values
    
    fig = go.Figure(data=go.Scatterpolar(
        r=values + [values[0]],
        theta=categories + [categories[0]],
        fill='toself',
        line_color='#003366'
    ))
    
    fig.update_layout(polar=dict(radialaxis=dict(visible=True, range=[0, 100])))
    st.plotly_chart(fig)
    st.write("This radar chart visualizes the balanced ADMET profile required for a successful drug candidate.")

# --- Module 4: Case Studies ---
elif app_mode == "Case Study Hub":
    st.header("AI in Drug Discovery Case Study Hub")
    
    with st.expander("Case 1: AlphaFold (DeepMind)"):
        st.write("Revolutionized target identification by providing high-accuracy 3D protein structure predictions.")
    
    with st.expander("Case 2: INS018_055 (Insilico Medicine)"):
        st.write("The first AI-discovered drug for Idiopathic Pulmonary Fibrosis to reach Phase II clinical trials.")
