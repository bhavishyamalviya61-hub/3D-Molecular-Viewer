
import streamlit as st
import py3Dmol
import streamlit.components.v1 as components

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

st.set_page_config(
    page_title="3D Molecular Viewer",
    page_icon="🧪",
    layout="wide"
)

st.title("🧪 3D Molecular Structure Viewer")
st.write(
    "Explore molecular structures and their chemical properties."
)

# Molecule name and SMILES database
MOLECULES = {
    "Water (H2O)": "O",
    "Methane (CH4)": "C",
    "Ammonia (NH3)": "N",
    "Carbon dioxide (CO2)": "O=C=O",
    "Ethanol (C2H6O)": "CCO",
    "Benzene (C6H6)": "c1ccccc1",
    "Acetic acid (C2H4O2)": "CC(=O)O",
    "Glucose (C6H12O6)": "C(C1C(C(C(C(O1)O)O)O)O)O"
}

selected = st.selectbox(
    "Choose a molecule",
    list(MOLECULES.keys())
)

if st.button("Generate 3D Molecule", type="primary"):
    with st.spinner("Generating molecular structure..."):
        try:
            smiles = MOLECULES[selected]

            # Create molecule from SMILES
            mol = Chem.MolFromSmiles(smiles)

            if mol is None:
                raise ValueError("Could not read molecule.")

            # Add hydrogen atoms
            mol = Chem.AddHs(mol)

            # Generate 3D coordinates
            result = AllChem.EmbedMolecule(
                mol,
                randomSeed=42
            )

            if result != 0:
                raise ValueError(
                    "Could not generate 3D coordinates."
                )

            # Optimize geometry using force-field methods
            if AllChem.MMFFHasAllMoleculeParams(mol):
                AllChem.MMFFOptimizeMolecule(mol)
            else:
                AllChem.UFFOptimizeMolecule(mol)

            # Convert molecule to 3D SDF format
            mol_block = Chem.MolToMolBlock(mol)

            # Calculate molecular properties
            formula = rdMolDescriptors.CalcMolFormula(mol)
            weight = Descriptors.MolWt(mol)
            atom_count = mol.GetNumAtoms()
            bond_count = mol.GetNumBonds()

            st.subheader(selected)

            col1, col2, col3, col4 = st.columns(4)

            col1.metric("Formula", formula)
            col2.metric("Molecular Weight", f"{weight:.3f} g/mol")
            col3.metric("Atoms", atom_count)
            col4.metric("Bonds", bond_count)

            # Create interactive 3D viewer
            viewer = py3Dmol.view(
                width=900,
                height=520
            )

            viewer.addModel(mol_block, "sdf")
            viewer.setStyle(
                {"stick": {"radius": 0.15},
                 "sphere": {"scale": 0.25}}
            )
            viewer.setBackgroundColor("white")
            viewer.zoomTo()

            components.html(
                viewer._make_html(),
                height=540,
                scrolling=False
            )

            st.success("3D molecular structure generated!")

            with st.expander("View chemical information"):
                st.write("**SMILES:**", smiles)
                st.write("**Formula:**", formula)
                st.write("**Molecular weight:**", f"{weight:.3f} g/mol")
                st.write("**Number of atoms:**", atom_count)
                st.write("**Number of bonds:**", bond_count)

            st.caption(
                "The displayed structure is a computed 3D model. "
                "Geometry is optimized with a molecular mechanics "
                "force field, not a full quantum-chemical calculation."
            )

        except Exception as error:
            st.error(f"Error generating molecule: {error}")

st.divider()
st.caption(
    "Computational Chemistry Mini Project | Python + RDKit + py3Dmol"
)
