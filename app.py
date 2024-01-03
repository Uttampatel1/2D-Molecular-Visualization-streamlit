import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors

def generate_2d_structure_image(smiles_input, image_id):
    mol = Chem.MolFromSmiles(smiles_input)
    if mol is not None:
        return image_id, Draw.MolToImage(mol, size=(400, 400))  # Adjusted size for side-by-side display
    else:
        return image_id, None

def get_molecular_formula_and_weight(smiles_input):
    mol = Chem.MolFromSmiles(smiles_input)
    if mol is not None:
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
        weight = Descriptors.MolWt(mol)
        return formula, weight
    else:
        return None, None

def main():
    st.set_page_config(page_title="2D Molecular Visualization", page_icon="🧪", layout="wide")
    st.title("2D Molecular Visualization 🌐")

    # Larger input box
    smiles_input = st.text_area("Enter SMILES strings 📝:", height=100)
    image_id = 0  # Counter for image IDs

    if smiles_input:
        smiles_inputs = smiles_input.split('\n')
        smiles_inputs = [smiles.strip() for smiles in smiles_inputs]
        smiles_inputs = [smiles.replace('"','') for smiles in smiles_inputs]
        smiles_inputs = [smiles.replace("\n","") for smiles in smiles_inputs ]
        smiles_inputs = [smiles.upper() for smiles in smiles_inputs if smiles]
        
        smiles_list = smiles_inputs
        
        # Display images and information in columns
        col1, col2, col3, col4, col5 = st.columns(5)

        # Button to generate images
    if st.button("Generate Images"):
        if smiles_input:
            for i, smiles in enumerate(smiles_list):
                image_id, img_2d = generate_2d_structure_image(smiles, image_id)
                with col1 if i % 5 == 0 else col2 if i % 5 == 1 else col3 if i % 5 == 2 else col4 if i % 5 == 3 else col5:
                    st.subheader(f"ID : {i+1}")
                    if img_2d is not None:
                        st.image(img_2d, use_column_width=True)
                        formula, weight = get_molecular_formula_and_weight(smiles)
                        if formula and weight:
                            st.write(f"SMILES: {smiles}")
                            st.write(f"Molecular Formula: {formula} 🧪")
                            st.write(f"Molecular Weight: {weight:.2f} g/mol ⚖️")
                        else:
                            st.warning(f"Unable to calculate formula and weight for {smiles}")
                    else:
                        st.warning(f"Invalid SMILES string: {smiles} ⚠️")
        else:
            st.warning(f"Input SMILES string ⚠️")

if __name__ == "__main__":
    main()
