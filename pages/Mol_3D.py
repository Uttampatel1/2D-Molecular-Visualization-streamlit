import streamlit as st
import streamlit.components.v1 as components

def main():
    st.set_page_config(page_title="3D Molecular Visualization", page_icon="🧪")
    st.title("3D Molecular Visualization 🌐")

    # Larger input box
    smiles_input = st.text_area("Enter SMILES strings 📝:", height=100)
    image_id = 0  # Counter for image IDs

    if smiles_input:
        smiles_inputs = smiles_input.split('\n')
        smiles_inputs = [smiles.strip() for smiles in smiles_inputs]
        smiles_inputs = [smiles.replace('"','') for smiles in smiles_inputs]
        smiles_inputs = [smiles.replace("\n","") for smiles in smiles_inputs ]
        smiles_inputs = [smiles for smiles in smiles_inputs if smiles]
        
        smiles_list = smiles_inputs
    
    if st.button("Generate 3D Images"):
        if smiles_input:
            for i, smiles in enumerate(smiles_list):
                # Dynamic update of the embedded MolView iframe based on user input
                st.subheader(f"ID : {i+1}")
                try:
                    molview_url = f"https://embed.molview.org/v1/?mode=balls&smiles={smiles}"
                    components.iframe(molview_url, width=300, height=300, scrolling=True)
                    st.write(f"SMILES: {smiles}")
                except:
                    st.warning(f"Invalid SMILES string: {smiles} ⚠️")
                

if __name__ == "__main__":
    main()
