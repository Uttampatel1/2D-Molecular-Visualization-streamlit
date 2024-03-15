from flask import Flask, render_template, request , jsonify
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from io import BytesIO
import base64

app = Flask(__name__)

def generate_3d_molfile_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        # Generate 3D coordinates
        AllChem.EmbedMolecule(mol, randomSeed=42)  # You can change the randomSeed value

        # Generate Molfile
        mol_block = Chem.MolToMolBlock(mol)
        return mol_block
    else:
        return None

def generate_molecule_image(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        img = Draw.MolToImage(mol, size=(500, 500))
        img_buffer = BytesIO()
        img.save(img_buffer, format="PNG")
        img_str = base64.b64encode(img_buffer.getvalue()).decode("utf-8")
        return img_str
    else:
        return None

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        smiles_input = request.form['smiles']

        # Validate SMILES input
        if not Chem.MolFromSmiles(smiles_input):
            return render_template('index.html', error='Invalid SMILES input')

        # Generate 2D structure
        mol_file = generate_3d_molfile_from_smiles(smiles_input)
        img_str = generate_molecule_image(smiles_input)

        if img_str is not None:
            # Get molecule information
            mol = Chem.MolFromSmiles(smiles_input)
            mol_info = {
                'Name': '',
                'Formula': rdMolDescriptors.CalcMolFormula(mol),
                'Exact Mass': Descriptors.ExactMolWt(mol),
                'Molecular Weight': Descriptors.MolWt(mol),
                'Canonical SMILES': Chem.MolToSmiles(mol, isomericSmiles=False),
                'Isomeric SMILES': Chem.MolToSmiles(mol, isomericSmiles=True)
            }

            return render_template('index.html', smiles=smiles_input, mol_img=img_str, mol_info=mol_info , mol_file=mol_file)
        else:
            return render_template('index.html', error='Failed to generate 2D structure for the given SMILES')

    return render_template('index.html')

if __name__ == '__main__':

    app.run(debug=False)
