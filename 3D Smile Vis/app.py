from flask import Flask, render_template, request, jsonify
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors, rdMolDescriptors
from io import BytesIO
import base64

app = Flask(__name__, static_folder='./static', template_folder='./templates')


def generate_3d_molfile_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        AllChem.EmbedMolecule(mol, randomSeed=42)
        mol_block = Chem.MolToMolBlock(mol)
        return mol_block
    else:
        return None


def generate_molecule_image(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        img = Draw.MolToImage(mol, size=(300, 250))
        img_buffer = BytesIO()
        img.save(img_buffer, format="PNG")
        img_str = base64.b64encode(img_buffer.getvalue()).decode("utf-8")
        return img_str
    else:
        return None


@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        # smiles_inputs = request.form.getlist('smiles')
        smiles_strings = request.form['smiles']

        # Split the input string into a list of SMILES strings
        smiles_inputs = smiles_strings.split(',')

        # Trim leading and trailing whitespaces from each SMILES string
        smiles_inputs = [smiles.strip() for smiles in smiles_inputs]

        print(smiles_inputs)

        error_smiles = []
        # Validate SMILES input
        for smiles_input in smiles_inputs:
            if not Chem.MolFromSmiles(smiles_input):
                # return render_template('index.html', error='Invalid SMILES input')
                error_smiles.append(smiles_input)
                smiles_inputs.remove(smiles_input)

        mol_data = []
        for smiles_input in smiles_inputs:
            mol_file = generate_3d_molfile_from_smiles(smiles_input)
            img_str = generate_molecule_image(smiles_input)

            if img_str is not None:
                mol = Chem.MolFromSmiles(smiles_input)
                mol_info = {
                    'Name': '',
                    'Formula': rdMolDescriptors.CalcMolFormula(mol),
                    'Exact Mass': round(Descriptors.ExactMolWt(mol), 2),
                    'Molecular Weight': round(Descriptors.MolWt(mol), 2),
                    'Canonical SMILES': Chem.MolToSmiles(mol, isomericSmiles=False),
                    'Isomeric SMILES': Chem.MolToSmiles(mol, isomericSmiles=True)
                }
                mol_data.append({'smiles': smiles_input, 'mol_img': img_str, 'mol_info': mol_info, 'mol_file': mol_file})
            else:
                return render_template('index.html', error='Failed to generate 2D structure for the given SMILES')

        print("Error String are : ",error_smiles)
        return render_template('index.html', mol_data=mol_data)

    return render_template('index.html')


if __name__ == '__main__':
    app.run(debug=False)
