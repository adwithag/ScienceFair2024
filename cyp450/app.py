from shiny import App, render, ui
import pandas as pd
import matplotlib.pyplot as plt
import io
import numpy as np
import scipy
from rdkit.Chem import PandasTools
from rdkit import DataStructs
from rdkit.Chem import AllChem, MolFromSmiles, Draw
import xgboost as xgb
import pickle
from PIL import Image
from pathlib import Path
import tempfile

# list of cytochrome p-450 enzymes
cyps = ['CYP2C19', 'CYP2D6', 'CYP3A4', 'CYP1A2', 'CYP2C9', 'CYP23A4']

ml_models = {}
for cyp in cyps:
    with open("./CYP450/ML/model_"+cyp, 'rb') as file:
        ml_models[cyp] = pickle.load(file)

def style_negative(x):
    """Style negative values in red."""
    color = 'red' if x == "Yes" else 'green'
    return f'color: {color}'


# Sample function to generate an image
def customPredict(s, ml_models, cyps) -> pd.DataFrame:
    """ Predicts if a given molcule inhibits CYP or not.

     Args:
         s: SMILES representation of a molecule.
         ml_models: ML Models for each CYP Enzyme
         cyps : List of CYP450 enzymes
     Returns:
         A Dataframe the contains CYP Enzymes and Predictions.
    """
    cyps = ['CYP2C19', 'CYP2D6', 'CYP3A4', 'CYP1A2', 'CYP2C9', 'CYP23A4']
    mol = MolFromSmiles(s)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    arr = np.zeros((1,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    preds = ['Yes' if ml_models[i].predict(np.stack([arr]))[0] == 1 else 'No' for i in cyps]
    # pred = xgb_model.predict(np.stack([arr]))[0]

    return pd.DataFrame({'CYP450' : cyps, 'Inhibition' : preds})


# def mol_to_image(s):
#     """ Transforms SMILES Representation to an Image.

#      Args:
#          s: SMILES representation of a molecule.
         

#      Returns:
#          Image of a given molecule.
#     """
#     mol = MolFromSmiles(s)
#     img = Draw.MolToImage(mol)
#     return img

def mol_to_image(s):
    """Transforms SMILES representation to an image and saves it to a temporary file."""
    mol = MolFromSmiles(s)
    if mol is None:
        return None  # Handle invalid SMILES
    
    # Create a temporary file and get the file path
    with tempfile.NamedTemporaryFile(delete=False, suffix=".png") as temp_file:
        temp_file_path = temp_file.name  # Get the file path

    # Generate and save the image
    img = Draw.MolToImage(mol)
    img.save(temp_file_path, format="PNG")  # Save image

    return temp_file_path 

# UI definition
# app_ui = ui.page_fluid(
#     ui.h2('Smart Screening for Cytochrome P450 Inhibitors: AI Meets Drug Discovery'),
#     ui.h4('Adwaitha Gadachanda'),
#     ui.sidebar(
#         ui.input_text("user_text", "Enter a string:"),
#         ui.input_text(id="text_input", label="Enter SMILES Rep:", value="NC(=O)C1=CN([C@@H]2O[C@@H](COP(=O)([O-])OP(=O)([O-])OC[C@@H]3O[C@H](n4cnc5c(N)ncnc54)[C@@H](OP(=O)([O-])[O-])[C@@H]3O)[C@H](O)[C@@H]2O)C=CC1.[Na+].[Na+].[Na+].[Na+]")
#     )
    
#     # ui.layout_sidebar(
#     #     ui.sidebar(
#     #         ui.output_table("predictions")
#     #     ),
#     #     ui.output_image("image_display")
#     # )
# )

# Server logic
def server(input, output, session):
    
    @output
    # @render.ui
    # def predictions():
    #     pred_df = customPredict(input.smiles(), ml_models, cyps)
    #     pred_df = pred_df.reset_index(drop=True).style.map(style_negative, subset=['Inhibition'])
    #     return ui.HTML(pred_df.to_html(escape=False))  # Convert styled table to HTML

    # @render.table(index=False)
    # def predictions():
    #     pred_df = customPredict(input.smiles(), ml_models, cyps)
    #     pred_df = pred_df.style.map(style_negative, subset=['Inhibition'])
    #     return pred_df
    @render.data_frame
    def predictions():
        pred_df = customPredict(input.smiles(), ml_models, cyps)
        # pred_df = pred_df.style.map(style_negative, subset=['Inhibition'])
        return render.DataTable(pred_df,
            styles=[  
                # Center the text of each cell (using Bootstrap utility class) 
                # {   
                #     "class": "text-center",  
                # },  
                # Bold the first column 
                {  
                    "cols": [0],  
                    "style": {"font-weight": "bold"},  
                },  
                # Highlight the penguin colors 
                {
                  "rows": pred_df[pred_df["Inhibition"] == "Yes"].index.tolist(),
                  "cols": ["Inhibition"],
                  "style": {"backgroundColor": "white", "color": "red"}
              },
              {
                  "rows": pred_df[pred_df["Inhibition"] == "No"].index.tolist(),
                  "cols": ["Inhibition"],
                  "style": {"backgroundColor": "white", "color": "green"}
              },
              
            ]
            )

    @output
    @render.image
    def image_display():
        img_path = mol_to_image(input.smiles())
        if img_path is None:
            return None  # No image if input is invalid
        return {"src": img_path, "width": "100%",
        "style": "border: 1px solid black"}
    # def image_display():
    #     img = mol_to_image(input.smiles())  # Get RDKit Image

    #     temp_file = tempfile.NamedTemporaryFile(delete=False, suffix=".png")
    #     img.save(temp_file.name, format="png")
    #     # buf = io.BytesIO()
    #     # img.save(buf, format="PNG")  # Save RDKit image to buffer
    #     # buf.seek(0)
    #     # img.close()
    #     # data = buf.getvalue()
    #     # img_buf = Image.open(buf)
        
    #     # # Ensure the image is in RGB format
    #     # if img_buf.mode != "RGB":
    #     #     img_buf = img_buf.convert("RGB")
        
    #     # # Save the image to a temporary file
    #     # temp_path = Path("temp_image.png")
    #     # img_buf.save(temp_path)
    #     # f"data:image/png;base64,{data.decode('latin1')}"
    #     # return {
    #     #     "src": str(temp_path),
    #     #     "alt": "Generated plot",
    #     #     "style": "border: 1px solid black"
    #     # }
    #     return {
    #         "src": temp_file, "width": "100%",
    #         "alt": "Generated plot",
    #         "style": "border: 1px solid black"
    #     }
    #     # return {"src": buf, "width": "100%"}

app_ui = ui.page_fluid(
    ui.br(),
    ui.h1('Smart Screening for Cytochrome P450 Inhibitors: AI Meets Drug Discovery'),
    ui.h4('Adwitha Gadhachanda'),
    ui.br(),
    ui.input_text(id="smiles", label="Enter SMILES Rep:", value="NC(=O)C1=CN([C@@H]2O[C@@H](COP(=O)([O-])OP(=O)([O-])OC[C@@H]3O[C@H](n4cnc5c(N)ncnc54)[C@@H](OP(=O)([O-])[O-])[C@@H]3O)[C@H](O)[C@@H]2O)C=CC1.[Na+].[Na+].[Na+].[Na+]"),
    ui.br(),
    ui.layout_sidebar(
        ui.sidebar(
            ui.h4('Molecule'),
            ui.output_image("image_display"),
            width = 400
        ),
        ui.h4('Inhibits CYP450 or Not'),
        # ui.output_ui("predictions", style='width: 300px; border: 1px solid black;')
        ui.output_data_frame("predictions")
        
    )
)

# def server(input, output, session):
#     @render.text
#     def greeting():
#         return f"Hello, {input.name()}! You are {input.age()} years old and live in {input.city()}."

app = App(app_ui, server)
