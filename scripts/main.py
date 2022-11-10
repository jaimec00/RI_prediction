# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
from functions import get_dndc_predictions, graph
import pandas as pd
import os
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
def main():

    rosetta_pdb_list = sorted([file for file in os.listdir(rosetta_path)])
    manual_pdb_list = sorted([file for file in os.listdir(pdb_path)])

    values = {"Experimental dn/dc": protein_info.sort_values(by="pdb file").loc[protein_info["pdb file"].isin(manual_pdb_list), "Experimental RI"].values, 
                "Paper Prediction": protein_info.sort_values(by="pdb file").loc[protein_info["pdb file"].isin(manual_pdb_list), "Predicted RI"].values,
                "Paper Prediction + CF": protein_info.sort_values(by="pdb file").loc[protein_info["pdb file"].isin(manual_pdb_list), "Predicted RI + cf"].values, 
                "Rosetta Centroid Prediction": get_dndc_predictions(rosetta_pdb_list, rosetta_path), 
                "Manual Centroid Prediction": get_dndc_predictions(manual_pdb_list, pdb_path)}

    print(values["Experimental dn/dc"])
    graph(values)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# INFO 
rosetta_path = r"C:\Users\hejac\OneDrive - UNC-Wilmington\wang_lab\crystallin_proj\ri_structure\structures\centroid_pdbs"
pdb_path = r"C:\Users\hejac\OneDrive - UNC-Wilmington\wang_lab\crystallin_proj\ri_structure\structures\pdbs"
protein_info = pd.read_csv(r"C:\Users\hejac\OneDrive - UNC-Wilmington\wang_lab\crystallin_proj\ri_structure\info\protein_info.csv").set_index("Protein").drop("J2").drop("Lysozyme")
# ------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------