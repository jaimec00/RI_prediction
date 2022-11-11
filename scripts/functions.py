# Function to extract PDB info
from collections import defaultdict
from protein_objects import Atom, Amino_Acid, Protein
from matplotlib import pyplot as plt
import numpy as np 
import os


def get_pdb_info(pdb_file, name):
    # This will have all the atoms in the protein and their info
    protein_info = defaultdict(list)
    if "cen" in pdb_file: 
        # Read PDB
        with open(pdb_file, "r") as pdb:
            # Extract info from "ATOM" + "CEN" lines and add to protein_info
            for line in pdb:
                if "ATOM" in line[:4] and "CEN" in line[12:17]:
                    pos = f"{line[21:23].strip()}_{line[23:31].strip()}"
                    aa = line[17:20].strip()
                    atom = line[12:17].strip()
                    x = float(line[31:39].strip())
                    y = float(line[38:47].strip())
                    z = float(line[47:56].strip())
                    protein_info[pos].append((aa, atom, x, y, z))
    else:
        # Read PDB
        with open(pdb_file, "r") as pdb:
            # Extract info from "ATOM" + "CEN" lines and add to protein_info
            for line in pdb:
                if "ATOM" in line[:4]:
                    pos = f"{line[21:23].strip()}_{line[23:31].strip()}"
                    aa = line[17:20].strip()
                    atom = line[12:17].strip()
                    x = float(line[31:39].strip())
                    y = float(line[38:47].strip())
                    z = float(line[47:56].strip())
                    protein_info[pos].append((aa, atom, x, y, z))
    
    # Create Atom, Amino_Acid, and Protein objects from info  
    amino_acids = list()
    for residue in protein_info.keys():    
        
        atoms = list()
        for atom in protein_info[residue]:
            atoms.append(Atom(atom[1], atom[2], atom[3], atom[4]))
        
        amino_acids.append(Amino_Acid(protein_info[residue][0][0], atoms, residue))
    
    protein = Protein(name, amino_acids)
    
    return protein


def get_dndc_predictions(pdb_list, path):
    dn_dc_list = np.empty(0)
    for pdb in pdb_list:
        name = pdb
        if "_cen" in name:
            name = name[:-8] + name[-4:] 
        protein = get_pdb_info(os.path.join(path, pdb), name)
        
        # Calculate the dn_dc and the correction factor
        dn_dc = protein.get_dn_dc()
        cf = protein.get_cf()
        
        dn_dc_list = np.append(dn_dc_list, dn_dc+cf)
        
    return dn_dc_list
    
def graph(values, names):
    
    colors = ["black", "red","blue","green", "orange"]
    true  = values["Experimental dn/dc"]

    for color, prediction in enumerate(values.keys()):
        
        slope, intercept = np.polyfit(true, values[prediction], 1)
        correlation = np.corrcoef(true, values[prediction])

        plt.scatter(true, values[prediction], c=colors[color], linewidths=0.25, alpha=0.75)
        plt.plot(true, slope*np.array(true) + intercept, c=colors[color], linewidth=0.75,label=f"{prediction} | R: {round(correlation[1][0], 2)}", alpha=0.75)
    
    plt.ylabel("Predicted dn/dc")
    plt.xlabel("Experimental dn/dc")
        
    # plt.plot([0,2], [0,2], c="black", label = "Experimental dn/dc")
    for name, point in enumerate(true):
        plt.text(x=point - 0.0005, y=point+ 0.0005, s=names[name], fontsize=7.5)

    plt.axis([0.197, 0.22, 0.19, 0.22])
    plt.legend(loc="upper left", fontsize="x-small")

    plt.show()
