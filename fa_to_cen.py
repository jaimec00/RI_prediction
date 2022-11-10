from pyrosetta import *
from pyrosetta.toolbox.cleaning import cleanATOM
pyrosetta.init()
import os

pdb_path = "/home/hejaca00/ri_structure/pdbs"
clean_path = "/home/hejaca00/ri_structure/clean_pdbs"


for pdb in os.listdir(pdb_path):

	cleanATOM(os.path.join(pdb_path, pdb), f"{clean_path}/{pdb}.clean.pdb")
	pose = pyrosetta.pose_from_pdb(f"{clean_path}/{pdb}.clean.pdb")

	xml_obj = pyrosetta.rosetta.protocols.rosetta_scripts.RosettaScriptsParser()
	protocol = xml_obj.generate_mover('fa_to_cen.xml')

	protocol.apply(pose)
	pose.dump_pdb(f"/home/hejaca00/ri_structure/centroid_pdbs/{pdb[:-4]}_cen.pdb")
