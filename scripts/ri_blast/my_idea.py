# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
import pandas as pd
from Bio import pairwise2
import numpy as np
from consensus import get_consensus
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# get protein info
aa_info = pd.read_csv("c:\\users\\hejac\\OneDrive - UNC-Wilmington\\wang_lab\\crystallin_proj\\ri_structure\\my_idea\\info\\aa_info.csv").set_index("AA")
protein_info = pd.read_csv("c:\\users\\hejac\\OneDrive - UNC-Wilmington\\wang_lab\\crystallin_proj\\ri_structure\\my_idea\\info\\protein_info.csv").set_index("Protein")
# ------------------------------------------------------------------------------
# get Consensus
protein = "yS"
consensus = get_consensus(f"{protein}_consensus.txt")
consensus_str = "".join(consensus.letters)

# get MSA info
with open(f"./alignments/{protein}_alignment.fa", "r") as f:
    seqs = []
    seq_inside = ""
    for line in f:
        if line[0] != ">":
            seq_inside = "".join([seq_inside, line[:-1]])
        elif seq_inside != "":
            seqs.append(seq_inside)
            seq_inside = ""
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
sequence_one = protein_info.loc[protein, "Sequence"]
sequence = [aa_info.index[aa_info["OneLetter"]==aa][0] for aa in sequence_one]
# ------------------------------------------------------------------------------
messy_alignment = pairwise2.align.globalxd(sequenceA=consensus_str, sequenceB=sequence_one,
                                                openA=-100000000,
                                                extendA=-100000000,
                                                openB=-6,
                                                extendB=-0.5)


alignment = pairwise2.format_alignment(*messy_alignment[0])

sequence_aligned = alignment[len(consensus_str)*2+2:len(consensus_str)*3+2]
sequence_aligned_three = [aa_info.index[aa_info["OneLetter"]==aa][0] for aa in sequence_aligned]
# ------------------------------------------------------------------------------
# get index of aligned sequence
sequence_idx = []
for idx, position in enumerate(sequence_aligned):
    if position != "-":
        sequence_idx.append(idx)



aa_all = []
entropies = []
for seq in seqs:
    for idx, aa in enumerate(seq):
        if idx in sequence_idx:
            aa_all.append(aa_info.index[aa_info["OneLetter"]==aa][0])
            entropies.append((max(consensus.entropy)-consensus.loc[idx, "entropy"])/max(consensus.entropy))

sequence_dndc_orig = sum(np.array([aa_info.loc[aa, "dn/dc"] for aa in sequence]) * np.array([aa_info.loc[aa, "MW"] for aa in sequence])) / sum([aa_info.loc[aa, "MW"] for aa in sequence])
sequence_dndc_new = sum(np.array([aa_info.loc[aa, "dn/dc"] for aa in aa_all]) * np.array([aa_info.loc[aa, "MW"] for aa in aa_all]) ) / sum([aa_info.loc[aa, "MW"] for aa in aa_all])

correction = sum(np.array([aa_info.loc[aa, "dn/dc"] for aa in sequence_aligned_three]) * np.array([aa_info.loc[aa, "MW"] for aa in sequence_aligned_three]) * np.array(consensus.entropy)) / sum([aa_info.loc[aa, "MW"] for aa in aa_all])


print(sequence_dndc_orig, sequence_dndc_new, correction)
print(protein_info.loc[protein, "Experimental RI"])