import pandas as pd

def get_consensus(file):
    aa_info = pd.read_csv("./info/aa_info.csv").set_index("AA")

    with open(file,"r") as f:
        seq = []
        for line in f:
            if line.split()[0].isdigit:
                seq.append(tuple(line.split()))

    letters = [aa[1] for aa in seq[2:]]
    letters_three = [aa_info.index[aa_info["OneLetter"] == aa][0] for aa in letters]
    entropy = [abs(float(aa[2])) for aa in seq[2:]]
    consensus = pd.DataFrame({"letters": letters, "letters_three":letters_three, "entropy": entropy})

    return consensus