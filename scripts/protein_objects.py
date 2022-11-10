import pandas as pd

# INFO 
path = r"C:\Users\hejac\OneDrive - UNC-Wilmington\wang_lab\crystallin_proj\ri_structure"
aa_info = pd.read_csv(path + r"\info\aa_info.csv").set_index("AA")

# Define Classes
class Atom:
    def __init__(self, name, x, y, z, radius=0): 
        self.name = name
        if "O" in name:
            self.weight = 15.999
        elif "C" in name:
            self.weight = 12.011
        elif "N" in name:
            self.weight = 14.007

        else:
            self.weight=0
        self.x = x
        self.y = y
        self.z = z
        self.radius = radius


    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name
# ----------------------------------------------------------------------------------------------------------------------
class Amino_Acid:
    def __init__(self, aa, atoms, pos) -> None:
        self.name = aa
        self.mw = aa_info.loc[aa, "MW"]
        self.specific_volume = aa_info.loc[aa, "Specific Volume"]
        self.molar_refractivity = aa_info.loc[aa, "Molar Refraction"]
        self.refractivity_per_gram = aa_info.loc[aa, "Refraction per gram"]
        self.ri = aa_info.loc[aa, "Refractive Index"]
        self.dn_dc = aa_info.loc[aa, "dn/dc"]
        self.atoms = atoms
        
        if next(atom.name for atom in atoms) == "CEN":
            self.centroid = next(atom for atom in atoms if atom.name == "CEN")
        else:
            if aa not in  ["GLY", "SER", "GLU", "ALA"]:
                backbone = ["N", "CA", "C", "O", "CB", "H"]
            else:
                backbone = []
            
            # Center of mass
            # x_avg = sum([atom.x*atom.weight for atom in atoms if atom.name not in backbone and "H" not in atom.name]) / sum([atom.weight for atom in atoms if atom not in backbone and "H" not in atom.name])# len([atom for atom in atoms if atom.name not in backbone and "H" not in atom.name])
            # y_avg = sum([atom.y*atom.weight for atom in atoms if atom.name not in backbone and "H" not in atom.name]) / sum([atom.weight for atom in atoms if atom not in backbone and "H" not in atom.name]) # len([atom for atom in atoms if atom.name not in backbone and "H" not in atom.name])
            # z_avg = sum([atom.z*atom.weight for atom in atoms if atom.name not in backbone and "H" not in atom.name]) / sum([atom.weight for atom in atoms if atom not in backbone and "H" not in atom.name]) # len([atom for atom in atoms if atom.name not in backbone and "H" not in atom.name])

            # Average position
            x_avg = sum([atom.x for atom in atoms if atom.name not in backbone and "H" not in atom.name]) / len([atom for atom in atoms if atom.name not in backbone and "H" not in atom.name])
            y_avg = sum([atom.y for atom in atoms if atom.name not in backbone and "H" not in atom.name]) / len([atom for atom in atoms if atom.name not in backbone and "H" not in atom.name])
            z_avg = sum([atom.z for atom in atoms if atom.name not in backbone and "H" not in atom.name]) / len([atom for atom in atoms if atom.name not in backbone and "H" not in atom.name])

            all_dists = [((atom.x - x_avg)**2 + (atom.y - y_avg)**2 + (atom.z - z_avg)**2)**(1/2) for atom in atoms if atom.name not in backbone and "H" not in atom.name]
            radius = sum(all_dists) / len(all_dists)
            self.centroid = Atom("CEN", x_avg, y_avg, z_avg, radius)

        self.pos = pos
        if self.name in ["HIS","PHE","TRY","TYR"]:
            self.cutoff = 7 # angstroms
        elif self.name == "ARG":
            self.cutoff = 6  # angstroms
    
    def __repr__(self):
        return self.name
# ----------------------------------------------------------------------------------------------------------------------
class Protein:
    def __init__(self, name, sequence) -> None:
        self.name = name
        self.sequence = sequence
        self.refractivity_per_gram = sum([aa.refractivity_per_gram * aa.mw for aa in sequence] ) / sum([aa.mw for aa in sequence])
        self.specific_vol = sum([aa.specific_volume * aa.mw for aa in sequence]) / sum([aa.mw for aa in sequence])
        self.mw = sum([aa.mw for aa in sequence])
        self.centroid = sequence[0].centroid

    def get_dn_dc(self):
        temp = 25
        wavelength = 589.3
        dn_dc_raw = sum([aa.mw * aa.dn_dc for aa in self.sequence]) / sum([aa.mw for aa in self.sequence])
        dn_dc = dn_dc_raw * (0.94 + (20000/(wavelength**2))) * (1 + (25 - temp)*(0.0005 / 30))

        return dn_dc

    def get_dn_dc_ii(self):
        ref_idx_water = 1.334
        np = ((2*self.refractivity_per_gram + self.specific_vol) / (self.specific_vol - self.refractivity_per_gram)) ** (1/2)
        dn_dc = (3/2)*self.specific_vol*ref_idx_water*((np**2-ref_idx_water**2) / (np**2 + 2*ref_idx_water**2) )

        return dn_dc

    def get_cf(self):
        polar_aa = [aa for aa in self.sequence if aa.name in ["HIS","PHE","TRY","TYR","ARG"]]
        cf = 0
        
        for idx, polar_i in enumerate(polar_aa):
            for polar_j in polar_aa[idx+1:]:
        
                dist = (((polar_i.centroid.x-polar_j.centroid.x)**2 + (polar_i.centroid.y-polar_j.centroid.y)**2 + (polar_i.centroid.z-polar_j.centroid.z)**2)**(1/2)) - (polar_i.centroid.radius+polar_j.centroid.radius)

                if polar_i.cutoff > dist or polar_i.cutoff > dist:
                    if polar_i.name == "ARG" or polar_j.name == "ARG":
                        delta = -0.08
                    else:
                        delta = 0.11
                    cf += ((polar_i.dn_dc) * (polar_j.dn_dc))  / dist**(3)

        return cf
