import pandas as pd
import re
import numpy as np

# Atomic weights
atom_weight = {
    "C": {"mono": 12.00000000000, "avg": 12.011137185},
    "H": {"mono": 1.00782503223, "avg": 1.00794},
    "O": {"mono": 15.99491461956, "avg": 15.9994},
    "N": {"mono": 14.00307400443, "avg": 14.0067},
    "S": {"mono": 31.9720711744, "avg": 32.065},
    "P": {"mono": 30.97376199842, "avg": 30.973761},
    "Se": {"mono": 79.9165218, "avg": 78.971}
}

# Data for amino acids
amino_acids_info = {
    "A": {"three_letter": "Ala", "full_name": "Alanine", "formula": "C3H7NO2"},
    "C": {"three_letter": "Cys", "full_name": "Cysteine", "formula": "C3H7NO2S"},
    "D": {"three_letter": "Asp", "full_name": "Aspartic acid", "formula": "C4H7NO4"},
    "E": {"three_letter": "Glu", "full_name": "Glutamic acid", "formula": "C5H9NO4"},
    "F": {"three_letter": "Phe", "full_name": "Phenylalanine", "formula": "C9H11NO2"},
    "G": {"three_letter": "Gly", "full_name": "Glycine", "formula": "C2H5NO2"},
    "H": {"three_letter": "His", "full_name": "Histidine", "formula": "C6H9N3O2"},
    "I": {"three_letter": "Ile", "full_name": "Isoleucine", "formula": "C6H13NO2"},
    "K": {"three_letter": "Lys", "full_name": "Lysine", "formula": "C6H14N2O2"},
    "L": {"three_letter": "Leu", "full_name": "Leucine", "formula": "C6H13NO2"},
    "M": {"three_letter": "Met", "full_name": "Methionine", "formula": "C5H11NO2S"},
    "N": {"three_letter": "Asn", "full_name": "Asparagine", "formula": "C4H8N2O3"},
    "O": {"three_letter": "Pyl", "full_name": "Pyrrolysine", "formula": "C12H21N3O3"},
    "P": {"three_letter": "Pro", "full_name": "Proline", "formula": "C5H9NO2"},
    "Q": {"three_letter": "Gln", "full_name": "Glutamine", "formula": "C5H10N2O3"},
    "R": {"three_letter": "Arg", "full_name": "Arginine", "formula": "C6H14N4O2"},
    "S": {"three_letter": "Ser", "full_name": "Serine", "formula": "C3H7NO3"},
    "T": {"three_letter": "Thr", "full_name": "Threonine", "formula": "C4H9NO3"},
    "U": {"three_letter": "Sec", "full_name": "Selenocysteine", "formula": "C3H7NO2Se"},
    "V": {"three_letter": "Val", "full_name": "Valine", "formula": "C5H11NO2"},
    "W": {"three_letter": "Trp", "full_name": "Tryptophan", "formula": "C11H12N2O2"},
    "Y": {"three_letter": "Tyr", "full_name": "Tyrosine", "formula": "C9H11NO3"}
}

def string_to_formula(formula_string):
    parsed = re.findall('([A-Z][a-z]*)(\d*)', formula_string)
    formula = {}
    for elem, count in parsed:
        count = int(count) if count else 1
        formula[elem] = count
    return formula

def peptide_to_formula(sequence, C_term=None, N_term=None):
    total_formula = {} 
    for aa in sequence:
        try:
            formula_string = amino_acids_info[aa]["formula"]
            formula = string_to_formula(formula_string)
            for elem, count in formula.items():
                total_formula[elem] = total_formula.get(elem, 0) + count
        except KeyError:
            print(f"Invalid amino acid: {aa}")

    # Consider the loss of a water molecule (H2O) for each peptide bond formation, taking into account the modifications for C-terminal and N-terminal
    total_formula['H'] -= len(sequence) * 2
    total_formula['O'] -= len(sequence)
    # Consider the modifications for N-terminal and C-terminal
    if C_term == None: # No modification
        total_formula['H'] += 1
        total_formula['O'] += 1
    elif C_term == "NH2": # Amidation
        total_formula['H'] += 2
        total_formula['N'] += 1
    elif C_term == "OMe": # Methylation
        total_formula['H'] += 3
        total_formula['O'] += 1
        total_formula['C'] += 1
    else:
        raise ValueError("C_term type not supported now!")
    
    if N_term == None:
        total_formula['H'] += 1   
    else:
        raise ValueError("N_term type not supported now!")

    return total_formula

def formula_to_MW(formula, mass_type="mono"):
    weight = 0.0
    for elem, count in formula.items():
        elem, count = str(elem), int(count)
        weight += atom_weight[elem][mass_type] * count
    return weight

def calc_mw(peptide, C_term=None, N_term=None, mass_type="mono"):
    return formula_to_MW(peptide_to_formula(peptide, C_term=C_term, N_term=N_term), mass_type = mass_type)

def calc_mz(peptide,C_term=None, N_term=None, charge=1, mass_type="mono"):
    weight = calc_mw(peptide,C_term=C_term, N_term=N_term, mass_type = mass_type)
    return (weight + charge*atom_weight["H"][mass_type]) / charge

def mz_to_mw(mz, charge, mass_type="mono"):
    return mz * charge - charge * atom_weight["H"][mass_type]

if __name__ == "__main__":
    type = 'avg'
    peptide = "YCRTYWRYGRLRRRCYRRR"
    print(f"Peptide: {peptide}")
    print(f"Formula: {peptide_to_formula(peptide)}")
