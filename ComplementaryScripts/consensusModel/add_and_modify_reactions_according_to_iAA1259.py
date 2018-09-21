# -*- coding: utf-8 -*-
"""
This file adds reactions and apply selected changes from iAA1259 based on the mapping of reactions and metabolites. 
Author: Snorre Sulheim
Date: 21.09.2018

# Description

# Files required

"""

import cobra
import pandas as pd
from collections import defaultdict
import re
import logging
import csv

# from compare_iMK1208_content import fix_iAA1259

S4_FN = "../../ComplementaryData/curation/iAA1259_suppl_S4.csv" # New reactions
S5_FN = "../../ComplementaryData/curation/iAA1259_suppl_S5.csv" # New metabolites
scoGEM_FN = "../../ModelFiles/xml/scoGEM.xml"
iAA1259_FN = r"..\\..\\ComplementaryData\\models\\iAA1259.xml"


metabolite_mappings = {
    "3cvobz_c": "CPD__45__16467_c",
    "6adxfut": "CPD__45__13324_c"
}

reaction_mappings = {
    "CHDHR": "RXN__45__12345",
    "ADXFUTSNT": "RXN__45__15264",
    "ADXFUTDA": "RXN__45__14910",
    "CYOO": "CYTBD2", # Check gene annotation
    "AMMQT9r": "AMMQT9", #Bounds changed and genes

}

def add_reactions(iAA1259_model, reaction_mapping_fn, metabolite_mapping_fn, add_new_metabolites = True):
    pass

def modify_reactions():
    pass


def get_suppl_sheets():
    S4_df = pd.read_csv(S4_FN, sep = ";")

def map_S4():
    S4_df = pd.read_csv(S4_FN, sep = ";")
    print(S4_df)
    scoGEM = cobra.io.read_sbml_model(scoGEM_FN)
    iAA1259 = cobra.io.read_sbml_model(iAA1259_FN)
    for r in iAA1259.reactions:
        r.id = r.id.replace("_LPAREN_e_RPAREN_", "(e)").replace("_LPAREN_c_RPAREN_", "(c)").replace("_DASH_", "-")

    for index, row in S4_df.iterrows():
        r_id = row["Reaction ID"]
        iAA1259_r = iAA1259.reactions.get_by_id(r_id)
        print(iAA1259_r.id, iAA1259_r.annotation)

        try:
            scoGEM.reactions.get_by_id(r_id)
        except KeyError:
            pass
        else:
            print(r_id, "  Exists!!")
    


if __name__ == '__main__':
    map_S4()
    # with open(S4_FN, "r") as f:
    #     reader = csv.reader(f, delimiter = ";")
    #     for row in reader:
    #         print(len(row), row)
