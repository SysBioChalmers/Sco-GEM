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
from ast import literal_eval

# from compare_iMK1208_content import fix_iAA1259
from add_reactions_from_sco4 import check_metabolites

S4_FN = "../../ComplementaryData/curation/iAA1259_suppl_S4.csv" # New reactions
S5_FN = "../../ComplementaryData/curation/iAA1259_suppl_S5.csv" # New metabolites
scoGEM_FN = "../../ModelFiles/xml/scoGEM.xml"
iAA1259_PATH = "../../ComplementaryData/models/iAA1259.xml"


METABOLITE_MAPPINGS = {
    "3cvobz_c": "CPD__45__16467_c",
    "6adxfut_c": "6a6doxf_c",
    "dad__5_c": "dad_5_c",
    "celb_e": "cellb_e",
    #"6adxfut_c"
}

reaction_mappings = {
# Not used, but was used as notes when idenentifying reactions that should be added
    "CHDHR": "RXN__45__12345",
    "ADXFUTSNT": "RXN__45__15264",
    "ADXFUTDA": "RXN__45__14910",
    "CYOO": "CYTBD2", # Check gene annotation
    "AMMQT9r": "AMMQT9", #Bounds changed and genes
}

# Changed reactions
"""These is a manually curated dict of reactions that are modified in iAA1259.
There were a few reactions where the number of hydrogens were changed, but this resulted in unbalanced 
equations and were not included"""
CHANGED_REACTIONS = {
    "ABTDG":    "genes:SCO0737",
    "4HGLSD":   "genes:SCO5520 or SCO3835",
    "P5CD":     "genes:SCO5520 or SCO3835",
    "PHCD":     "genes:SCO5520 or SCO3835",
    "PUTA3":    "genes:SCO5520 or SCO3835",
    "UDPGD":    "genes:SCO3052 or SCO0382",
    "AOXSr2":   "genes:SCO1243",
    "BACCL":    "genes:SCO4927 or SCO4655",
    "CFL":      "remove:",
    "DHFUTALS": "remove:",
    "DHNANT4":  "genes:SCO4491 and SCO4556",
    "GLUTRS":   "genes:SCO5547 or SCO5498 or SCO5499",
    "GLBRAN2":  "genes:SCO5440 or SCO7332 or SCO0765 or SCO0554 or SCO7637",
    "FERO":     "genes:SCO3439 and SCO3440",
    "HACD1":    "genes:SCO6732 or SCO0984",
    "HACD2":    "genes:SCO6732 or SCO0984",
    "HACD3":    "genes:SCO6732 or SCO0984",
    "HACD4":    "genes:SCO6732 or SCO0984",
    "HACD5":    "genes:SCO6732 or SCO0984",
    "HACD6":    "genes:SCO6732 or SCO0984",
    "HACD7":    "genes:SCO6732 or SCO0984",
    "HACD8":    "genes:SCO6732 or SCO0984",
    "CYO2b":    "genes:SCO1934 and SCO2156 and (SCO7234 or SCO2155) and SCO2151 and SCO1930 and SCO2154",
    "CYTBD2":   "genes:SCO7234 and SCO7235 and SCO7236 and SCO1934 and SCO7120",
    "NADH17b":   "genes:(SCO4562 or SCO4599) and (SCO4563 or SCO4600) and SCO4564 and (SCO3392 or SCO4565) and SCO4566 and (SCO4567 or SCO6560) and SCO4568 and (SCO4569 or SCO4602) and (SCO4570 or SCO4603) and (SCO4571 or SCO4604) and (SCO4572 or SCO4605) and (SCO4573 or SCO4606 or SCO6954) and (SCO4574 or SCO4607) and (SCO4575 or SCO4608 or SCO6956)",
    "ATPM":     "bounds:(2.64, 1000)",
    "AMMQT9":   "bounds:(-1000,1000);genes:SCO4556 or SCO5940", 
    }



def add_reactions(iAA1259_model, scoGEM, reaction_mapping_fn, add_new_metabolites = True):
    apply_metabolite_mapping(iAA1259_model)

    # Get reaction_mapping
    reaction_mapping_df = pd.read_csv(reaction_mapping_fn, sep = ";")

    # Get new reactions
    new_reactions_id_list = list(reaction_mapping_df[reaction_mapping_df["Add"]]["Reaction ID"])

    # scoGEM reactions
    scoGEM_reaction_ids = [r.id for r in scoGEM.reactions]
    not_added_reactions = []
    N_r = len(scoGEM.reactions)
    N_m = len(scoGEM.metabolites)
    N_g = len(scoGEM.genes)
    i = 0
    for new_reaction_id in new_reactions_id_list:
        new_reaction = iAA1259_model.reactions.get_by_id(new_reaction_id)
        new_reaction.annotation["origin"] = "iAA1259"
        if new_reaction.id in scoGEM_reaction_ids:
            logging.info("Reaction {0} is already in scoGEM".format(new_reaction.id))
            # check_extra_reaction_annotations(scoGEM.reactions.get_by_id(new_reaction.id), new_reaction)        
        else:
            has_all, missing_metabolites, _ = check_metabolites(new_reaction, scoGEM, origin = "iAA1259")
            if has_all or add_new_metabolites:
                scoGEM.add_reaction(new_reaction)
                logging.info("Added reaction {0}:{1}".format(new_reaction.id, new_reaction.name))
                if len(missing_metabolites):
                    logging.info("\t and added new metabolites: {0}".format(", ".join([m.id for m in missing_metabolites])))
                i += 1
            else:
                print(new_reaction_id, [x.id for x in missing_metabolites])
                not_added_reactions.append(new_reaction_id)

    print("Added {0} reactions".format(i))
    print("Previously: {0} metabolites, {1} reactions, {2} genes".format(N_m, N_r, N_g))
    print("Now: {0} metabolites, {1} reactions, {2} genes".format(len(scoGEM.metabolites), len(scoGEM.reactions), len(scoGEM.genes)))
    return scoGEM

def modify_reactions(scoGEM):
    for reaction_id, change_string in CHANGED_REACTIONS.items():
        key, value = change_string.split(":")
        if key == "remove":
            scoGEM.reactions.get_by_id(reaction_id).remove_from_model()
        else:
            if key == "genes":
                scoGEM.reactions.get_by_id(reaction_id).gene_reaction_rule = value
            elif key == "bounds":
                bounds = literal_eval(value)
                scoGEM.reactions.get_by_id(reaction_id).bounds = bounds
            else:
                raise NotImplementedError
    return scoGEM

def change_biomass(iAA1259_model, scoGEM):
    iAA1259_biomass = iAA1259_model.reactions.Biomass

    changed_mets = {}
    for m, s in iAA1259_biomass.metabolites.items():
        try:
            coeff = scoGEM.reactions.BIOMASS_SCO.get_coefficient(m.id)
        except KeyError:
            coeff = 0

        if s != coeff:
            changed_mets[m] = s - coeff
            logging.info("{0}: {1}, {2}".format(m.id, s, coeff))
    
    scoGEM.reactions.BIOMASS_SCO.add_metabolites(changed_mets)
    scoGEM.reactions.BIOMASS_SCO_tRNA.add_metabolites(changed_mets)
    
    scoGEM.reactions.BIOMASS_SCO.annotation["origin"] = "iAA1259"
    scoGEM.reactions.BIOMASS_SCO.name = "S. coelicolor biomass objective function - with 75.79 GAM estimate"
    scoGEM.reactions.BIOMASS_SCO_tRNA.annotation["origin"] = "iAA1259"
    print("Changed biomass according to iAA1259")
    return scoGEM


def get_suppl_sheets():
    S4_df = pd.read_csv(S4_FN, sep = ";")

def fix_iAA1259(model):
    for r in model.reactions:
        r.id = r.id.replace("_LPAREN_e_RPAREN_", "_e").replace("_LPAREN_c_RPAREN_", "_c").replace("_DASH_", "__")
    for m in model.metabolites:
        m.id = m.id.replace("_LPAREN_e_RPAREN_", "_e").replace("_LPAREN_c_RPAREN_", "_c").replace("_DASH_", "__")


def apply_metabolite_mapping(iAA1259_model):
    for old_id, new_id in METABOLITE_MAPPINGS.items():
        m = iAA1259_model.metabolites.get_by_id(old_id)
        m.id = new_id

def map_S4():
    S4_df = pd.read_csv(S4_FN, sep = ";")
    print(S4_df)
    scoGEM = cobra.io.read_sbml_model(scoGEM_FN)
    iAA1259 = cobra.io.read_sbml_model(iAA1259_PATH)
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
    iAA1259_model = cobra.io.read_sbml_model(iAA1259_PATH)
    scoGEM = cobra.io.read_sbml_model(scoGEM_FN)
    fix_iAA1259(iAA1259_model)
    add_reactions(iAA1259_model, scoGEM, S4_FN, True)
    change_biomass(iAA1259_model, scoGEM)
    # map_S4()

    # with open(S4_FN, "r") as f:
    #     reader = csv.reader(f, delimiter = ";")
    #     for row in reader:
    #         print(len(row), row)
