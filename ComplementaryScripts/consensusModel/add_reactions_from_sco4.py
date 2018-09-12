# -*- coding: utf-8 -*-
"""
This file adds reactions from Sco4 based on the mapping of reactions and metabolites. 
Author: Snorre Sulheim
Date: 10.09.2018

# Description

# Files required

"""
import cobra
import pandas as pd

def insert_ascii(reaction_id):
    return reaction_id.replace("-","__45__").replace("(","__40__").replace(")","__41__").replace(".", "__46__").replace("+", "__43__")

def fix_compartment(met_id):
    return met_id.replace("[", "_").replace("]", "")
    

def add_reactions(sco4_model, scoGEM, reaction_mapping_fn, metabolite_mapping_fn):
    # Apply metabolite mapping to Sco4
    apply_metabolite_mapping(sco4_model, metabolite_mapping_fn)
    print(sco4_model.metabolites.atp_c)

    #Read reaction mapping
    reaction_mapping_df = pd.read_csv(reaction_mapping_fn, header = 0, sep = ";")
    # print(reaction_mapping_df)
    # Get new reactions
    new_reaction_id_list = list(reaction_mapping_df[reaction_mapping_df["Add"]]["Sco4_ID"])

    scoGEM_reaction_ids = [r.id for r in scoGEM.reactions]
    print(len(scoGEM.reactions))
    print(len(scoGEM.metabolites))
    all_missing = []
    for reaction_id in new_reaction_id_list:
        ascii_reaction_id = insert_ascii(reaction_id)
        new_reaction = sco4_model.reactions.get_by_id(ascii_reaction_id)
        if new_reaction.id in scoGEM_reaction_ids:
            print("Reaction {0} is already in scoGEM".format(new_reaction.id))
            continue

        has_all, missing = check_metabolites(new_reaction, scoGEM)
        all_missing += missing
        if has_all:
            scoGEM.add_reaction(new_reaction)
            print("Added reaction {0}:{1}".format(new_reaction.id, new_reaction.name))
    print(len(scoGEM.reactions))
    print(len(scoGEM.metabolites))
    print(set(all_missing))

        

def check_metabolites(reaction, scoGEM):
    scoGEM_mets = [x.id for x in scoGEM.metabolites]
    has_all = True
    missing = []
    for met, coeff in reaction.metabolites.items():
        if not met.id in scoGEM_mets:
            print("Missing ", met.id, "in reaction ", reaction.id)
            has_all = False
            missing.append(met.id)
    return has_all, missing

def apply_metabolite_mapping(sco4_model, metabolite_mapping_fn):
    # Read metabolite mapping
    metabolite_mapping_df = pd.read_csv(metabolite_mapping_fn, header = 0)

    # Apply metabolite mapping
    for index, row in metabolite_mapping_df.iterrows():
        if not pd.isnull(row["iKS1317_ID"]) and not pd.isnull(row["Sco4_ID"]):
            if row["iKS1317_ID"] != row["Sco4_ID"]:
                ascii_met_id = insert_ascii(row["Sco4_ID"])
                try:
                    sco4_model.metabolites.get_by_id(ascii_met_id).id = fix_compartment(row["iKS1317_ID"])
                except ValueError:
                    replace_metabolite(ascii_met_id, fix_compartment(row["iKS1317_ID"]), sco4_model)

                    # print("Mapped {0} to {1}".format(row["Sco4_ID"], row["iKS1317_ID"])

def replace_metabolite(old_id, new_id, sco4_model):
    old_m = sco4_model.metabolites.get_by_id(old_id)
    new_m = sco4_model.metabolites.get_by_id(new_id)
    for r in old_m.reactions:
        coeff = r.get_coefficient(old_id)
        r.add_metabolites({old_m: -coeff, new_m: coeff})


if __name__ == '__main__':
    SCO4_PATH = "../../ComplementaryData/models/Sco4.xml"
    sco4_model = cobra.io.read_sbml_model(SCO4_PATH)

    MODEL_PATH = "../../ModelFiles/xml/scoGEM.xml"
    scoGEM = cobra.io.read_sbml_model(MODEL_PATH)

    reaction_mapping_fn = "../../ComplementaryData/curation/rxns_iKS1317_vs_Sco4.csv"
    metabolite_mapping_fn = "../../ComplementaryData/curation/mets_iKS1317_vs_Sco4.csv"
    
    add_reactions(sco4_model, scoGEM, reaction_mapping_fn, metabolite_mapping_fn)

