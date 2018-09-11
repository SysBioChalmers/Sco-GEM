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


def add_reactions(sco4_model, scoGEM, reaction_mappping_fn, metabolite_mapping_fn):
    # Apply metabolite mapping to Sco4
    apply_metabolite_mapping(sco4_model, metabolite_mapping_fn)


    #Read reaction mapping
    reaction_mapping_df = pd.read_csv(reaction_mappping_fn, header = 0)

    # Get new reactions
    new_reaction_id_list = list(reaction_mappping_df[reaction_mappping_df["comment"] == "Sco4 specific"]["Sco4_ID"])

    for reaction_id in new_reaction_id_list:
        new_reaction = sco4_model.reactions.get_by_id(reaction_id)
        scoGEM.add_reaction(new_reaction)
        print("Added reaction {0}:{1}".format(new_reaction.id, new_reaction.name))
        
def apply_metabolite_mapping(sco4_model, metabolite_mappping_fn):
    # Read metabolite mapping
    metabolite_mapping_df = pd.read_csv(metabolite_mappping_fn, header = 0)

    # Apply metabolite mapping
    for index, row in metabolite_mapping_df.iterrows():
        if row["comment"] == "identical":
            if row["iKS1317_ID"] != row["Sco4_ID"]:
                sco4_met = sco4_model.metabolites.get_by_id(row["Sco4_ID"])
                sco4_met.id = row["iKS1317_ID"]
                print("Mapped {0} to {1}".format(row["Sco4_ID"], row["iKS1317_ID"]))


if __name__ == '__main__':
    SCO4_PATH = "../../ComplementaryData/models/Sco4.xml"
    sco4_model = cobra.io.read_sbml_model(SCO4_PATH)

    MODEL_PATH = "../../ModelFiles/xml/scoGEM.xml"
    scoGEM = cobra.io.read_sbml_model(MODEL_PATH)

    reaction_mappping_fn = "../../ComplementaryData/curation/rxns_iKS1317_vs_Sco4.csv"
    metabolite_mappping_fn = "../../ComplementaryData/curation/mets_iKS1317_vs_Sco4.csv"
    
    add_reactions(sco4_model, scoGEM, reaction_mappping_fn, metabolite_mappping_fn)

