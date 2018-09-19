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
from collections import defaultdict
import re
import logging

def add_reactions(sco4_model, scoGEM, reaction_mapping_fn, metabolite_mapping_fn, add_new_metabolites = True):
    # Apply metabolite mapping to Sco4
    apply_metabolite_mapping(sco4_model, metabolite_mapping_fn)

    #Read reaction mapping
    reaction_mapping_df = pd.read_csv(reaction_mapping_fn, header = 0, sep = ";")

    # Get new reactions
    new_reactions_id_list = list(reaction_mapping_df[reaction_mapping_df["Add"]]["Sco4_ID"])

    scoGEM_reaction_ids = [r.id for r in scoGEM.reactions]
    N_r = len(scoGEM.reactions)
    N_m = len(scoGEM.metabolites)
    N_g = len(scoGEM.genes)

    all_missing = []
    missing_reaction_dict = defaultdict(list)
    not_added_reactions = []
    i = 0
    for reaction_id in new_reactions_id_list:
        ascii_reaction_id = insert_ascii(reaction_id)
        new_reaction = sco4_model.reactions.get_by_id(ascii_reaction_id)
        new_reaction.annotation["origin"] = "Sco4"
        if new_reaction.id in scoGEM_reaction_ids:
            logging.info("Reaction {0} is already in scoGEM".format(new_reaction.id))
            check_extra_reaction_annotations(scoGEM.reactions.get_by_id(new_reaction.id), new_reaction)
        else:
            has_all, missing_metabolites, missing_reactions = check_metabolites(new_reaction, scoGEM)
            all_missing += missing_metabolites
            for key, value in missing_reactions.items():
                missing_reaction_dict[key] += value
            if has_all or add_new_metabolites:
                scoGEM.add_reaction(new_reaction)
                logging.info("Added reaction {0}:{1}".format(new_reaction.id, new_reaction.name))
                if len(missing_metabolites):
                    logging.info("\t and added new metabolites: {0}".format(", ".join([m.id for m in missing_metabolites])))
                i += 1
            else:
                not_added_reactions.append(reaction_id)

    # evaluate_missing_reactions(missing_reaction_dict, reaction_mapping_df)
    print("Added {0} reactions".format(i))
    print("Previously: {0} metabolites, {1} reactions, {2} genes".format(N_m, N_r, N_g))
    print("Now: {0} metabolites, {1} reactions, {2} genes".format(len(scoGEM.metabolites), len(scoGEM.reactions), len(scoGEM.genes)))
    return scoGEM




def insert_ascii(reaction_id):
    return reaction_id.replace("-","__45__").replace("(","__40__").replace(")","__41__").replace(".", "__46__").replace("+", "__43__")

def replace_ascii(reaction_id):
    return reaction_id.replace("__45__", "-").replace("__40__", "(").replace("__41__", ")").replace( "__46__", ".").replace("__43__", "+")

def add_kegg_annotations(model_sco4):
    for r in model_sco4.reactions:
        kegg_id = re.findall(r"R\d{5}", r.id)
        if len(kegg_id) == 0:
            pass
        elif len(kegg_id) == 1:
            r.annotation["kegg.reaction"] = kegg_id[0]
        else:
            print("2 kegg annotations in ID: ", kegg_id)

def fix_compartment(met_id):
    return met_id.replace("[", "_").replace("]", "")
    
def print_new_reactions(new_reactions_id_list, sco4_model):
    all_new_list = []
    add_kegg_annotations(sco4_model)
    for reaction_id in new_reactions_id_list:
        temp = []
        ascii_reaction_id = insert_ascii(reaction_id)
        new_reaction = sco4_model.reactions.get_by_id(ascii_reaction_id)
        try:
            kegg = new_reaction.annotation["kegg.reaction"]
        except:
            kegg = ""
        temp = [new_reaction.id, new_reaction.name, kegg]
        all_new_list.append(temp)
    df = pd.DataFrame(all_new_list, columns = ["Reaction ID", "Reaction name", "KEGG annotation"])
    print(df)
    df.to_csv("new_reactions.csv", index = False, sep = ";")



def check_extra_reaction_annotations(scoGEM_reaction, sco4_reaction, print_new_annotations = False):
    # Genes
    scoGEM_genes = [g.id for g in scoGEM_reaction.genes]
    for gene in sco4_reaction.genes:
        if not gene.id in scoGEM_genes:
            print("Reaction {0} miss gene annotation {1}".format(sco4_reaction.id, gene.id))
    
    # Annotations
    for key, value in sco4_reaction.annotation.items():
        try:
            scoGEM_reaction.annotation[key]
        except KeyError:
            scoGEM_reaction.annotation[key] = value
            logging.info("Added annotation {0}: {1} to {2}".format(key, value, scoGEM_reaction.id))
        else:
            sco_rxn_anno = scoGEM_reaction.annotation[key]
            if isinstance(sco_rxn_anno, list):
                if not value in sco_rxn_anno:
                    scoGEM_reaction.annotation[key] = sco_rxn_anno.append(value)
                    logging.info("Appended annotation {0} to {1}".format(value, scoGEM_reaction.id))
            else:
                if value != sco_rxn_anno:
                    scoGEM_reaction.annotation[key] = [sco_rxn_anno, value]  
                    logging.info("Appended annotation {0} to {1}".format(value, scoGEM_reaction.id))

def evaluate_missing_reactions(missing_reaction_dict, reaction_mapping_df):
    for key, value in missing_reaction_dict.items():
        matched_reactions = []
        for r_id in value:
            ascii_free = replace_ascii(r_id)
            matched = reaction_mapping_df[reaction_mapping_df["Sco4_ID"] == ascii_free]
            if matched.shape[0] == 1:
                matched = matched.iloc[0]
            else:
                # print("Missing ", ascii_free)   
                continue
            if not pd.isnull(matched["iKS1317_ID"]):
                matched_reactions.append(matched["iKS1317_ID"])
        
        # print(key, list(set(matched_reactions)))

      

        
def add_met_annotations(scoGEM_met, met):
    for key, value in met.annotation.items():
        try:
            scoGEM_met.annotation[key]
        except KeyError:
            scoGEM_met.annotation[key] = value
            logging.info("Added annotation {0}: {1} to {2}".format(key, value, met.id))
        else:
            sco_met_anno = scoGEM_met.annotation[key]
            if isinstance(sco_met_anno, list):
                if not value in sco_met_anno:
                    scoGEM_met.annotation[key] = sco_met_anno.append(value)
                    logging.info("Appended annotation {0} to {1}".format(value, met.id))
            else:
                if value != sco_met_anno:
                    scoGEM_met.annotation[key] = [sco_met_anno, value]  
                    logging.info("Appended annotation {0} to {1}".format(value, met.id))
                  
def check_metabolites(reaction, scoGEM):
    scoGEM_mets = [x.id for x in scoGEM.metabolites]
    has_all = True
    missing = []
    missing_reactions = {}
    for met, coeff in reaction.metabolites.items():

        if met.id in scoGEM_mets:
            add_met_annotations(scoGEM.metabolites.get_by_id(met.id), met)
        else:
            logging.info("Missing {0} in reaction {1}".format(met.id, reaction.id))
            has_all = False
            met.annotation["origin"] = "Sco4"
            missing.append(met)
            missing_reactions[met.id] = [r.id for r in met.reactions if r.id != reaction.id]

    return has_all, missing, missing_reactions

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
    

    logging.basicConfig(filename='add_reactions_from_sco4.log', level=logging.INFO)
    scoGEM = add_reactions(sco4_model, scoGEM, reaction_mapping_fn, metabolite_mapping_fn)
    cobra.io.write_sbml_model(scoGEM, MODEL_PATH)

