#!/usr/bin/env python
"""
Author: Snorre Sulheim

# Description
A variety fo functions not used in the pipeline directly, but to cureate a draft mapping
of reactions and metabolites in Sco4 to iKS1317
"""
import sys
import os
import cobra
import re
import numpy as np
import pandas as pd
from collections import defaultdict
import fix_sco4_issues
import fix_iKS1317_issues



# GLOBAL CONSTANTS
SCO4_PATH = r"../../ComplementaryData/models/Sco4.xml"
iKS1317_PATH = r"../../ModelFiles/xml/scoGEM.xml"

ORGANISM = "sco" # Streptomyces coelicolor
CHANGED_REACTIONS_FN = r"../..//ComplementaryData/curation/iKS1317_changed_reaction_ids.csv"

def load_models(sco4_path, iKS1317_path):
    model_sco4 = cobra.io.read_sbml_model(sco4_path)
    model_iKS1317 = cobra.io.read_sbml_model(iKS1317_path)
    return model_sco4, model_iKS1317

def replace_ascii_codes(new_reactions):
    ascii_free_id_dict = {}
    for r in new_reactions:
        ascii_free_id = re.sub(r"__(\d{2})__(\d{2})__(\d{2})__", lambda x: chr(int(x.group(1))), r.id)
        ascii_free_id = re.sub(r"__(\d{2})__(\d{2})__", lambda x: chr(int(x.group(1))), r.id)
        ascii_free_id = re.sub(r"__(\d{2})__", lambda x: chr(int(x.group(1))), r.id)
        # ascii_free_id = ascii_free_id.lstrip("__")
        ascii_free_id = re.sub(r"RXN__", "RXN", ascii_free_id)
        ascii_free_id = re.sub(r"CPD__", "CPD", ascii_free_id)
        ascii_free_id_dict[r.id] =  ascii_free_id
        # print(r.id, "\t", ascii_free_id)
    return ascii_free_id_dict

def discard_reactions_with_identical_id(model_iKS1317, model_sco4):
    iKS1317_reaction_ids = [x.id for x in model_iKS1317.reactions]
    new_reactions = []
    for r in model_sco4.reactions:
        if r.id in iKS1317_reaction_ids:
            pass
        else:
            new_reactions.append(r)
    return new_reactions

def discard_reactions_based_on_biocyc_annotation(reaction_list, model_iKS1317):
    for r in model_iKS1317.reactions:
        try:
            biocyc_id = r.annotation["biocyc"]
        except:
            continue
        else:
            if 0:
                # if biocyc_id in 
                pass

def make_reaction_dataframe(new_reactions):
    reaction_df = pd.DataFrame([(str(r.id), r, str(r.name)) for r in new_reactions], columns = ["Reaction ID", "Reaction", "Name"])
    return reaction_df

def make_metabolite_dataframe(new_metabolites):
    metabolite_df = pd.DataFrame([(str(r.id), r, str(r.name)) for r in new_metabolites], columns = ["metabolite ID", "metabolite", "Name"])
    return metabolite_df


def add_metabolite_annotations(new_metabolites):
    database_dicts = [({}, "kegg.compound"), ({}, "biocyc"), ({}, "chebi"), ({}, "pubchem.compound")]
    for m in new_metabolites:
        for dct, dct_key in database_dicts:
            try:
                db_id = m.annotation[dct_key]
            except KeyError:
                pass
            else:
                if dct_key == "chebi":
                    db_id = db_id.replace("CHEBI:", "")
                dct[m.id] = db_id
    return [x for x, i in database_dicts]

def add_reaction_annotations(new_reactions):
    ### BIOCYC
    kegg_dict = {}
    rhea_dict = {}
    biocyc_dict = {}

    for r in new_reactions:
        try:
            rhea = r.annotation["rhea"]
        except KeyError:
            pass
        else:
            rhea_dict[r.id]= rhea

        try:
            kegg = r.annotation["kegg.reaction"]
        except KeyError:
            # print(r.annotation)
            pass
        else:
            kegg_dict[r.id] = str(kegg)

        try:
            biocyc_id = r.annotation["biocyc"]
        except KeyError:
            pass
        else:
            biocyc_dict[r.id] = biocyc_id

    return kegg_dict, rhea_dict, biocyc_dict

def map_metabolites_to_iKS1317(metabolite_df, model_iKS1317):
    mapping_dict = defaultdict(list)
    reason_dict = defaultdict(list)
    unmatched_metabolites = [r.id for r in model_iKS1317.metabolites]
    for (index, row) in metabolite_df.iterrows():
        #Map based on identical id
        if row["metabolite ID"] in unmatched_metabolites:
            unmatched_metabolites.remove(row["metabolite ID"])
            mapping_dict[index].append(row["metabolite ID"])
            reason_dict[index].append("Identical ID")

        # Map based on converted ID
        elif row["Converted ID"] in unmatched_metabolites:
            unmatched_metabolites.remove(row["Converted ID"])
            mapping_dict[index].append(row["Converted ID"])
            reason_dict[index].append("Converted ID")
        
    # Map based on KEGG        
    for m in model_iKS1317.metabolites:
        if m.id in unmatched_metabolites:
            try:
                match = metabolite_df["KEGG ID"] == m.annotation["kegg.compound"]
            except KeyError:
                continue
            else:
                if match.any():
                    kegg_match = metabolite_df.index[match].tolist()
                    for x in kegg_match:
                        if metabolite_df.loc[x]["Converted ID"].rsplit("_")[1] == m.compartment: 
                            reason_dict[x].append("KEGG")
                            mapping_dict[x].append(m.id)

    # Map based on biocyc
    for r in model_iKS1317.metabolites:
        if r.id in unmatched_metabolites:
            try:
                match = metabolite_df["biocyc ID"] == r.annotation["biocyc"]
            except KeyError:
                continue
            else:
                if match.any():
                    print("### {0} ####".format(r.id))
                    biocyc_match = metabolite_df.index[match].tolist()
                    for x in biocyc_match:
                        if metabolite_df.loc[x]["Converted ID"].rsplit("_")[1] == m.compartment: 
                            mapping_dict[x].append(r.id)
                            reason_dict[x].append("biocyc")
     # Map based on chebi
    for m in model_iKS1317.metabolites:
        if m.id in unmatched_metabolites:
            try:
                match = metabolite_df["chebi ID"] == m.annotation["chebi"]
            except KeyError:
                continue
            else:
                if match.any():
                    print("### {0} ####".format(m.id))
                    chebi_match = metabolite_df.index[match].tolist()
                    for x in chebi_match:
                        if metabolite_df.loc[x]["Converted ID"].rsplit("_")[1] == m.compartment: 
                            mapping_dict[x].append(m.id)
                            reason_dict[x].append("chebi")

    metabolite_df["Matched metabolites"] = metabolite_df.index.map(mapping_dict.get)
    metabolite_df["Reason"] = metabolite_df.index.map(reason_dict.get)
    return metabolite_df, unmatched_metabolites



def map_reactions_to_iKS1317(reaction_df, model_iKS1317):
    mapping_dict = defaultdict(list)
    reason_dict = defaultdict(list)

    unmatched_reactions = [r.id.lower() for r in model_iKS1317.reactions]

    change_reaction_ids = pd.read_csv(CHANGED_REACTIONS_FN, header = None, sep = ",", index_col = 0).to_dict()[1]

    for (index, row) in reaction_df.iterrows():
        #Map based on identical id
        if row["Reaction ID"].lower() in unmatched_reactions:
            unmatched_reactions.remove(row["Reaction ID"].lower())
            mapping_dict[index].append(row["Reaction ID"])
            reason_dict[index].append("Identical ID")

        # Map based on converted ID
        elif row["Converted ID"].lower() in unmatched_reactions:
            unmatched_reactions.remove(row["Converted ID"].lower())
            mapping_dict[index].append(row["Converted ID"])
            reason_dict[index].append("Converted ID")
        
        # Map based on changed reaction ids
        elif row["Converted ID"] in change_reaction_ids.keys():
            unmatched_reactions.remove(change_reaction_ids[row["Converted ID"]].lower())
            mapping_dict[index].append(change_reaction_ids[row["Converted ID"]])
            reason_dict[index].append("Changed ID")


    # Map based on biocyc
    for r in model_iKS1317.reactions:
        if r.id.lower() in unmatched_reactions:
            try:
                match = reaction_df["Ascii free ID"] == r.annotation["biocyc"]
            except KeyError:
                continue
            else:
                if match.any():
                    print("### {0} ####".format(r.id))
                    biocyc_match = reaction_df.index[match].tolist()
                    for x in biocyc_match:
                        mapping_dict[x].append(r.id)
                        reason_dict[x].append("biocyc")

    # Map based on KEGG        
    for r in model_iKS1317.reactions:
        if r.id.lower() in unmatched_reactions:
            try:
                match = reaction_df["KEGG ID"] == r.annotation["kegg.reaction"]
            except KeyError:
                continue
            else:
                if match.any():
                    kegg_match = reaction_df.index[match].tolist()
                    for x in kegg_match:
                        reason_dict[x].append("KEGG")
                        mapping_dict[x].append(r.id)
  
    reaction_df["Matched reactions"] = reaction_df.index.map(mapping_dict.get)
    reaction_df["Add"] = [not x for x in reaction_df["Matched reactions"]]
    reaction_df["Reason"] = reaction_df.index.map(reason_dict.get)
    return reaction_df, unmatched_reactions

def map_to_biocyc(reaction_df):
    # Map based on biocyc
    try:
        biocyc = pythoncyc.get_organism(ORGANISM)
    except:
        print("Pythoncyc doesn't work. Have you started pathwaytools? \n")
    else:
        for index, row in reaction_df.iterrows():
            biocyc_match = biocyc[row["Ascii free ID"]]

def add_biocyc_annotations(model_sco4):
    for m in model_sco4.metabolites:
        biocyc = re.findall(r"CPD__45__\d{3,5}", m.id)
        if len(biocyc) == 0:
            pass
        elif len(biocyc) == 1:
            m.annotation["biocyc"] = biocyc[0].replace("__45__", "-")
        else:
            print("2 kegg annotations in ID: ", biocyc)

    return model_sco4

def add_kegg_annotations(model_sco4):
    for r in model_sco4.reactions:
        kegg_id = re.findall(r"R\d{5}", r.id)
        if len(kegg_id) == 0:
            pass
        elif len(kegg_id) == 1:
            r.annotation["kegg.reaction"] = kegg_id[0]
        else:
            print("2 kegg annotations in ID: ", kegg_id)

    for m in model_sco4.metabolites:
        kegg_id = re.findall(r"C\d{5}", m.id)
        if len(kegg_id) == 0:
            pass
        elif len(kegg_id) == 1:
            m.annotation["kegg.compound"] = kegg_id[0]
        else:
            print("2 kegg annotations in ID: ", kegg_id)

    return model_sco4

def run_reaction_mapping_to_iKS1317(apply_iKS1317_fix = False, apply_sco4_fix = False):
    model_sco4, model_iKS1317 = load_models(SCO4_PATH, iKS1317_PATH)

    # Apply model changes
    if apply_iKS1317_fix:
        fix_iKS1317_issues.fix(model_iKS1317)
    if apply_sco4_fix:
        fix_sco4_issues.fix(model_sco4)
    

    # new_reactions = discard_reactions_with_identical_id(model_iKS1317, model_sco4)
    model_sco4 = add_kegg_annotations(model_sco4)
    # reaction_df = make_dataframe(new_reactions)
    
    ### ???
    reaction_df = make_reaction_dataframe(list(model_sco4.reactions))
    # print(reaction_df.head())
    # reaction_ascii_dict = replace_ascii_codes(new_reactions)
    reaction_ascii_dict = replace_ascii_codes(model_sco4.reactions)

    reaction_df["Ascii free ID"] = reaction_df["Reaction ID"].map(reaction_ascii_dict)
    reaction_df["Converted ID"] = [x.replace("(e)", "_e").replace("(c)", "_c").replace("-L", "__L").replace("-D", "__D").replace("-R", "__R") for x in reaction_df["Ascii free ID"]]
    kegg_dict, rhea_dict, biocyc_dict = add_reaction_annotations(list(model_sco4.reactions))#new_reactions)
    reaction_df["KEGG ID"] = reaction_df["Reaction ID"].map(kegg_dict)
    reaction_df["Rhea ID"] = reaction_df["Reaction ID"].map(rhea_dict)

    mapped_reaction_df, unmatched_reactions = map_reactions_to_iKS1317(reaction_df, model_iKS1317)
    mapped_reaction_df["Matched reactions"] = mapped_reaction_df["Matched reactions"].apply(list_to_str)
    mapped_reaction_df["Reason"] = mapped_reaction_df["Reason"].apply(list_to_str)
    mapped_reaction_df.to_csv(os.path.join(os.getcwd(), "mapping.csv"), columns = ["Reaction ID", "KEGG ID", "Matched reactions", "Reason", "Add"], index = False, sep = "\t")
    print(unmatched_reactions)    


def run_metabolite_mapping_to_iKS1317():
    model_sco4, model_iKS1317 = load_models(SCO4_PATH, iKS1317_PATH)
    # new_metabolites = discard_metabolites_with_identical_id(model_sco4, model_iKS1317)
    add_kegg_annotations(model_sco4)
    add_biocyc_annotations(model_sco4)
    metabolites = [x for x in model_sco4.metabolites]
    metabolite_df = make_metabolite_dataframe(metabolites)
    metabolite_ascii_dict = replace_ascii_codes(model_sco4.metabolites)
    metabolite_df["Ascii free ID"] = metabolite_df["metabolite ID"].map(metabolite_ascii_dict)
    metabolite_df["Converted ID"] = [x.replace("-D", "__D").replace("-R", "__R").replace("-L", "__L").replace("-S", "__S").replace("-", "__").replace("CPD__", "CPD-") for x in metabolite_df["Ascii free ID"]]
    # metabolite_df["Converted ID"] = [x.replace("-B", "_B").replace("-C", "_C").replace("-", "__").replace("CPD__", "CPD-") for x in metabolite_df["Ascii free ID"]]
    kegg_dict, biocyc_dict, chebi_dict, pubchem_dict = add_metabolite_annotations(metabolites)
    metabolite_df["KEGG ID"] = metabolite_df["metabolite ID"].map(kegg_dict)
    metabolite_df["biocyc ID"] = metabolite_df["metabolite ID"].map(biocyc_dict)
    metabolite_df["chebi ID"] = metabolite_df["metabolite ID"].map(chebi_dict)
    metabolite_df["pubchem ID"] = metabolite_df["metabolite ID"].map(pubchem_dict)

    mapped_metabolite_df, unmapped_mets = map_metabolites_to_iKS1317(metabolite_df, model_iKS1317)
    mapped_metabolite_df["Matched metabolites"] = mapped_metabolite_df["Matched metabolites"].apply(list_to_str)
    mapped_metabolite_df["Reason"] = mapped_metabolite_df["Reason"].apply(list_to_str)
    mapped_metabolite_df.to_csv(os.path.join(os.getcwd(), "metabolite_mapping.csv"), 
        columns = ["metabolite ID", "Converted ID", "KEGG ID", "biocyc ID", "chebi ID", "Matched metabolites", "Reason"], index = False, sep = "\t")
    

    # print(len(new_metabolites))
    # print([x.id for x in new_metabolites], sep = "\n")

def discard_metabolites_with_identical_id(model_sco4, model_iKS1317):
    new_metabolites = []
    iKS1317_metabolite_ids = [m.id for m in model_iKS1317.metabolites]

    for m in model_sco4.metabolites:
        if m.id in iKS1317_metabolite_ids:
            try:
                m_kegg = m.annotation["kegg.compound"]
            except KeyError:
                pass
            else:
                iKS1317_metabolite = model_iKS1317.metabolites.get_by_id(m.id)
                try:
                    m_iKS1317_kegg = iKS1317_metabolite.annotation["kegg.compound"]
                except KeyError:
                    pass
                else:
                    if isinstance(m_kegg, str) and len(m_kegg):
                        if isinstance(m_iKS1317_kegg, str) and len(m_iKS1317_kegg):
                            if m_iKS1317_kegg != m_kegg:
                                print("Same metabolite id, different KEGG: ", m.id, m_kegg, m_iKS1317_kegg)
                                
        else:
            new_metabolites.append(m)
    return new_metabolites


def list_to_str(lst):
    if isinstance(lst, list):
        return ", ".join(set(lst))
    else:
        return ""



def check_biocyc():
    model_sco4, model_iKS1317 = load_models(SCO4_PATH, iKS1317_PATH)
    add_biocyc_annotations(model_sco4)
    biocyc_sco4 = []
    biocyc_iks1317 = []
    for m in model_sco4.metabolites:
        try:
            biocyc = m.annotation["biocyc"]
        except:
            continue
        else:
            # biocyc = biocyc.replace("biocyc:", "")
            biocyc_sco4.append(biocyc)

    for m in model_iKS1317.metabolites:
        try:
            biocyc = m.annotation["biocyc"]
        except:
            continue
        else:
            # biocyc = biocyc.replace("biocyc:", "")
            biocyc_iks1317.append(biocyc)

    print(set(biocyc_iks1317).intersection(biocyc_sco4))
    print(biocyc_iks1317[:10], len(biocyc_iks1317))
    print(biocyc_sco4[:10], len(biocyc_sco4))

    m_ids = [m.id.replace("_c", "").replace("__45__", "-") for m in model_sco4.metabolites]
    print(set(biocyc_iks1317).intersection(m_ids))

def check_chebi():
    model_sco4, model_iKS1317 = load_models(SCO4_PATH, iKS1317_PATH)
    chebi_sco4 = []
    chebi_iks1317 = []
    for m in model_sco4.metabolites:
        try:
            chebi = m.annotation["chebi"]
        except:
            continue
        else:
            chebi = chebi.replace("CHEBI:", "")
            chebi_sco4.append(chebi)

    for m in model_iKS1317.metabolites:
        try:
            chebi = m.annotation["chebi"]
        except:
            continue
        else:
            # chebi = chebi.replace("CHEBI:", "")
            chebi_iks1317.append(chebi)

    print(set(chebi_iks1317).intersection(chebi_sco4))
    print(chebi_iks1317[:10], len(chebi_iks1317))
    print(chebi_sco4[:10], len(chebi_sco4))

def update_scoGEM_mapping():
    model_sco4, model_iKS1317 = load_models(SCO4_PATH, iKS1317_PATH)

    scoGEM_csv = r"..\..\ComplementaryData\curation\rxns_iKS1317_vs_Sco4.csv"
    df_scoGEM = pd.read_csv(scoGEM_csv, sep = ";")#, keep_default_na=False)
    df_scoGEM["Add"] = False
    
    df_mapping = pd.read_csv(r"../../ComplementaryData/curation/Sco4_reaction_mapping.csv", sep = ";")
    print(df_mapping.columns.values)
    print(df_scoGEM.columns.values)
    reaction_ascii_dict = replace_ascii_codes(model_sco4.reactions)
    df_mapping["Converted ID"] = df_mapping["Reaction ID"].map(reaction_ascii_dict)


    droplist_index = []
    droplist_iKS1317_ids = []
    g = 0
    n = 0
    j = 0
    for index, row in df_scoGEM.iterrows():
        if row["comment"] == "Sco4 specific":
            # print(row)
            matched = df_mapping[df_mapping["Converted ID"] == row["Sco4_ID"]]
            # print(matched)
            if matched.empty:
                droplist_index.append(index)
                print("Could not find ", row["Sco4_ID"])
            else:
                if matched.shape[0] == 1:
                    matched = matched.iloc[0]
                else:
                    print("2 matches!")
                    break

                if pd.isnull(matched["Matched reactions"]):
                    df_scoGEM.loc[index, "Add"] = matched["Add"]
                    grr = model_sco4.reactions.get_by_id(matched["Reaction ID"]).gene_reaction_rule
                    if not len(grr) and row["Sco4_ID"][:2] != "EX":
                        g += 1
                    # print(grr, end = "\t")
                    # print("No match, Add {0}: {1}".format(row["Sco4_ID"], matched["Add"]))
                else:
                    j += 1
                    df_scoGEM.loc[index, "iKS1317_ID"] = matched["Matched reactions"]
                    droplist_iKS1317_ids.append(str(matched["Matched reactions"]))
                    matched["Comment"] = matched["Comment"].replace(", ", " and ")
                    df_scoGEM.loc[index, "comment"] = matched["Comment"]
                    df_scoGEM.loc[index, "Add"] = matched["Add"]
                    if df_scoGEM.loc[index,"Add"]:
                        print("Matched {0} with {1}, Add: {2}".format(row["Sco4_ID"], row["iKS1317_ID"], row["Add"]))
                if df_scoGEM.loc[index,"Add"]:
                    n += 1

    
    for index, row in df_scoGEM.iterrows():
        if row["comment"] == "iKS1317 specific":
            if row["iKS1317_ID"] in droplist_iKS1317_ids:
                droplist_index.append(index)
                print("Dropping", row["iKS1317_ID"])
    df_scoGEM.drop(droplist_index, inplace = True)

    # Drop matched reactions
    print(g, n)
    df_scoGEM.to_csv(scoGEM_csv, columns = ["iKS1317_ID", "Sco4_ID", "comment", "Add"], index = False, sep = ",")


def print_missing_genes():
    model_sco4 = cobra.io.read_sbml_model(SCO4_PATH)
    reaction_ascii_dict = replace_ascii_codes(model_sco4.reactions)
    scoGEM_csv = r"..\..\ComplementaryData\curation\rxns_iKS1317_vs_Sco4.csv"
    df_scoGEM = pd.read_csv(scoGEM_csv, sep = ";")
    df_scoGEM = df_scoGEM.replace("TRUE", True)
    df_scoGEM = df_scoGEM.replace("FALSE", False)
    df_mapping = pd.read_csv(r"../../ComplementaryData/curation/Sco4_reaction_mapping.csv", sep = ";")
    df_mapping["Converted ID"] = df_mapping["Reaction ID"].map(reaction_ascii_dict)
    # print(df_mapping)
    g = 0
    for index, row in df_scoGEM.iterrows():
        if row["Add"] == True:
            matched = df_mapping[df_mapping["Converted ID"] == row["Sco4_ID"]].iloc[0]
            grr = model_sco4.reactions.get_by_id(matched["Reaction ID"]).gene_reaction_rule
            if not len(grr) and row["Sco4_ID"][:2] != "EX":
                print(row["Sco4_ID"])
                g += 1
    print(g)

    
if __name__ == '__main__':
    # run_reaction_mapping_to_iKS1317(apply_sco4_fix = True)
    # update_scoGEM_mapping()
    print_missing_genes()
    # run_metabolite_mapping_to_iKS1317()
    # check_biocyc()

