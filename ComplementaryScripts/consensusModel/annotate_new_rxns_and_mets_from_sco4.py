import cobra
import pandas as pd
import requests
from bs4 import BeautifulSoup as BS
import logging

BIGG_REACTION_URL_BASE = "http://bigg.ucsd.edu/universal/reactions/"
BIGG_METABOLITE_URL_BASE = "http://bigg.ucsd.edu/universal/metabolites/"



def in_BiGG(BiGG_id, id_type = "Reaction"):

    # If last characheter is lower, check if the main reaction exists
    if id_type.lower() == "reaction":
        url_base = BIGG_REACTION_URL_BASE
        if BiGG_id[-1].islower():
            general_BiGG_id = BiGG_id[:-1]
            bigg_ids = [BiGG_id, general_BiGG_id]
        else:
            bigg_ids = [BiGG_id]
    else:
        url_base = BIGG_METABOLITE_URL_BASE
        bigg_ids = [BiGG_id]

    in_BiGG = False    
    for b_id in bigg_ids:
        url = url_base + BiGG_id.strip()
        html_text = requests.get(url).text
        soup = BS(html_text, "html.parser")

        if soup.title.contents[0][:3] == "404":
            pass
        else:
            in_BiGG = True
    return in_BiGG

def add_rxn_annotations(model, rxn_annotation_fn):
    reaction_annotation_df = pd.read_csv(rxn_annotation_fn, sep = ";", keep_default_na = False)
    
    for index, row in reaction_annotation_df.iterrows():
        # Check annotation
        try:
            reaction = model.reactions.get_by_id(row["Reaction ID"])
        except KeyError:
            logging.info("Can't find reaction {0} with proper id: {1}".format(row["Reaction ID"], row["Fixed ID"]))
            continue

        if not len(row["bigg"]):
            # The BiGG id is created and should not be an existing BiGG id
            # Check that it isn't used
            if in_BiGG(row["New ID"]):
                raise ValueError("""The BiGG id {0} is already defined in BiGG,
                 but it is defined as new in spreadsheet. Check that it is correct""".format(row["New ID"]))
        
        try:
            reaction.id = row["New ID"]
        except ValueError as e:
            print(reaction.id, end = "\t")
            print(e)
            continue

        
        if len(row["KEGG annotation"]):
            reaction.annotation["kegg.reaction"] = row["KEGG annotation"]
        if len(row["metanetx"]):
            reaction.annotation["metanetx.reaction"] = row["metanetx"]
        
        reaction.annotation["biocyc"] = row["Fixed ID"]
        reaction.name = row["Reaction name"]
        reaction.annotation["bigg.reaction"] = row["New ID"]
        logging.info("Changed reaction id from {0} to {1}".format(row["Reaction ID"], row["New ID"]))

def insert_ascii(reaction_id):
    return reaction_id.replace("-","__45__").replace("(","__40__").replace(")","__41__").replace(".", "__46__").replace("+", "__43__")

def add_met_annotations(model, met_annotation_fn):
    met_annotation_df = pd.read_csv(met_annotation_fn, sep = ";", keep_default_na = False)
    
    for index, row in met_annotation_df.iterrows():
        # Check annotation
        fixed_id = insert_ascii(row["Metabolite ID"])
        try:
            metabolite = model.metabolites.get_by_id(fixed_id)
        except KeyError:
            print("Can't find metabolite {1} with proper id: {0}".format(row["Metabolite ID"], fixed_id))
            continue

        if not len(row["BIGG"]):
            # The BiGG id is created and should not be an existing BiGG id
            # Check that it isn't used
            if in_BiGG(row["New ID"]):
                raise ValueError("""The BiGG id {0} is already defined in BiGG,
                 but it is defined as new in spreadsheet. Check that it is correct""".format(row["New ID"]))
        
        new_m_id = "{0}_{1}".format(row["New ID"], metabolite.id.rsplit("_", maxsplit = 1)[-1])
        try:
            metabolite.id = new_m_id
        except ValueError as e:
            print(metabolite.id, end = "\t")
            print(e)
            continue
        metabolite.annotation["biocyc"] = row["Metabolite ID"]
        metabolite.annotation["kegg.compound"] = row["KEGG"]
        metabolite.annotation["metanetx.chemical"] = row["MNX"]
        metabolite.name = row["Name"]
        metabolite.annotation["bigg.metabolite"] = row["New ID"]
        logging.info("Changed metabolite id from {0} to {1}".format(row["Metabolite ID"], row["New ID"]))


if __name__ == '__main__':
    fn = "../../ComplementaryData/curation/added_sco4_reactions.csv"
    model = cobra.io.read_sbml_model("../../ModelFiles/xml/scoGEM.xml")
    fnm = "../../ComplementaryData/curation/added_sco4_metabolites.csv"
    # add_rxn_annotations(model, fn)
    add_met_annotations(model, fnm)

    