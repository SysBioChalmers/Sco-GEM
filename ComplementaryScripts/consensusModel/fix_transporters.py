# -*- coding: utf-8 -*-
"""
This files applies changes to transport reactions in Sco-GEM suggested by Tjasa Kumelj.
The changes include:
1) Change gene annotation of existing reactions
2) Add new transport reactions to existing compounds
3) Add new transport and exchange reactions to new compounds
4))

The issue is described in the corresponding pull request (#70):
https://github.com/SysBioChalmers/Sco-GEM/pull/70

Author: Snorre Sulheim
Created: 08.07.2018

"""

import cobra
import pandas as pd
import logging

def fix_transporters(model, existing_reactions_fn, new_transport_reactions_fn, 
                     new_transport_reactions_new_mets_fn, new_metabolites_fn):
    modify_existing_transporters(model, existing_reactions_fn)
    new_transport_reactions_to_existing_metabolites(model, new_transport_reactions_fn)
    new_transport_reactions_to_new_metabolites(model, new_transport_reactions_new_mets_fn, new_metabolites_fn)

def modify_existing_transporters(model, existing_reactions_fn):
    df = pd.read_csv(existing_reactions_fn, header = 0, sep = ",")
    for (i, row) in df.iterrows():
        r_id = row["RxnID"].strip()
        try:
            r = model.reactions.get_by_id(r_id)
        except KeyError:
            logging.info("Reaction {0} is not in the model".format(r_id))
            continue
        r.gene_reaction_rule = row["Gene association"]

def new_transport_reactions_to_existing_metabolites(model, new_transport_reactions_fn):
    df = pd.read_csv(new_transport_reactions_fn, header = 0, sep = ",", index_col = False)
    for i, row in df.iterrows():
        r_id = row["RxnID"].strip()
        new_reaction = cobra.Reaction(r_id)
        new_reaction.gene_reaction_rule = row["Gene association"]
        model.add_reaction(new_reaction)
        reaction_string = row["Rxn_Formula"].replace("[","_").replace("]","")
        new_reaction.build_reaction_from_string(reaction_string)
        new_reaction.annotation["SBO"] = "SBO:0000655"
        new_reaction.annotation["subsystem"] = row["Subsystem"]
        new_reaction.name = row["Reaction name"]

        _check_and_fix_new_extracellular_metabolites(model, new_reaction)

def new_transport_reactions_to_new_metabolites(model, new_transport_reactions_new_mets_fn, new_mets_fn):
    new_mets_df = pd.read_csv(new_mets_fn, index_col = False, sep = ",", header = 0)
    df = pd.read_csv(new_transport_reactions_new_mets_fn, index_col = False, sep = ",", header = 0, quotechar = '"')

    # Add all new metabolites
    new_metabolites = []
    for i, row in new_mets_df.iterrows():
        for compartment in ["c", "e"]:
            m_id = "{0}_{1}".format(row["metabolite ID"], compartment)
            try:
                model.metabolites.get_by_id(m_id)
            except KeyError:
                m = cobra.Metabolite(m_id)
                m.name = row["KEGG_name"]
                m.compartment = compartment
                m.formula = row["chemical formula"]
                m.annotation["kegg.compound"] = row["KEGG ID"]
                m.annotation["SBO"] = "SBO:0000247"
                m.annotation["metabetx.chemical"] = row["metanetx ID"]
                m.charge = 0
                new_metabolites.append(m)
            else:
                logging.info("Metabolite {0} is alredy in the model".format(m_id))
                continue
            
    model.add_metabolites(new_metabolites)

    # Add new reactions
    for i, row in df.iterrows():
        r_id = row["RxnID"].strip()
        new_reaction = cobra.Reaction(r_id)
        new_reaction.gene_reaction_rule = row["Gene annotation"]
        model.add_reaction(new_reaction)
        reaction_string = row["Rxn_Formula"].replace("[","_").replace("]","")
        new_reaction.build_reaction_from_string(reaction_string)
        new_reaction.annotation["SBO"] = "SBO:0000655"
        new_reaction.annotation["subsystem"] = row["Subsystem"]
        new_reaction.name = row["Reaction name"]        
        _check_and_fix_new_extracellular_metabolites(model, new_reaction)






def _check_and_fix_new_extracellular_metabolites(model, reaction):
    for m, i in reaction.metabolites.items():
        # Check if this metabolite was added now and is an extracellular metabolite
        if not m.compartment:
            logging.info("Fixing metabolite {0}".format(m.id))

            # Find information from intracellular version
            try:
                m_c = model.metabolites.get_by_id(m.id.replace("_e", "_c"))
            except KeyError:
                logging.info("No intracellular version of {0} in the model".format(m.id))
            else:
                m.compartment = m.id[-1]
                m.name = m_c.name
                m.charge = m_c.charge
                m.formula = m_c.formula
                m.annotation = m_c.annotation

            # Add exchange reaction
            exchange_reaction = model.add_boundary(m)
            exchange_reaction.annotation["subsystem"] = "Exchange"
            exchange_reaction.annotation["SBO"]= "SBO:0000627"







if __name__ == '__main__':
    existing_reactions_fn = "../../ComplementaryData/curation/transport_reactions/updated_grRules.csv"
    new_transport_reactions_to_existing_metabolites_fn = "../../ComplementaryData/curation/transport_reactions/newTransportRxns.csv"
    new_transport_reactions_to_new_metabolites_fn = "../../ComplementaryData/curation/transport_reactions/NewTransportRxns_newMets.csv"
    new_metabolites_fn = "../../ComplementaryData/curation/transport_reactions/new_mets_annotation.csv"

    model = cobra.io.read_sbml_model("../../ModelFiles/xml/scoGEM.xml")
    # modify_existing_transporters(model, existing_reactions_fn)

    # new_transport_reactions_to_existing_metabolites(model, new_transport_reactions_to_existing_metabolites_fn)
    # new_transport_reactions_to_new_metabolites(model, new_transport_reactions_to_new_metabolites_fn, new_metabolites_fn)

    
    fix_transporters(model, existing_reactions_fn, new_transport_reactions_to_existing_metabolites_fn, 
                     new_transport_reactions_to_new_metabolites_fn, new_metabolites_fn)
    cobra.io.write_sbml_model(model, "../../test.xml")
    # cobra.io.save_yaml_model(model, "../../test.yml")
