#!/usr/bin/env python
"""
This file reconstructs scoGEM, the genome-scale model for Streptomyces coelicolor A3(2).
Author: Snorre Sulheim
Created: 27.08.2018
Modified: 30.08.2018


# Description
The scoGEM community model of Streptomyces coelicolor is constructed using the three

"""


import cobra
import logging

import fix_iKS1317_issues
import fix_sco4_issues
import add_reactions_from_sco4
import add_missing_gene_annotations_sco4
import add_and_modify_reactions_according_to_iAA1259

SAVE_PATH = "../../ModelFiles/xml/scoGEM.xml"
iKS1317_PATH = "../../ComplementaryData/models/iKS1317.xml"

SCO4_PATH = "../../ComplementaryData/models/Sco4.xml"
SCO4_REACTION_MAPPING_FN = "../../ComplementaryData/curation/rxns_iKS1317_vs_Sco4.csv"
SCO4_METABOLITE_MAPPING_FN =  "../../ComplementaryData/curation/mets_iKS1317_vs_Sco4.csv"

iAA1259_PATH = "../../ComplementaryData/models/iAA1259.xml"
iAA1259_NEW_REACTIONS_FN = "../../ComplementaryData/curation/iAA1259_suppl_S4.csv" # New reactions


def reconstruct_scoGEM(model_fn, save_fn = None):
    scoGEM = cobra.io.read_sbml_model(model_fn)
    scoGEM.name = "scoGEM"
    scoGEM.id = "scoGEM"
    
    if save_fn is None:
        save_fn = model_fn


    # Part 1: Fix known issues in models
    ## 1a) Issues in iKS1317
    fix_iKS1317_issues.fix(scoGEM)
    
    ## 1b) Issues in Sco4 v4.00
    sco4_model = cobra.io.read_sbml_model(SCO4_PATH)
    fix_sco4_issues.fix(sco4_model)

    ## 1c) Add missing / changeds gene annotations in iMK1208 identifed in Sco4 / and by Snorre 21.09.2018
    add_missing_gene_annotations_sco4.add_gene_annotations(scoGEM)

    # Part 2: Add reactions from Sco4
    scoGEM = add_reactions_from_sco4.add_reactions(sco4_model, scoGEM, SCO4_REACTION_MAPPING_FN, SCO4_METABOLITE_MAPPING_FN)


    # Part 3: Add and modify reactions according to iAA1259
    iAA1259_model = cobra.io.read_sbml_model(iAA1259_PATH)
    add_and_modify_reactions_according_to_iAA1259.fix_iAA1259(iAA1259_model)
    scoGEM = add_and_modify_reactions_according_to_iAA1259.add_reactions(iAA1259_model, scoGEM, iAA1259_NEW_REACTIONS_FN)
    scoGEM = add_and_modify_reactions_according_to_iAA1259.modify_reactions(scoGEM)
    # Change biomass
    scoGEM = add_and_modify_reactions_according_to_iAA1259.change_biomass(iAA1259_model, scoGEM)
    

    # Save model
    ## Version number
    cobra.io.write_sbml_model(scoGEM, save_fn)


if __name__ == '__main__':
    logging.basicConfig(filename='reconstruct_scoGEM.log', level=logging.INFO)
    reconstruct_scoGEM(iKS1317_PATH, SAVE_PATH)
