# -*- coding: utf-8 -*-
"""
This file reconstructs scoGEM, the genome-scale model for Streptomyces coelicolor A3(2).
Author: Snorre Sulheim
Created: 27.08.2018
Updated: 12.12.2018


# Description
The scoGEM community model of Streptomyces coelicolor is constructed using this script

"""


import cobra
import logging
from pathlib import Path

from consensusModel import fix_iKS1317_issues
from consensusModel import fix_sco4_issues
from consensusModel import add_missing_gene_annotations_sco4
from consensusModel import add_reactions_from_sco4
from consensusModel import annotate_new_rxns_and_mets_from_sco4
from consensusModel import add_and_modify_reactions_according_to_iAA1259
from consensusModel import fix_issue12_reversibility
from consensusModel import fix_issue33_annotation_bugs
from consensusModel import redox_pseudometabolite
from consensusModel import fix_SBO_terms
from consensusModel import fix_biomass
from consensusModel import feat_annotations
from reversibility import reversibility
import export


SAVE_PATH = "../ModelFiles/xml/scoGEM.xml"
iKS1317_PATH = "../ComplementaryData/models/iKS1317.xml"

REPO_DIR = Path(__file__).parent.parent

SCO4_PATH = str(REPO_DIR / "ComplementaryData/models/Sco4.xml")
SCO4_REACTION_MAPPING_FN = "../ComplementaryData/curation/rxns_iKS1317_vs_Sco4.csv"
SCO4_METABOLITE_MAPPING_FN =  "../ComplementaryData/curation/mets_iKS1317_vs_Sco4.csv"
SCO4_REACTION_ANNOTATION_FN = "../ComplementaryData/curation/added_sco4_reactions.csv"
SCO4_METABOLITE_ANNOTATION_FN = "../ComplementaryData/curation/added_sco4_metabolites.csv"

iAA1259_PATH = "../ComplementaryData/models/iAA1259.xml"
iAA1259_NEW_REACTIONS_FN = "../ComplementaryData/curation/iAA1259_suppl_S4.csv" # New reactions

MET_TO_METANETX_FN = str(REPO_DIR / "ComplementaryData" / "curation" /"metanetx_to_change.csv")
RXN_TO_METANETX_FN = str(REPO_DIR / "ComplementaryData" / "curation" /"metanetx_reaction_annotations_to_change.csv")
MET_TO_CHEBI_FN = str(REPO_DIR / "ComplementaryData" / "curation" /"chebi_annotation.csv")
NEW_BIOMASS_DATA_FN = "../ComplementaryData/biomass/biomass_scaled.txt"
EQUILIBRATOR_FN_1 = "../ComplementaryData/curation/reversibility/eQuilibrator_reversibility.csv"
EQUILIBRATOR_FN_2 = "../ComplementaryData/curation/reversibility/eQuilibrator_reversibility_lethals.csv"

DOI_ANNOTATIONS_FN = str(REPO_DIR / "ComplementaryData" / "annotations" / "reaction_notes_and_references.csv")

def reconstruct_scoGEM(model_fn, save_fn = None, write_requirements = True):
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

    ## 1c) Add missing / changed gene annotations in iMK1208 identifed in Sco4 / and by Snorre 21.09.2018
    add_missing_gene_annotations_sco4.add_gene_annotations(scoGEM)

    # Part 2: Add reactions from Sco4
    scoGEM = add_reactions_from_sco4.add_reactions(sco4_model, scoGEM, SCO4_REACTION_MAPPING_FN, SCO4_METABOLITE_MAPPING_FN)

    ## 2b) Rename metabolites added from Sco4 to BIGGish Ids
    annotate_new_rxns_and_mets_from_sco4.add_rxn_annotations(scoGEM, SCO4_REACTION_ANNOTATION_FN, False)
    annotate_new_rxns_and_mets_from_sco4.add_met_annotations(scoGEM, SCO4_METABOLITE_ANNOTATION_FN, False)

    # Part 3: Add and modify reactions according to iAA1259
    iAA1259_model = cobra.io.read_sbml_model(iAA1259_PATH)
    add_and_modify_reactions_according_to_iAA1259.fix_iAA1259(iAA1259_model)
    scoGEM = add_and_modify_reactions_according_to_iAA1259.add_reactions(iAA1259_model, scoGEM, iAA1259_NEW_REACTIONS_FN)
    scoGEM = add_and_modify_reactions_according_to_iAA1259.modify_reactions(scoGEM)
    # Change biomass
    scoGEM = add_and_modify_reactions_according_to_iAA1259.change_biomass(iAA1259_model, scoGEM)

    # Part 4
    fix_issue12_reversibility.fix(scoGEM)
    fix_issue33_annotation_bugs.fix(scoGEM)
    redox_pseudometabolite.run(scoGEM)
    fix_SBO_terms.add_SBO(scoGEM)
    fix_issue33_annotation_bugs.fix_metanetx_metabolite_annotations(scoGEM, MET_TO_METANETX_FN)
    fix_biomass.fix_biomass(scoGEM, NEW_BIOMASS_DATA_FN)
    fix_issue33_annotation_bugs.apply_new_chebi_annotations(scoGEM, MET_TO_CHEBI_FN)
    fix_issue33_annotation_bugs.fix_c_c_in_metabolite_ids(scoGEM)
    fix_issue33_annotation_bugs.fix_metanetx_reaction_annotations(scoGEM, RXN_TO_METANETX_FN)

    # Additional annotations 
    feat_annotations.add_doi_annotations(scoGEM, DOI_ANNOTATIONS_FN)


    # Part 5
    scoGEM = reversibility.change_bounds_according_to_eQuilibrator(scoGEM, EQUILIBRATOR_FN_1, EQUILIBRATOR_FN_2)

    # Save model
    export.export(scoGEM, formats = ["xml", "yml"], write_requirements = write_requirements)

if __name__ == '__main__':
    logging.basicConfig(filename='reconstruct_scoGEM.log', level=logging.INFO)
    reconstruct_scoGEM(iKS1317_PATH, SAVE_PATH)
    print("Finished reconstructing scoGEM")
