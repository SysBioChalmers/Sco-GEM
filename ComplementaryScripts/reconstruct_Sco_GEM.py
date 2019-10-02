# -*- coding: utf-8 -*-
"""
This file reconstructs Sco_GEM, the genome-scale model for Streptomyces coelicolor A3(2).
Author: Snorre Sulheim
Created: 27.08.2018
Updated: 19.12.2018
Email: snorre.sulheim@sintef.no


# Description
The Sco_GEM community model of Streptomyces coelicolor is constructed using this script.
The reconstruction pipeline is based on the iKS1317 model. Below is a coarse list of the different steps of the pipeline

1.  Load iKS1317 and fix some issues
2a. Add reactions and metabolites from Sco4
2b. Update IDs and some annotations of Sco4
3.  Add reactions, metabolites, update biomass and gene-reaction-rules according to iAA1259
4.  Fix several small tasks / issues:
    - add SBO terms
    - add pseudometabolites
    - update metanetx annotations
    - update the biomass according to proteomics data
5.  Change reaction bounds according to dG values from eQuilibrator


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
from consensusModel import fix_transporters
from consensusModel import feat_annotations
from consensusModel import feat_subsystem_annotation
from consensusModel import issue_82_delete_reactions
from reversibility import reversibility
import export

REPO_DIR = Path(__file__).parent.parent.resolve()
COMPLEMENTARY_DATA = REPO_DIR / "ComplementaryData"
print("Repo dir: ", REPO_DIR)

SAVE_PATH = str(REPO_DIR/ "ModelFiles/xml/Sco-GEM.xml")
iKS1317_PATH = str(COMPLEMENTARY_DATA / "models/iKS1317.xml")


SCO4_PATH = str(COMPLEMENTARY_DATA / "models/Sco4.xml")
SCO4_REACTION_MAPPING_FN = str(COMPLEMENTARY_DATA / "curation/rxns_iKS1317_vs_Sco4.csv")
SCO4_METABOLITE_MAPPING_FN =  str(COMPLEMENTARY_DATA /"curation/mets_iKS1317_vs_Sco4.csv")
SCO4_REACTION_ANNOTATION_FN = str(COMPLEMENTARY_DATA /"curation/added_sco4_reactions.csv")
SCO4_METABOLITE_ANNOTATION_FN = str(COMPLEMENTARY_DATA /"curation/added_sco4_metabolites.csv")

iAA1259_PATH = str(COMPLEMENTARY_DATA / "models/iAA1259.xml")
iAA1259_NEW_REACTIONS_FN = str(COMPLEMENTARY_DATA / "curation/iAA1259_suppl_S4.csv") # New reactions

MET_TO_METANETX_FN = str(COMPLEMENTARY_DATA / "curation" /"metanetx_to_change.csv")
RXN_TO_METANETX_FN = str(COMPLEMENTARY_DATA / "curation" /"metanetx_reaction_annotations_to_change.csv")
MET_TO_CHEBI_FN = str(COMPLEMENTARY_DATA / "curation" /"chebi_annotation.csv")
NEW_BIOMASS_DATA_FN = str(COMPLEMENTARY_DATA / "biomass/biomass_scaled.txt")
EQUILIBRATOR_FN_1 = str(COMPLEMENTARY_DATA / "curation/reversibility/eQuilibrator_reversibility.csv")
EQUILIBRATOR_FN_2 = str(COMPLEMENTARY_DATA / "curation/reversibility/eQuilibrator_reversibility_lethals.csv")
ATP_DRIVEN_REACTIONS_REVERSIBILITY_FN = str(COMPLEMENTARY_DATA / "curation/reversibility/reversibility_ATP_driven_reactions.csv")

MODIFIED_TRANSPORT_REACTIONS_FN = str(COMPLEMENTARY_DATA / "curation" / "transport_reactions" / "updated_grRules.csv")
NEW_TRANSPORT_REACTIONS_FN = str(COMPLEMENTARY_DATA / "curation" / "transport_reactions" / "newTransportRxns.csv")
NEW_TRANSPORT_REACTIONS_TO_NEW_METABOLITES_FN = str(COMPLEMENTARY_DATA / "curation" / "transport_reactions" / "newTransportRxns_newMets.csv")
NEW_METABOLITES_TO_NEW_TRANSPORT_REACTIONS = str(COMPLEMENTARY_DATA / "curation" / "transport_reactions" / "new_mets_annotation.csv")

DOI_ANNOTATIONS_FN = str(COMPLEMENTARY_DATA / "annotations" / "reaction_notes_and_references.csv")
GENE_ANNOTATIONS_FN = str(COMPLEMENTARY_DATA / "annotations" / "genes.csv")
SUBSYSTEM_ANNOTATION_FN = str(COMPLEMENTARY_DATA / "curation" / "pathway_and_subsystem" / "subsystem_curation.csv")

# Settings
SOLVER = "glpk"

def reconstruct_Sco_GEM(model_fn, save_fn = None, write_requirements = True):
    Sco_GEM = cobra.io.read_sbml_model(model_fn)
    Sco_GEM.name = "Sco-GEM"
    Sco_GEM.id = "Sco-GEM"
    Sco_GEM.solver = SOLVER
    
    if save_fn is None:
        save_fn = model_fn

    # Part 1: Fix known issues in models
    ## 1a) Issues in iKS1317
    fix_iKS1317_issues.fix(Sco_GEM)

    ## 1b) Issues in Sco4 v4.00
    sco4_model = cobra.io.read_sbml_model(SCO4_PATH)
    sco4_model.solver = SOLVER
    fix_sco4_issues.fix(sco4_model)

    ## 1c) Add missing / changed gene annotations in iMK1208 identifed in Sco4 / and by Snorre 21.09.2018
    add_missing_gene_annotations_sco4.add_gene_annotations(Sco_GEM)

    # Part 2: Add reactions from Sco4
    Sco_GEM = add_reactions_from_sco4.add_reactions(sco4_model, Sco_GEM, SCO4_REACTION_MAPPING_FN, SCO4_METABOLITE_MAPPING_FN)

    ## 2b) Rename metabolites added from Sco4 to BIGGish Ids
    annotate_new_rxns_and_mets_from_sco4.add_rxn_annotations(Sco_GEM, SCO4_REACTION_ANNOTATION_FN, False)
    annotate_new_rxns_and_mets_from_sco4.add_met_annotations(Sco_GEM, SCO4_METABOLITE_ANNOTATION_FN, False)

    # Part 3: Add and modify reactions according to iAA1259
    iAA1259_model = cobra.io.read_sbml_model(iAA1259_PATH)
    iAA1259_model.solver = SOLVER
    add_and_modify_reactions_according_to_iAA1259.fix_iAA1259(iAA1259_model)
    Sco_GEM = add_and_modify_reactions_according_to_iAA1259.add_reactions(iAA1259_model, Sco_GEM, iAA1259_NEW_REACTIONS_FN)
    Sco_GEM = add_and_modify_reactions_according_to_iAA1259.modify_reactions(Sco_GEM)
    # Change biomass
    Sco_GEM = add_and_modify_reactions_according_to_iAA1259.change_biomass(iAA1259_model, Sco_GEM)

    # Part 4
    fix_issue12_reversibility.fix(Sco_GEM)
    fix_issue33_annotation_bugs.fix(Sco_GEM)
    redox_pseudometabolite.run(Sco_GEM)
    fix_SBO_terms.add_SBO(Sco_GEM)
    fix_issue33_annotation_bugs.fix_metanetx_metabolite_annotations(Sco_GEM, MET_TO_METANETX_FN)
    fix_biomass.fix_biomass(Sco_GEM, NEW_BIOMASS_DATA_FN)
    fix_issue33_annotation_bugs.apply_new_chebi_annotations(Sco_GEM, MET_TO_CHEBI_FN)
    fix_issue33_annotation_bugs.fix_c_c_in_metabolite_ids(Sco_GEM)
    fix_issue33_annotation_bugs.fix_metanetx_reaction_annotations(Sco_GEM, RXN_TO_METANETX_FN)


    # Part 5
    reversibility.change_bounds_according_to_eQuilibrator(Sco_GEM, EQUILIBRATOR_FN_1, EQUILIBRATOR_FN_2)
    reversibility.change_lower_bound_on_CPKS_reactions(Sco_GEM)
    reversibility.change_bounds_on_ATP_driven_reactions(Sco_GEM, ATP_DRIVEN_REACTIONS_REVERSIBILITY_FN)



    # Additional annotations 
    feat_annotations.add_doi_annotations(Sco_GEM, DOI_ANNOTATIONS_FN)
    feat_annotations.add_gene_annotations(Sco_GEM, GENE_ANNOTATIONS_FN)
    feat_subsystem_annotation.update_subsystem_annotations(Sco_GEM, SUBSYSTEM_ANNOTATION_FN)

    # Issue 82 Delete reactions without gene associations
    issue_82_delete_reactions.delete_reactions(Sco_GEM)

    # Issue 85 cpk exchange reaction
    add_and_modify_reactions_according_to_iAA1259.add_exchange_reaction_for_ycpk(Sco_GEM)
   
    #Part 6 - Add transport reactions
    fix_transporters.fix_transporters(Sco_GEM,MODIFIED_TRANSPORT_REACTIONS_FN, NEW_TRANSPORT_REACTIONS_FN,
                                      NEW_TRANSPORT_REACTIONS_TO_NEW_METABOLITES_FN, NEW_METABOLITES_TO_NEW_TRANSPORT_REACTIONS)

    # Save model
    export.export(Sco_GEM, formats = ["xml", "yml"], write_requirements = write_requirements)

if __name__ == '__main__':
    logging.basicConfig(filename='reconstruct_Sco_GEM.log', level=logging.INFO)
    reconstruct_Sco_GEM(iKS1317_PATH, SAVE_PATH)
    print("Finished reconstructing Sco-GEM")
