#!/usr/bin/env python
"""
This file reconstructs scoGEM, the genome-scale model for Streptomyces coelicolor A3(2).
Author: Snorre Sulheim
Created: 27.08.2018


# Description
The scoGEM community model of Streptomyces coelicolor is constructed using the three

"""


import cobra
import logging

import fix_iKS1317_issues
import fix_sco4_issues
import add_missing_gene_annotations_sco4

SAVE_PATH = "../../ModelFiles/xml/scoGEM.xml"
iKS1317_PATH = "../../ComplementaryData/models/iKS1317.xml"
SCO4_PATH = "../../ComplementaryData/models/Sco4.xml"

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

    ## 1c) Add missing / changed gene annotations in iMK1208 identifed in Sco4 / and by Snorre 21.09.2018
    add_missing_gene_annotations_sco4.add_gene_annotations(scoGEM)

    # Save model
    ## Version number
    cobra.io.write_sbml_model(scoGEM, save_fn)


if __name__ == '__main__':
    logging.basicConfig(filename='reconstruct_scoGEM.log', level=logging.INFO)
    reconstruct_scoGEM(iKS1317_PATH, SAVE_PATH)
