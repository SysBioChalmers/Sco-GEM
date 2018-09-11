#!/usr/bin/env python
"""
This file reconstructs scogem, the genome-scale model for Streptomyces coelicolor A3(2).
Author: Snorre Sulheim
Created: 27.08.2018
Modified: 30.08.2018


# Description
ToDo

"""


import cobra

import fix_iKS1317_issues
import fix_sco4_issues

MODEL_PATH = "../../ModelFiles/model.xml"
SCO4_PATH = "../../ComplementaryData/models/Sco4.xml"


def reconstruct_scoGEM(model_fn, save_fn = None):
    model = cobra.io.read_sbml_model(model_fn)
    if save_fn is None:
        save_fn = model_fn


    # Part 1: Fix known issues in models
    ## 1a) Issues in iKS1317
    fix_iKS1317_issues.fix(model)
    
    ## 1b) Issues in Sco4 v4.00
    Sco4 = fix_sco4_issues.fix(SCO4_PATH)

    # Part 2: Add reactions from Sco4
    

    # Save model
    ## Version number
    cobra.io.write_sbml_model(model, save_fn)


if __name__ == '__main__':
    reconstruct_scoGEM(MODEL_PATH)
