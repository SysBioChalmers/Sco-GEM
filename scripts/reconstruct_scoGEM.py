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

MODEL_PATH = "C:../model.xml"

def reconstruct_scoGEM(model_fn, save_fn = None):
    model = cobra.io.read_sbml_model(model_fn)
    if save_fn is None:
        save_fn = model_fn


    # Part 1: Fix known issues in iKS1317
    fix_iKS1317_issues.fix(model)


    # Save model
    ## Version number
    cobra.io.write_sbml_model(model, save_fn)


if __name__ == '__main__':
    reconstruct_scoGEM(MODEL_PATH)
