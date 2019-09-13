# -*- coding: utf-8 -*-
"""
This script removed reactions without any gene associations and which also create unrealistic loops 
when doing random sampling, this is described in issue #82

Author: Snorre Sulheim 
Created: 25.06.2019


"""

import pandas as pd
import cobra
from collections import defaultdict
import logging

LIST_OF_REACTIONS_TO_DELETE = ["ARAB2r",
                               "ABT2DG",
                               "ABTD",
                               "XYLK2",
                               "XYLR2",
                               "XYLR3",
                               "XYLTD_D",
                               "D5KGPA2",
                               "GLCOAS",
                               "LYXE",
                               "OCDOR",
                               "HPCOR",
                               "DM_ACT_c",
                               ]

def delete_reactions(model, reactions = LIST_OF_REACTIONS_TO_DELETE):
    model.remove_reactions(reactions, remove_orphans = True)

if __name__ == '__main__':
    model = cobra.io.read_sbml_model("../../ModelFiles/xml/scoGEM.xml")
    model.optimize()
    print(model.summary())
    delete_reactions(model)
    model.optimize()
    print(model.summary())
