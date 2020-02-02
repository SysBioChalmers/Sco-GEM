# -*- coding: utf-8 -*-
"""
This file fixes some of the known errors in Sco4 4.0.0
Author: Snorre Sulheim
Date: 10.09.2018

The changes were suggested by Eduard Kerkhoven, and are the following.
# Deleted as they are duplicates: 
- RXN0-5224
- METHYLGLUTACONYL-COA-HYDRATASE-RXN
- GLU6PDEHYDROG-RXN
- RXN-15856
- 1.14.13.84-RXN_NADPH
- R03998
- R03999
- R09692_NADPH
- RXN-9930
- 1.17.1.1-RXN_NADH
- R09692_NADH

# Gene annotations:
- MGCH should have gene SCO4930

"""
import cobra

SCO4_PATH = "../../ComplementaryData/models/Sco4.xml"

DELETE_REACTION_LIST = ["RXN0-5224",
                        "METHYLGLUTACONYL-COA-HYDRATASE-RXN",
                        "GLU6PDEHYDROG-RXN",
                        "RXN-15856",
                        "1.14.13.84-RXN_NADPH",
                        "R03998",
                        "R03999",
                        "R09692_NADPH",
                        "RXN-9930",
                        "1.17.1.1-RXN_NADH",
                        "R09692_NADH"]


def fix(sco4_model, model_fn = None, save = False):
    """
    Deleted reactions from Sco4 that were identified as wrong during Sco4 model development

    Note: This doesn't work with cobra > 0.16 because they then handle ascii conversion upon reading the model
    """
    converted_ids = convert_to_ascii_codes(DELETE_REACTION_LIST)
    for reaction_id in converted_ids:
        sco4_model.reactions.get_by_id(reaction_id).remove_from_model()

    if save:
        if not model_fn:
            model_fn = "sco4.xml"
        cobra.io.write_sbml_model(sco4_model, model_fn)
        
def convert_to_ascii_codes(reaction_list):
    new_ids = []
    for r_id in reaction_list:
        new_ids.append(r_id.replace("-", "__45__").replace(".", "__46__"))
    return new_ids


if __name__ == '__main__':
    model = cobra.io.read_sbml_model(SCO4_PATH)
    model = fix(model)