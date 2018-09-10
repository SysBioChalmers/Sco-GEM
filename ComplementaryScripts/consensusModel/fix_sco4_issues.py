# -*- coding: <encoding name> -*-
"""
This file fixes some of the known errors in Sco4 4.0.0
Author: Snorre Sulheim
Date: 27.08.2018

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

SCO4_FN = "../../ComplementaryData/models/Sco4.xml"

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


def fix_sco4_issues(model_fn = SCO4_FN):
    model = cobra.io.read_sbml_model(model_fn)

    for reaction_id in DELETE_REACTION_LIST:
        model.reactions.get_by_id(reaction_id).remove_from_model()
    return model


if __name__ == '__main__':
    model = fix_sco4_issues(SCO4_FN)
    cobra.io.write_sbml_model(model, SCO4_FN)