#!/usr/bin/env python
"""
Author: Snorre Sulheim

# Description
This file fixes an issue related to reaction reversibility, see issue 12 for more info.

"""
import logging
import cobra

def fix(scoGEM):
    # Make BFBP irreversible
    scoGEM.reactions.BFBP.bounds = (0, 1000)

    # Change direction of ILETA, LEUTA, VALTA, NOR_syn and FNOR
    reaction_ids = ["ILETA", "LEUTA", "VALTA", "NOR_syn", "FNOR"]
    for r_id in reaction_ids:
        reaction = scoGEM.reactions.get_by_id(r_id)
        old_bound = reaction.bounds
        old_string = reaction.reaction
        
        reaction.bounds = (0, 1000)
        reverse_reaction_dict = {}
        for m, coeff in reaction.metabolites.items():
            reverse_reaction_dict[m] = -2*coeff
        reaction.add_metabolites(reverse_reaction_dict)

        logging.info("{0}: Changed bounds from {1} to {2}".format(r_id, old_bound, reaction.bounds))
        logging.info("{0}: Changed reaction direction, i.e. from {1} to {2}".format(r_id, old_string, reaction.reaction))
        


if __name__ == '__main__':
    scoGEM_FN = "../../ModelFiles/xml/scoGEM.xml"
    scoGEM = cobra.io.read_sbml_model(scoGEM_FN)
    fix(scoGEM)
