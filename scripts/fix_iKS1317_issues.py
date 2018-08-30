#!/usr/bin/env python
"""
This file fixes some of the errors in iKS1317
Author: Snorre Sulheim
Date: 27.08.2018

"""

import cobra

def fix(model):


    # Mass and charge balance
    model.reactions.NOR_syn.add_metabolites({model.metabolites.fdxox_c:-12, model.metabolites.fdxrd_c:12})

    

    # Ec code, kegg annotations and genes
    model.reactions.get_by_id("3OXCOAT").annotation["ec-code"] = "2.3.1.16" # Was 2.3.1.174
    model.reactions.get_by_id("3OXCOAT").gene_reaction_rule = "SCO1324 or SCO6027 or SCO6701 or SCO6967"
    model.reactions.MMSYNB.gene_reaction_rule = "s0001"

    # Fix id og mmy acetoacetyl and malonyl ACP
    model.metabolites.malACPmmy_c.charge = -1
    model.metabolites.actACPmmy_c.id = "temp"
    model.metabolites.malACPmmy_c.id = "actACPmmy_c"
    model.metabolites.actACPmmy_c.id = "malACPmmy_c"


    # Charge and reaction balancing
    model.reactions.OAADC.add_metabolites({model.metabolites.h_c: -1})
    model.reactions.SEPHCHCS.add_metabolites({model.metabolites.h_c: -1})
    model.reactions.DIOP5OR.add_metabolites({model.metabolites.h_c: 1})




    ## Phenylalanine metabolism
    model.reactions.PACCOAE.annotation["ec-code"] = "1.14.13.149"
    model.reactions.PACCOAE.annotation["kegg.reaction"] = "R09838"
    model.reactions.REPHACCOAI.annotation["ec-code"] = "5.3.3.18"
    model.reactions.REPHACCOAI.annotation["kegg.reaction"] = "R09837"
    model.metabolites.get_by_id("2oxpaccoa_c").annotation["kegg.compound"] = "C19975"
    model.metabolites.get_by_id("2oxpaccoa_c").annotation["chebi"] = "63252"
    model.metabolites.get_by_id("rephaccoa_c").annotation["kegg.compound"] = "C20062"



if __name__ == '__main__':
    iKS1317_PATH = "C:/Users/snorres/git/gem_sco/iKS1317.xml"
    model = cobra.io.read_sbml_model(iKS1317_PATH)
    fix(model)