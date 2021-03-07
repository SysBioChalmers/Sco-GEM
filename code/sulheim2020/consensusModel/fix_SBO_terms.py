# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
Author: Snorre Sulheim / Eduard Kerkhoven

# Description
This file adds SBO terms to reactions and metabolites. This is kind of a
copy of a matlab-script Eduard used to to the same
"""

import cobra
import logging

"""

Reactions SBO terms
===================

SBO term    | name                  
------------------------------------
SBO:0000627 | exchange reaction                 
SBO:0000628 | demand reaction       
SBO:0000629 | biomass production    
SBO:0000630 | ATP maintenance      
SBO:0000395 | encapsulating process 
SBO:0000655 | transport reaction    
SBO:0000176 | biochemical reaction
SBO:0000632 | sink reaction
SBO:0000631 | pseudoreaction


Metabolite SBO terms
====================
SBO term    | name                  
------------------------------
SBO:0000649 | biomass
SBO:0000247 | simple chemical

Genes
======
All genes SBO:0000243


"""
METABOLITE_SBO_TERMS = {
    "biomass":         "SBO:0000649",
    "simple chemical": "SBO:0000247",

}

REACTION_SBO_TERMS = {
    "exchange reaction":     "SBO:0000627",
    "demand reaction":       "SBO:0000628",
    "biomass production":    "SBO:0000629",
    "ATP mainteinance":      "SBO:0000630",
    "encapsulating process": "SBO:0000395",
    "transport reaction":    "SBO:0000655",
    "biochemical reaction":  "SBO:0000176",
    "sink reaction":         "SBO:0000632",
    "pseudoreaction":        "SBO:0000631",
}
GENE_SBO_TERM = "SBO:0000243"




def add_SBO(scoGEM):
    # Metabolites
    for m in scoGEM.metabolites:
        if m.name in ["biomass", "DNA", "RNA", "protein", "carbohydrate", "cell wall", "lipid"]:
            m.annotation["SBO"] = METABOLITE_SBO_TERMS["biomass"]
        else:
            m.annotation["SBO"] = METABOLITE_SBO_TERMS["simple chemical"]

    # Reactions
    all_reactions = [r.id for r in scoGEM.reactions]
    for r in scoGEM.exchanges:
        r.annotation["SBO"] = REACTION_SBO_TERMS["exchange reaction"]
        all_reactions.remove(r.id)

    for r in scoGEM.demands:
        r.annotation["SBO"] = REACTION_SBO_TERMS["demand reaction"]
        all_reactions.remove(r.id)

    for r in scoGEM.sinks:
        r.annotation["SBO"] = REACTION_SBO_TERMS["sink reaction"]
        all_reactions.remove(r.id)

    scoGEM.reactions.ATPM.annotation["SBO"] = REACTION_SBO_TERMS["ATP mainteinance"]
    all_reactions.remove("ATPM")

    for r_id in all_reactions:
        r = scoGEM.reactions.get_by_id(r_id)
        if "BIOMASS_SCO" in r.id:
            r.annotation["SBO"] = REACTION_SBO_TERMS["biomass production"]
        elif r.id in ["PSEUDO_DONOR_NADH", "PSEUDO_DONOR_NADPH", "PSEUDO_ACCEPTOR_NAD", "PSEUDO_ACCEPTOR_NADP"]:
            r.annotation["SBO"] = REACTION_SBO_TERMS["pseudoreaction"]
        elif "pseudoreaction" in r.name.lower():
            r.annotation["SBO"] = REACTION_SBO_TERMS["encapsulating process"]
        else:

            if len(r.compartments) == 2:
                r.annotation["SBO"] = REACTION_SBO_TERMS["transport reaction"]
            else:
                r.annotation["SBO"] = REACTION_SBO_TERMS["biochemical reaction"]

    # Genes
    for g in scoGEM.genes:
        g.annotation["SBO"] = GENE_SBO_TERM

    logging.info("Added SBO terms to genes, reactions and metabolites")
    print("Added SBO terms")




if __name__ == '__main__':
    scoGEM_FN = "../../ModelFiles/xml/scoGEM.xml"
    scoGEM = cobra.io.read_sbml_model(scoGEM_FN)
    add_SBO(scoGEM)