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

    model.reactions.PGMT.annotation["ec-code"] = "5.4.2.2"
    model.reactions.PPM.annotation["ec-code"] = "5.4.2.2"
    model.reactions.ME1.annotation["ec-code"] = "1.1.1.38"

    model.reactions.GLUDyi.annotation["ec-code"] = "1.4.1.4"
    model.reactions.GLUDyi.gene_reaction_rule = "SCO4683"

    model.reactions.GLUSx.gene_reaction_rule = "(SCO1977 or SCO2025) and SCO2026"
    model.reactions.GLUSy.gene_reaction_rule = "(SCO1977 or SCO2025) and SCO2026"
    model.reactions.GLUN.annotation["ec-code"] = "3.5.1.2"


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


    # Genes without reaction annotated
    model.reactions.ALLTAMH2.gene_reaction_rule = "SCO3072"
    model.genes.remove("SCO5184")
    model.genes.remove("SCO6640")
    model.genes.remove("SCO5188")
    
    # SCO6297 - Add reaction R02395
    add_R02395(model)

    # ADD R02301
    add_R02301(model)

    # Metabolites without reactions
    model.metabolites.get_by_id("3hmp_c").remove_from_model()
    model.metabolites.get_by_id("malylcoa_c").remove_from_model()
    model.metabolites.get_by_id("c78dhguantp_c").remove_from_model()

def add_R02301(model):
    reaction = cobra.Reaction("FTHFCL_1")
    reaction.name = "5-Formyltetrahydrofolate cyclo-ligase"
    metabolite_dict = {model.metabolites.atp_c: -1,
                       model.metabolites.get_by_id("5fthf_c"): -1,
                       model.metabolites.adp_c: 1,
                       model.metabolites.pi_c: 1,
                       model.metabolites.methf_c: 1}

    reaction.add_metabolites(metabolite_dict)
    reaction.annotation["ec-code"] = "6.3.3.2"
    reaction.annotation["kegg.reaction"] = "R02301"
    reaction.annotation["origin"] = "KEGG"
    reaction.gene_reaction_rule = "SCO3183"
    reaction.bounds = (0, 1000)
    model.add_reaction(reaction)
    # return model


def add_R02395(model):
    carnitine_c = cobra.Metabolite("carnitine_c", "C7H16NO3", charge = 0, compartment = "c")
    carnitine_c.name = "Carnitine"
    carnitine_c.annotation["kegg.compound"] = "C00487"
    carnitine_c.annotation["origin"] = "KEGG"
    
    dhdcarn_c = cobra.Metabolite("dhdcarn_c", "C7H14NO3", charge = 0, compartment = "c")
    dhdcarn_c.name = "3-Dehydrocarnitine"
    dhdcarn_c.annotation["kegg.compound"] = "C02636"
    dhdcarn_c.annotation["origin"] = "KEGG"

    model.add_metabolites([carnitine_c, dhdcarn_c])

    reaction = cobra.Reaction("CARNOX")
    reaction.name = "Carnitine:NAD+ 3-oxidoreductase"
    metabolite_dict = {model.metabolites.get_by_id("carnitine_c"): -1,
                       model.metabolites.get_by_id("dhdcarn_c"): 1,
                       model.metabolites.nad_c: -1,
                       model.metabolites.nadh_c: 1,
                       model.metabolites.h_c: 1}
    reaction.add_metabolites(metabolite_dict)
    reaction.annotation["ec-code"] = "1.1.1.108"
    reaction.annotation["kegg.reaction"] = "R02395"
    reaction.annotation["origin"] = "KEGG"
    reaction.gene_reaction_rule = "SCO6297"
    reaction.bounds = (-1000, 1000)


    model.add_reaction(reaction)

if __name__ == '__main__':
    scoGEM_PATH = "../../ModelFiles/scoGEM.xml"
    model = cobra.io.read_sbml_model(scoGEM_PATH)
    add_R02395(model)
    print(model.optimize)