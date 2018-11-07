# -*- coding: utf-8 -*-
"""
This script updates the biomass equations according to modifications described in 
issue #13. The new biomass function is created using the coefficients given in
ComplementaryData/biomass/biomass_scaled.txt

Author: Snorre Sulheim / Eduard Kerkhoven
Created: 05.11.2018

The biomass reaction is seperated into the follwoing pseudoreactions:
    - R: RNA_PSEUDO
    - B: BIOMASS_SCO
    - L: LIPID_PSEUDO
    - P: PROTEIN_PSEUDO
    - W: CELL_WALL_PSEUDO
    - C: CARBOHYDRATE_PSEUDO
    - M: MISC_PSEUDO
    - D: DNA_PSEUDO


"""
import pandas as pd
import cobra
from collections import defaultdict
import logging

from consensusModel.fix_SBO_terms import REACTION_SBO_TERMS, METABOLITE_SBO_TERMS


BIOMASS_PSEUDOREACTIONS = {
    "R": ("RNA_PSEUDO", "RNA pseudoreaction"),
    "B": ("BIOMASS_SCO", "S. coelicolor biomass objective function - with 75.79 GAM estimate"),
    "L": ("LIPID_PSEUDO", "Lipid pseudoreaction"),
    "P": ("PROTEIN_PSEUDO", "Protein pseudoreaction"),
    "W": ("CELL_WALL_PSEUDO", "Cell wall pseudoreaction"),
    "C": ("CARBOHYDRATE_PSEUDO", "Carbohydrate pseudoreaction"),
    "M": ("MISC_PSEUDO", "Misc pseudoreaction"),
    "D": ("DNA_PSEUDO", "DNA pseudoreaction")
}

# PSEUDOMETABOLITES_ANNOTATION = {
#     "lipid_c": ["chebi:18509", "kegg.compound:C01356"],
#     "protein_c": ["chebi:36080"]

# }
AMINO_ACID_DICT = {"glutamate": ("glu__L_c", "glutrna_c", "trnaglu_c"), 
                   "alanine": ("ala__L_c", "alatrna_c", "trnaala_c"), 
                   "arginine": ("arg__L_c", "argtrna_c", "trnaarg_c"),
                   "asparagine": ("asn__L_c", "asntrna_c", "trnaasn_c"),
                   "aspartate": ("asp__L_c", "asptrna_c", "trnaasp_c"),
                   "cysteine": ("cys__L_c", "cystrna_c", "trnacys_c"),
                   "Methionine": ("met__L_c", "mettrna_c", "trnamet_c"),
                    
                   "Glutamine": ("gln__L_c", "glntrna_c", "trnagln_c"),
                   "Glycine": ("gly_c", "glytrna_c", "trnagly_c"),
                   "Histidine": ("his__L_c", "histrna_c", "trnahis_c"),
                   "Isoleucine": ("ile__L_c", "iletrna_c", "trnaile_c"),
                   "Leucine": ("leu__L_c", "leutrna_c", "trnaleu_c"),
                   "Lysine": ("lys__L_c", "lystrna_c", "trnalys_c"),
                   "phenylalanine": ("phe__L_c", "phetrna_c", "trnaphe_c"),
                   "Proline": ("pro__L_c", "protrna_c", "trnapro_c"),
                   "Serine": ("ser__L_c", "sertrna_c", "trnaser_c"),
                   "Threonine": ("thr__L_c", "thrtrna_c", "trnathr_c"),
                   "Tryptophan": ("trp__L_c", "trptrna_c", "trnatrp_c"),
                   "Tyrosine": ("tyr__L_c", "tyrtrna_c", "trnatyr_c"),
                   "Valine": ("val__L_c", "valtrna_c", "trnaval_c"),
                    
                    # ""fmettrna_c", # N-formylmethionine
}

BIOMASS_DATA_FN = "../../ComplementaryData/biomass/biomass_scaled.txt"

def read_biomass_scaled(biomass_data_fn):
    df = pd.read_csv(biomass_data_fn, sep = "\t", header = 0)
    return df

def fix_biomass(model, biomass_data_fn, new_reactions_dict = BIOMASS_PSEUDOREACTIONS):
    # Delete existing biomass function
    model.reactions.get_by_id("BIOMASS_SCO").remove_from_model()
    model.reactions.get_by_id("BIOMASS_SCO_tRNA").remove_from_model()
        
    biomass_scaled_df = read_biomass_scaled(biomass_data_fn)
    
    # Add new biomass reactions
    create_pseudo_reactions(model, new_reactions_dict)

    reaction_metabolite_dict = defaultdict(dict)
    # Add metabolites according to biomass dataframe
    for index, row in biomass_scaled_df.iterrows():
        try:
            m = model.metabolites.get_by_id(row["met"])
        except KeyError:
            m = create_pseudometabolite(row)
            model.add_metabolites([m])
        
        reaction_metabolite_dict[row["rxn"]][m] = row["S"]

    for key, met_dict in reaction_metabolite_dict.items():
        r_id = new_reactions_dict[key][0]
        reaction = model.reactions.get_by_id(r_id)
        reaction.add_metabolites(met_dict)

    # Set new biomass reaction as objective
    model.reactions.get_by_id("BIOMASS_SCO").objective_coefficient = 1


    add_biomass_tRNA(model)


def create_pseudometabolite(row):
    print("Creating metabolite with id: ", row["met"])
    m = cobra.Metabolite(row["met"])
    m.name = row["name"]
    m.annotation["SBO"] = METABOLITE_SBO_TERMS["biomass"]
    m.compartment = "c"
    m.formula = "R"
    m.charge = 0
    return m




def create_pseudo_reactions(model, new_reactions_dict):
    for key, tupl in new_reactions_dict.items():
        r_id, r_name = tupl
        r = cobra.Reaction(r_id)
        r.name = r_name
        r.annotation["SBO"] = REACTION_SBO_TERMS["biomass production"]
        r.annotation["subsystem"] = 'Biomass and maintenance functions'
        r.lower_bound = 0
        r.upper_bound = 1000
        model.add_reaction(r)


def add_biomass_tRNA(model):
    # Make biomass _reaction
    biomass_trna = model.reactions.get_by_id("BIOMASS_SCO").copy()
    biomass_trna.id = "BIOMASS_SCO_tRNA"
    biomass_trna.name = "S. coelicolor biomass function with AA replaced by tRNA-AA"
    model.add_reaction(biomass_trna)
    # tRNA protein reaction and pseudometabolite
    m_p = model.metabolites.protein_c.copy()
    m_p.id = "protein_tRNA_c"
    m_p.name = "Protein pseudometabolite created with tRNA-AA instead of AA"

    r_p = model.reactions.PROTEIN_PSEUDO.copy()
    r_p.id = "PROTEIN_PSEUDO_tRNA"
    r_p.name = "Protein pseudoreaction with tRNA-AA instead of AA"
    model.add_reaction(r_p)

    # Change coefficients
    summed_coefficient = 0

    for key, tup in AMINO_ACID_DICT.items():
        base_met = model.metabolites.get_by_id(tup[0])
        trna_complex_met = model.metabolites.get_by_id(tup[1])
        trna_premodule_met = model.metabolites.get_by_id(tup[2])

        # Get Biomass coefficient for that metabolite
        biomass_coeff = r_p.get_coefficient(base_met.id)
        summed_coefficient += biomass_coeff

        # Add trnA metabolites and remove aminoacid
        r_p.add_metabolites({base_met: -biomass_coeff, trna_complex_met: biomass_coeff, trna_premodule_met: -biomass_coeff})

    # Remove excessive ATP, ADP, proton and PI and h2o from biomass with tRNA
    m_atp = model.metabolites.atp_c
    m_adp = model.metabolites.adp_c
    m_pi = model.metabolites.pi_c
    m_h = model.metabolites.h_c
    m_h2o = model.metabolites.h2o_c
    biomass_trna.add_metabolites({m_atp: -summed_coefficient, m_h2o: -summed_coefficient, m_adp: summed_coefficient, 
                                m_pi: summed_coefficient, m_h: summed_coefficient})

    # Replace protein metabolite in biomass
    biomass_trna.add_metabolites({model.metabolites.protein_c: 1, m_p:-1})

    # Replace protein metabolite un protein pseudoreaction with tRNA
    r_p.add_metabolites({model.metabolites.protein_c: -1, m_p:1})

    # model.add_reactions([biomass_trna, r_p])
    # model.reactions.BIOMASS_SCO_tRNA.objective_coefficient = 0





if __name__ == '__main__':
    model = cobra.io.read_sbml_model("../../ModelFiles/xml/scoGEM.xml")
    model.optimize()
    print(model.summary())

    fix_biomass(model, BIOMASS_DATA_FN)

    for key, tup in BIOMASS_PSEUDOREACTIONS.items():
        print(model.reactions.get_by_id(tup[0]))
    # cobra.io.write_sbml_model(model, "../../test.xml")

    
    print(model.optimize())

