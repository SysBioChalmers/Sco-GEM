# -*- coding: utf-8 -*-
import cobra
import logging

REACTIONS_TO_JOIN = [  
 # NADH reaction    NADPH reaction       New reaction id    KEGG ID        New Name
 ('3OCHOCDH_NADH',   '3OCHOCDH_NADPH',     '3OCHOCDH',       None,         None), 
 ('3OXCOADH_NADH',   '3OXCOADH_NADPH',     '3OXCOADH',       None,         None),
 ('44DPCDH_NADH',    '44DPCDH_NADPH',      '44DPCDH',        "R11914",     None),
 ('4HYDPRO_NADH',    '4HYDPRO_NADPH',      '4HYDPRO',        "R11428",     None),
 ('4NITROB_NADH',    '4NITROB_NADPH',      '4NITROB',        None,         None),
 ('AHLGAL_NADH',     'AHLGAL_NADPH',       'AHLGAL',         None,         None), 
 ('AHOPS_NADH',      'AHOPS_NADPH',        'AHOPS',          None,         None),
 ('CADHX_NADH',      'CADHX_NADPH',        'CADHX',          None,         None),
 ('DDALLO_NADH',     'DDALLO_NADPH',       'DDALLO',         None,         None),
 ('DPCOX_NADH',      'DPCOX_NADPH',        'DPCOX',          "R09671",     None), 
 ('GDP64HRD_NADH',   'GDP64HRD_NADPH',     'GDP64HRD',       None,         None),
 ('CDP4D6DG_NADH',   '4HDAPMO',            '4HDAPMO',        "R06892",     "4-hydroxyacetophenone monooxygenase"),
 ('PHYFLUDS_NADH',   'PHYFLUDS_NADPH',     'PHYFLUDS',       "R04787",     None),
 ('PHYTDES_NADH',    'PHYTDES_NADPH',      'PHYTDES',        "R09692",     None),
 ('SORBDH_NADH',     'SORBDH_NADPH',       'SORBDH',         None,         None),
 ('ZCARDS_NADH',     'ZCARDS_NADPH',       'ZCARDS',         "R04798",     None),
 ('ZCAROTDH2_NADH',  'ZCAROTDH2_NADPH',    'ZCAROTDH2',      "R04800",     None)
 ]

REACTIONS_TO_KEEP = [
# Reaction ID     KEGG ID    New reaction ID
("HPYRR2_NADH",    "R01388", "HPYRRx"),
("HPYRR_NADPH",    "R01392", "HPYRRy"),
('2ABZMO_NADH',    "R03998", "2ABZMOx"),
('2ABZMO_NADPH',   "R03999", "2ABZMOy"),
('ANTDIO_NADH',    "R00823", "ANTDIOx"),
('ANTDIO_NADPH',   "R00825", "ANTDIOy"),
('CARNMO_NADH',    "R11875", "CARNMOx"),
('CARNMO_NADPH',   "R11911", "CARNMOy"),
('ORNHYDX_NADH',   "R10790", "ORNHYDXx"),
('ORNHYDX_NADPH',  "R10789", "ORNHYDXy"),
('CDP4DGNO',       "R03391", "CDP4DGNOx"),
('CDP4D6DG_NADPH', "R03392", "CDP4DGNOy"),
]      



# Notes

def run(scoGEM):
    # Create donor and acceptor pseudometabolites
    donor, acceptor = create_pseudometabolites()

    # Create reactions converting nad / nadp to acceptor and nadh/nadph to donor
    create_pseudoreactions(scoGEM, donor, acceptor)

    # Join the reactions given in REACTIONS_TO_JOIN
    join_reactions(scoGEM, donor, acceptor)

    # Change ID and add KEGG annotation of the reactions to keep
    keep_reactions(scoGEM)

def create_pseudoreactions(scoGEM, donor, acceptor):
    pseudo_donor_nadh = cobra.Reaction("PSEUDO_DONOR_NADH")
    pseudo_donor_nadh.name = "Pseudoreaction converting nadh to donor"
    pseudo_donor_nadh.bounds = (-1000, 1000)
    pseudo_donor_nadh.add_metabolites({scoGEM.metabolites.nadh_c: -1, donor: 1})


    pseudo_donor_nadph = cobra.Reaction("PSEUDO_DONOR_NADPH")
    pseudo_donor_nadph.name = "Pseudoreaction converting nadph to donor"
    pseudo_donor_nadph.bounds = (-1000, 1000)
    pseudo_donor_nadph.add_metabolites({scoGEM.metabolites.nadph_c: -1, donor: 1})
    
    pseudo_acceptor_nad = cobra.Reaction("PSEUDO_ACCEPTOR_NAD")
    pseudo_acceptor_nad.name = "Pseudoreaction converting nad to acceptor"
    pseudo_acceptor_nad.bounds = (-1000, 1000)
    pseudo_acceptor_nad.add_metabolites({scoGEM.metabolites.nad_c: -1, acceptor: 1})

    pseudo_acceptor_nadp = cobra.Reaction("PSEUDO_ACCEPTOR_NADP")
    pseudo_acceptor_nadp.name = "Pseudoreaction converting nadp to acceptor"
    pseudo_acceptor_nadp.bounds = (-1000, 1000)
    pseudo_acceptor_nadp.add_metabolites({scoGEM.metabolites.nadp_c: -1, acceptor: 1})
    scoGEM.add_reactions([pseudo_donor_nadh, pseudo_donor_nadph, pseudo_acceptor_nad, pseudo_acceptor_nadp])
    logging.info("Added 4 pseudoreactions converting nadp/nad and nadh/nadph to acceptor / donor")



def create_pseudometabolites():
    # Acceptor
    acceptor = cobra.Metabolite("acceptor_c") #NAD+
    acceptor.name = "Hydrogen-acceptor"
    acceptor.annotation["kegg.compound"] = "C00028"
    acceptor.annotation["chebi"] = "15339"
    acceptor.annotation["biocyc"] = "Acceptor"
    acceptor.formula = "R"
    acceptor.charge = 0
    acceptor.compartment = "c"
    
    
    # Donor
    donor = cobra.Metabolite("donor_c") #NADH
    donor.name = "Hydrogen-donor" 
    donor.annotation["kegg.compound"] = "C00030"
    donor.annotation["chebi"] = "17499"
    donor.annotation["biocyc"] = "Donor-H2"
    donor.formula = "RH"
    donor.charge = -1
    donor.compartment = "c"
    logging.info("Created two pseudometabolites: acceptor_c and donor_c")

    return donor, acceptor

def join_reactions(scoGEM, donor, acceptor):
    new_reactions = []
    remove_reactions = []
    for reaction_tuple in REACTIONS_TO_JOIN:
        r1 = scoGEM.reactions.get_by_id(reaction_tuple[0])
        r2 = scoGEM.reactions.get_by_id(reaction_tuple[1])
        new_id = reaction_tuple[2]
        kegg_id = reaction_tuple[3]
        new_name = reaction_tuple[4]

        new_reaction = cobra.Reaction(new_id)

        # NAME
        if new_name:
            new_reaction.name = new_name
        else:
            new_reaction.name = r1.name.replace("(NADH)", "")


        # Add metabolites
        mets_to_add = {}
        for m, coeff in r1.metabolites.items():
            if m.id in ["nadp_c", "nad_c"]:
                mets_to_add[acceptor] = coeff
            elif m.id in ["nadph_c", "nadh_c"]:
                mets_to_add[donor] = coeff
            else:
                mets_to_add[m] = coeff
        new_reaction.add_metabolites(mets_to_add)

        # Add genes
        if r1.gene_reaction_rule == r2.gene_reaction_rule:
            new_reaction.gene_reaction_rule = r2.gene_reaction_rule
        else:
            logging.warning("Different reaction rules: \n 1) {0} \n 2) {1}".format(r1.gene_reaction_rule, r2.gene_reaction_rule))


        # Add annotations
        # KEGG ID    
        if kegg_id:
            new_reaction.annotation["kegg.reaction"] = kegg_id
        new_reaction.annotation["bigg.reaction"] = new_id
        new_reaction.annotation["biocyc"] = r1.annotation["biocyc"].rsplit("_", 1)[0]

        for key in ["SBO", "ec-code", "metanetx.reaction", "origin", "rhea"]:
            try:
                r1_value = r1.annotation[key]
            except KeyError:
                r1_value = None

            try:
                r2_value = r2.annotation[key]
            except KeyError:
                r2_value = None
            
            if (r1_value is None) and (r2_value is None):
                continue
            elif (r1_value is None):
                new_reaction.annotation[key] = r2_value
            
            elif (r2_value is None):
                new_reaction.annotation[key] = r1_value
            else:
                if r1_value == r2_value:
                    new_reaction.annotation[key] = r1_value
                else:
                    new_reaction.annotation[key] = [r1_value, r2_value]
        # Remove r1 and r2
        new_reactions.append(new_reaction)
        remove_reactions += [r1, r2]
        logging.info("Removed reactions {0} and {1}, replaced by {2}".format(r1.id, r2.id, new_id))

    scoGEM.remove_reactions(remove_reactions)
    scoGEM.add_reactions(new_reactions)


def keep_reactions(scoGEM):
    for (r_id, kegg_id, new_id) in REACTIONS_TO_KEEP:
        r = scoGEM.reactions.get_by_id(r_id)
        r.annotation["kegg.reaction"] = kegg_id
        r.annotation["bigg.reaction"] = new_id
        logging.info("Changed id of reaction {0} to {1}".format(r.id, new_id))
        r.id = new_id


if __name__ == '__main__':
    scoGEM = cobra.io.read_sbml_model("../../ModelFiles/xml/scoGEM.xml")
    print(scoGEM.optimize())
    run(scoGEM)
    print(scoGEM.optimize())
    # cobra.io.write_sbml_model(scoGEM, "../../test.xml")