"""

Author: Eduard Kerkhoven 
Created: 2021-01-19
Updated: 2021-01-20

Curations for release v1.3.0.

This script combines various curations on Sco-GEM v1.2.1 to produce v1.3.0.
Indicated is which issues are solved, more detailed explanation is given in the
relevant pull requests.

"""
from cobra.io import read_sbml_model, write_sbml_model
import export
import pandas as pd

def misc_reaction_curations(model):
    # Fixes Issue #87
    ICDHyr = model.reactions.get_by_id("ICDHyr")
    ICDHyr.lower_bound = 0
    # Fixes Issue #90
    FNOR = model.reactions.get_by_id("FNOR")
    FNOR.lower_bound = 0
    FNOR.annotation['doi'] = '10.1016/j.abb.2008.02.014'
    FNOR.notes['NOTES'] = 'See DOI, in non-photosynthetic bacteria this reaction is only in forward direction.'
    # Fixes Issue #33
    ZCAROTDH2 = model.reactions.get_by_id("ZCAROTDH2")
    ZCAROTDH2.annotation['biocyc'] = 'META:RXN-12412'
    # Fixes Issue #111
    SUCD9 = model.reactions.get_by_id("SUCD9")
    SUCD9.annotation['ec-code'] = '1.3.5.1'
    PFK_2 = model.reactions.get_by_id("PFK_2")
    PFK_2.annotation['ec-code'] = '2.7.1.144'
    APAT2r = model.reactions.get_by_id("APAT2r")
    APAT2r.annotation['ec-code'] = '2.6.1.-'
    HADPCOADH3 = model.reactions.get_by_id("HADPCOADH3")
    HADPCOADH3.annotation['ec-code'] = '1.1.1.35'  
    _45DOPA = model.reactions.get_by_id("45DOPA")
    _45DOPA.annotation['kegg.reaction'] = 'R08836'  

def remove_annotation_type(model):
    # Contributes to Issue #105
    for g in model.genes:
        g.annotation.pop('go', None)
        g.annotation.pop('pfam', None)

def list_annotations(model):
    # Contributes to Issue #105
    df = {'metabolite': [], 'pubchem_substance': []}
    for m in model.metabolites:
        if 'pubchem.substance' in m.annotation.keys():
            df['reaction'].append(m.id)
            df['pubchem_substance'].append(m.annotation['pubchem.substance'])
    df = pd.DataFrame(df)
    df.to_csv('../../data/curation/v130_pubchem_substance.csv', sep = ';', index = 0)

def correct_pubchem(model):
    # Contributes to Issue #105
    df = pd.read_csv('../../data/curation/v130_pubchem_substance.csv', sep = ';')
    m = model.metabolites.get_by_id('3oddecACP_c') # Has no replacement pubchem.compound, remove manually
    m.annotation.pop('pubchem.substance', None)
    for idx,id in enumerate(df.metabolite):
        m = model.metabolites.get_by_id(id)
        m.annotation.pop('pubchem.substance', None)
        m.annotation['pubchem.compound'] = str(df.pubchem_compound[idx])

if __name__ == '__main__':
    model = read_sbml_model("../../model/xml/Sco-GEM.xml")
    misc_reaction_curations(model)
    remove_annotation_type(model)  
    # list_annotations(model)  # Only needs to be run once to gather metabolite IDs and pubchem.substance annotations
    correct_pubchem(model)
    export.export(model, formats = ["xml", "yml"], write_requirements = 1, objective = "BIOMASS_SCO_tRNA")
