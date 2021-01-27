"""

Author: Eduard Kerkhoven 
Created: 2021-01-19
Updated: 2021-01-24

Curations for release v1.3.0.

This script combines various curations on Sco-GEM v1.2.1 to produce v1.3.0.
Indicated is which issues are solved, more detailed explanation is given in the
relevant pull requests.

"""
import sys
sys.path.append("..")

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
    # Fixes Issue #88
    G6PDH2r = model.reactions.get_by_id('G6PDH2r')
    G6PDH2r.annotation['doi'] = '10.1371/journal.pone.0084151'
    G6PDH2r.notes['NOTES'] = 'See DOI, G6PDH can only use NADPH, not NADH.'
    model.reactions.remove('G6PDH1b')
    # As mentioned in #119
    XYLabc = model.reactions.get_by_id('XYLabc')
    XYLabc.gene_reaction_rule = '(SCO2404 or SCO3667) and SCO6010 and SCO6011 and SCO6009'
    # Fixes Issue #100
    CAT = model.reactions.get_by_id("CAT")
    CAT.gene_reaction_rule = "SCO0379 or SCO0560 or SCO0666 or SCO6204 or SCO7590"
    # Fixes Issue #89
    model.reactions.remove('G6PI')
    model.reactions.remove('PGIA')
    model.metabolites.remove('g6p_A_c')
    model.metabolites.remove('g6p_B_c')
    model.reactions.remove('BFBP')
    model.reactions.remove('TALAb')
    model.reactions.remove('TKT2h')
    model.metabolites.remove('f6p_B_c')
    model.reactions.remove('FBA_1')
    model.metabolites.remove('fdp_B_c')

    # Fix issue #127
    accoa_res = model.metabolites.get_by_id("accoa_res_c")
    accoa = model.metabolites.get_by_id("accoa_c")

    IPPS = model.reactions.get_by_id("IPPS")
    n = IPPS.metabolites.pop(accoa_res)
    IPPS.add_metabolites({accoa: n})

    MMSAD3 = model.reactions.get_by_id("MMSAD3")
    n = MMSAD3.metabolites.pop(accoa_res)
    MMSAD3.add_metabolites({accoa: n})


    

def add_gene_annotation(model):
    # Fixes #44
    df = pd.read_csv('../../data/curation/v130_uniprot_proteome_UP000001973.tab', sep = '\t')
    df.index=df['id'].str.replace('.','')
    df=df.fillna('')
    for g in model.genes:
        g.annotation['sbo'] = 'SBO:0000243'
        if g.id in df.index:
            pos = df.index.get_loc(g.id)
            if df.uniprot[pos]:
                g.annotation['uniprot'] = df.uniprot[pos]
            if df.pfam[pos]:
                tmp = df.pfam[pos].rstrip(';')
                g.annotation['pfam'] = tmp.split(';')
            if df.panther[pos]:
                tmp = df.panther[pos].rstrip(';')
                g.annotation['panther'] = tmp.split(';')   
            if df.go[pos]:
                tmp = df.go[pos].rstrip(';')
                tmp = tmp.replace(' ','')
                g.annotation['go'] = tmp.split(';')   
            if df.refseq[pos]:
                tmp = df.refseq[pos].rstrip(';')
                g.annotation['refseq'] = tmp.split(';')   
            if df.name[pos] != '':
                g.name = df.name[pos]
            else:
                g.name = df.id[pos] # If left empty, SBML I/O inserts "G_geneId"
        else:
            print("No annotation data for gene: {0}".format(g.id))

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
    # list_annotations(model)  # Only needs to be run once to gather metabolite IDs and pubchem.substance annotations
    correct_pubchem(model)
    add_gene_annotation(model)
    export.export(model, formats = ["xml", "yml"], write_requirements = 1, objective = "BIOMASS_SCO_tRNA")
