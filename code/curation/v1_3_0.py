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
from pathlib import Path
import re

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
    
    # Fixes Issue #88, #106
    G6PDH2r = model.reactions.get_by_id('G6PDH2r')
    G6PDH2r.annotation['doi'] = '10.1371/journal.pone.0084151'
    G6PDH2r.notes['NOTES'] = 'See DOI, G6PDH can only use NADPH, not NADH.'
    model.reactions.remove('G6PDH1b')
    model.reactions.remove('GND2')

    # As mentioned in #119
    XYLabc = model.reactions.get_by_id('XYLabc')
    XYLabc.gene_reaction_rule = '(SCO2404 or SCO3667) and SCO6010 and SCO6011 and SCO6009'
    
    # Fixes Issue #100
    CAT = model.reactions.get_by_id("CAT")
    CAT.gene_reaction_rule = "SCO0379 or SCO0560 or SCO0666 or SCO6204 or SCO7590"
    # TODO: remove SCO2529 gene?
    
    # Fixes Issue #104
    model.reactions.remove('ABPYRATA')
    GABTA = model.reactions.get_by_id('GABTA')
    GABTA.annotation['kegg.reaction'] = 'R10178'

    # Fixes Issue #89
    ## Remove model reactions
    model.reactions.remove('G6PI')
    model.reactions.remove('PGIA')
    model.reactions.remove('GLUKB')
    model.reactions.remove('G6PBDH')
    model.reactions.remove('S6PGB')
    model.reactions.remove('BFBP')
    model.reactions.remove('TALAb')
    model.reactions.remove('TKT2h')
    model.reactions.remove('FBA_1')
    model.reactions.remove('FRUISOs')
    
    ## Replace specific isomers with generic versions
    g6p_B_c = model.metabolites.get_by_id('g6p_B_c')
    g6p_c = model.metabolites.get_by_id('g6p_c')
    TRE6PPP = model.reactions.get_by_id('TRE6PPP')
    TRE6PPP.add_metabolites({g6p_c: 1})
    TRE6PPP.subtract_metabolites({g6p_B_c: 1})
    
    g1p_B_c = model.metabolites.get_by_id('g1p_B_c')
    g1p_c = model.metabolites.get_by_id('g1p_c')
    TRE6PPP.add_metabolites({g1p_c: 1})
    TRE6PPP.subtract_metabolites({g1p_B_c: 1})

    MPL = model.reactions.get_by_id('MPL')
    MPL.add_metabolites({g1p_c: 1})
    MPL.subtract_metabolites({g1p_B_c: 1})

    TREPP = model.reactions.get_by_id('TREPP')
    TREPP.add_metabolites({g1p_c: 1})
    TREPP.subtract_metabolites({g1p_B_c: 1})

    ## Delete now orphaned metabolites
    model.metabolites.remove('g6p_A_c')
    model.metabolites.remove('g6p_B_c')
    model.metabolites.remove('f6p_B_c')  
    model.metabolites.remove('fdp_B_c')    
    model.metabolites.remove('fru_B_c')

	# Fixes Issue #105
    for r in model.reactions:
        r.notes.pop('GENE ASSOCIATION', None)
        r.notes.pop('GENE_ASSOCIATION', None)
        r.notes.pop('SUBSYSTEM', None)
        r.notes.pop('AUTHORS', None)
        r.notes.pop('GENE SYMBOL', None)
        r.notes.pop('GENE ASSOCIATION', None)
        r.annotation.pop('kegg.orthology', None)
        if 'Confidence Level' in list(r.notes.keys()):
            r.notes['CONFIDENCE LEVEL'] = r.notes.pop('Confidence Level')
        if 'biocyc' in list(r.annotation.keys()):
            oldValues = r.annotation['biocyc']
            if isinstance(oldValues, list):
                newValues = [re.sub(r'^(?!META:)','META:',value) for value in oldValues]
            else:
                newValues = re.sub(r'^(?!META:)',r'META:',oldValues)
            r.annotation['biocyc'] = newValues
        if 'doi' in list(r.annotation.keys()):
            if 'PMID' in r.annotation['doi']:
                pubmed = r.annotation.pop('doi')
                r.annotation['pubmed'] = re.sub('PMID: ','',pubmed)
    # Remove unused genes
    model.genes.remove('SCO3945')
    model.genes.remove('SCO3946')

def add_gene_annotation(model):
    # Fixes #44, #64
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
            else:
                g.annotation.pop('pfam', None)
            if df.panther[pos]:
                tmp = df.panther[pos].rstrip(';')
                g.annotation['panther'] = tmp.split(';')   
            if df.go[pos]:
                tmp = df.go[pos].rstrip(';')
                tmp = tmp.replace(' ','')
                g.annotation['go'] = tmp.split(';')
            else:
                g.annotation.pop('go', None)                
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

def correct_metabolite_annotations(model):
    # Contributes to Issues #105
    for m in model.metabolites:
        if 'biocyc' in list(m.annotation.keys()):
            oldValues = m.annotation['biocyc']
            if isinstance(oldValues, list):
                newValues = [re.sub(r'^(?!META:)','META:',value) for value in oldValues]
            else:
                newValues = re.sub(r'^(?!META:)',r'META:',oldValues)
            m.annotation['biocyc'] = newValues
        m.annotation.pop('obo.chebi',None)
        m.annotation.pop('3dmet',None)
        m.annotation.pop('cas',None)
        m.annotation.pop('kegg.drug',None)
        m.notes.pop('FORMULA',None)
        if 'metabetx.chemical' in list(m.annotation.keys()):
            m.annotation['metanetx.chemical'] = m.annotation.pop('metabetx.chemical')
    
    # Manually curate CHEBI that were delete above as obo.chebi annotation
    cellb_c = model.metabolites.get_by_id('cellb_c')
    cellb_c.annotation['chebi'] = ['CHEBI:36217', 'CHEBI:17057']
    cellb_e = model.metabolites.get_by_id('cellb_e')
    cellb_e.annotation['chebi'] = ['CHEBI:36217', 'CHEBI:17057']    
    pi_c = model.metabolites.get_by_id('pi_c')
    pi_c.annotation['chebi'] = ['CHEBI:35780', 'CHEBI:26078', 'CHEBI:39745', 'CHEBI:18367', 'CHEBI:43474', 'CHEBI:26020']

    # Fixes #62
    g3p_c = model.metabolites.get_by_id('g3p_c')
    g3p_c.annotation['kegg.compound'] = 'C00118'
    rmn_c = model.metabolites.get_by_id('rmn_c')
    rmn_c.annotation['kegg.compound'] = 'C02338'
    rmn_e = model.metabolites.get_by_id('rmn_e')
    rmn_e.annotation['kegg.compound'] = 'C02338'
    icit_c = model.metabolites.get_by_id('icit_c')
    icit_c.annotation['kegg.compound'] = 'C00451'
    orn_c = model.metabolites.get_by_id('orn_c')
    orn_c.annotation['kegg.compound'] = 'C00077'
    udcpp_c = model.metabolites.get_by_id('udcpp_c')
    udcpp_c.annotation['kegg.compound'] = 'C17556'  
    fe2_c = model.metabolites.get_by_id('fe2_c')
    fe2_c.annotation['kegg.compound'] = 'C14818'
    msh_c = model.metabolites.get_by_id('msh_c')
    msh_c.annotation['kegg.compound'] = 'C00051'
    mn2_c = model.metabolites.get_by_id('mn2_c')
    mn2_c.annotation['kegg.compound'] = 'C19610'
    mn2_e = model.metabolites.get_by_id('mn2_e')
    mn2_e.annotation['kegg.compound'] = 'C19610'
    ni2_e = model.metabolites.get_by_id('ni2_e')
    ni2_e.annotation['kegg.compound'] = 'C19609'
    ni2_c = model.metabolites.get_by_id('ni2_c')
    ni2_c.annotation['kegg.compound'] = 'C19609'    
    feroxB_c = model.metabolites.get_by_id('feroxB_c')
    feroxB_c.annotation.pop('kegg.compound', None)
    feroxB_e = model.metabolites.get_by_id('feroxB_e')
    feroxB_e.annotation.pop('kegg.compound', None)

def add_model_id(model):
    # Closes issue #128
    model.id = "Sco_GEM"

if __name__ == '__main__':
    model = export.get_latest_master_unversioned()
    misc_reaction_curations(model)
    # list_annotations(model)  # Only needs to be run once to gather metabolite IDs and pubchem.substance annotations
    correct_pubchem(model)
    add_gene_annotation(model)
    correct_metabolite_annotations(model)
    export.export(model, formats = "xml", write_requirements = 0, objective = "BIOMASS_SCO_tRNA")
    model = read_sbml_model("../../model/Sco-GEM.xml") # Extra round of I/O to consistently output single GO and PFAM annotations in YAML file
    export.export(model, formats = ["xml", "yml"], write_requirements = 1, objective = "BIOMASS_SCO_tRNA")
