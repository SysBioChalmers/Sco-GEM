"""
Curations for release v1.3.2.

This script combines various curations on Sco-GEM v1.3.1 to produce v1.3.2.
Indicated is which issues are solved, more detailed explanation is given in the
relevant pull requests.

"""
# Import required functions, add any that are necessary for your code.
import sys
import cobra
from dotenv import find_dotenv
import pandas as pd

# Find .env + define paths. REPO_PATH refers to the repository root folder, and
# can be used to conveniently define absolute paths.
REPO_PATH = find_dotenv()
REPO_PATH = REPO_PATH[:-5]

sys.path.append(REPO_PATH + "/code")
import export

def fix_stoichiometric_consistency(model):
    # Add formula for metabolites without this information
    # Based on a smart table generated from biocyc
    _add_biocyc_formulas(model)


    # Add manual curations
    fn = "../../data/curation/v1_3_2/v1_3_2_curated_metabolite_charge_formula.csv"
    df = pd.read_csv(fn, header = 0)
    for _, row in df.iterrows():
        m = model.metabolites.get_by_id(row["ID"])
        m.charge = row["Charge"]
        m.formula = row["Formula"]

    # Fix reactions containing donor / acceptor pseduo-metabolites
    _add_proton_to_donor_reactions(model)

    # Fix reactions related to Electron-transferring flavoproteins (ETF)
    _remove_proton_from_ETF_reactions(model)

    _edit_reactions(model)

def _add_biocyc_formulas(model):
    fn =  "../../data/curation/v1_3_2/v1_3_2_biocyc_mets_missing_formula.tsv"
    df = pd.read_csv(fn, header = 0, sep = "\t", index_col = 0)
    for met_id, row in df.iterrows():
        m = model.metabolites.get_by_id(met_id)
        m.formula = row['Chemical Formula']

def _add_proton_to_donor_reactions(model):
    """
    Many of the stoichiometric inconcistenscies derive from the donor/acceptor reactions. These donor / acceptor metabolites represent nad(p)/nad(p)h. These reactions often also feature a h+, but this seems to be omitted in these reactions, leading to a charge and mass balance of 
    """
    reaction_ids = ['3OCHOCDH', '3OXCOADH', '44DPCDH', '4HYDPRO', 'AHOPS',
                     'CADHX', 'DPCOX', 'PHYFLUDS', 'PHYTDES', 'SORBDH', 'ZCARDS', 'ZCAROTDH2']
    
    # Get metabolites 
    m_donor = model.metabolites.get_by_id('donor_c')
    m_h = model.metabolites.get_by_id('h_c')

    for r_id in reaction_ids:
        r = model.reactions.get_by_id(r_id)
        s_m_donor = r.metabolites[m_donor]
        r.add_metabolites({m_h:s_m_donor})
        if len(r.check_mass_balance()):
            print("Error!", r, r.check_mass_balance())

def _remove_proton_from_ETF_reactions(model):
    """
    According to biocyc these reactions consumes / produces a proton, but that doesn't match the charge and mass balance
    """
    reaction_ids = ['BCOADH', 'BSUZCOAD', 'CYC1CCDH', 'CYCOADH', 'MSUCOADH']
    m_h = model.metabolites.get_by_id('h_c')
    for r_id in reaction_ids:
        r = model.reactions.get_by_id(r_id)
        s = r.metabolites[m_h]
        r.add_metabolites({m_h:-s})



def _edit_reactions(model):
    m_h = model.metabolites.get_by_id('h_c')
    m_h2o = model.metabolites.get_by_id('h2o_c')
    m_acp = model.metabolites.get_by_id('ACP_c')
    m_apoACP = model.metabolites.get_by_id('apoACP_c')
    m_donor = model.metabolites.get_by_id('donor_c')
    m_acceptor = model.metabolites.get_by_id('acceptor_c')

    # 2AMUCSAD
    model.reactions.get_by_id('2AMUCSAD').add_metabolites({m_h:1})

    # PICDDHs
    r = model.reactions.get_by_id('PICDDHs')
    r.gene_reaction_rule = 's0001' # Spontaneous
    r.add_metabolites({m_h:1})

    # QUINHYDs
    r = model.reactions.get_by_id('QUINHYDs')
    r.gene_reaction_rule = 's0001' # Spontaneous
    r.add_metabolites({m_h:1})

    # ASCOROX
    r = model.reactions.get_by_id('ASCOROX')
    r.gene_reaction_rule = 's0001' # Spontaneous
    r.add_metabolites({m_h:-1})

    # ASCORPO
    r = model.reactions.get_by_id('ASCORPO')
    r.gene_reaction_rule = 's0001' # Spontaneous
    r.add_metabolites({m_h:1})

    # GUL14LACDH
    model.reactions.get_by_id('GUL14LACDH').add_metabolites({m_h:-2})

    # LACDHCYT
    model.reactions.get_by_id('LACDHCYT').add_metabolites({m_h:-2})

    # AMMQT9
    model.reactions.get_by_id('AMMQT9').add_metabolites({m_h:-1})

    # CELLUL_DEGe
    model.reactions.get_by_id('CELLUL_DEGe').add_metabolites({m_h2o:250})

    # ADXFUTNS
    model.reactions.get_by_id('ADXFUTNS').add_metabolites({m_h:-1})

    # CYSBLY
    model.reactions.get_by_id('CYSBLY').add_metabolites({m_h:-1})    

    # MCPST
    model.reactions.get_by_id('MCPST').add_metabolites({m_h:1})    

    # SRUBSYN
    r = model.reactions.get_by_id('SRUBSYN')
    r.add_metabolites({m_donor: -1, m_acceptor:1, m_h:1})

    # CHTNDG
    model.reactions.get_by_id('CHTNDG').add_metabolites({m_h:-1})    

    # O4OXOCYR
    model.reactions.get_by_id('O4OXOCYR').add_metabolites({m_h:1})    
    
    # ONTROOXR
    model.reactions.get_by_id('ONTROOXR').add_metabolites({m_h:1})    

    # KHK2
    r =  model.reactions.get_by_id('KHK2')
    r.annotation['biocyc'] = "META:LYXK-RXN"
    r.add_metabolites({m_h:1}) 

    # FORGLUIH
    model.reactions.get_by_id('FORGLUIH').add_metabolites({m_h:-1})

    # SAMAT
    model.reactions.get_by_id('SAMAT').add_metabolites({m_h:-1})   
    
    # MYCARSL 
    model.reactions.get_by_id('MYCARSL').add_metabolites({m_h:-1})

    # DPROLORED
    model.reactions.get_by_id('DPROLORED').add_metabolites({m_h:1})

    # DRIB2FUR
    model.reactions.get_by_id('DRIB2FUR').add_metabolites({m_h:1})
    
    # 25DBUT
    model.reactions.get_by_id('25DBUT').add_metabolites({m_h:-1})

    # 4HPROOR
    model.reactions.get_by_id('4HPROOR').add_metabolites({m_h:1})

    # AMMQT9
    model.reactions.get_by_id('AMMQT9').add_metabolites({m_h:2})

    # FeHemu
    model.reactions.get_by_id('FeHemu').add_metabolites({m_h2o:-1, m_h:5})

    # FEROXG1S
    # It is not clear what ionic state the feroxE is in, BioCyC uses 0, but the concensus in this model of 
    # other unloaded feroxamines is -3. Thus, we keep this and balance one reaction with 3 protons
    model.reactions.get_by_id('FEROXG1S').add_metabolites({m_h:3})



    # Delete HISpep - does not make sense and is anyway block
    model.remove_reactions("HISpep")

    # Delete DHNAOT9
    # This reaction is similar to DHNANT4, except that the reaction is not balanced
    # and dmmql9 is replaced by dmmq9. This is more likely similar to META:DMK-RXN
    model.remove_reactions("DHNAOT9")
    model.reactions.get_by_id('DHNANT4').annotation["biocyc"] = "META:DMK-RXN"
    


    # SCB11
    r = model.reactions.get_by_id('SCB11')
    r.add_metabolites({m_acp:1, m_apoACP:-1})

    r = model.reactions.get_by_id('SCB21b')
    r.add_metabolites({m_acp:1, m_apoACP:-1})

    r = model.reactions.get_by_id('SCB31')
    r.add_metabolites({m_acp:1, m_apoACP:-1})

def make_metabolite_charge_formula_csv(model, path):
    metabolite_ids = []
    metabolite_charge = []
    metabolite_formula = []

    for m in model.metabolites:
        metabolite_ids.append(m.id)
        metabolite_charge.append(m.charge)
        metabolite_formula.append(m.formula)

    df = pd.DataFrame()
    df["ID"] = metabolite_ids
    df["Charge"] = metabolite_charge
    df["Formula"] = metabolite_formula
    df.to_csv(path, index = False)

def print_all_reaction_balances(model):
    i = 0
    for r in model.reactions:
        try: 
            ann = r.annotation["biocyc"]
        except KeyError:
            ann = 'None'
        if len(r.metabolites)>1:
            mb = r.check_mass_balance()
            if len(mb):
                print(r.id, r.reaction, mb, ann)
                i+=1
    print(i)

def print_biocyc_ids_for_mets_without_formula(model):
    for m in model.metabolites:
        if not m.formula:
            try: 
                ann = m.annotation["biocyc"]
            except KeyError:
                print("{0}".format(m.id))
            else:
                print("{0},{1}".format(ann, m.id))

if __name__ == '__main__':
    if 1:
        # Load the latest model version, that your script aims to update
        model = export.get_earlier_model_unversioned('v1.3.1') 
        
        # Call the required subroutines or functions to curate the model
        fix_stoichiometric_consistency(model)
        
        # Export the model using the export() function
        export.export(model, formats = ["xml", "yml"], write_requirements = 1, objective = "BIOMASS_SCO_tRNA")
        cobra.io.validate_sbml_model(REPO_PATH + "/model/Sco-GEM.xml")

    else:
        # If you want to test curation without storing the models
        model = export.get_earlier_model_unversioned('v1.3.1')


        fix_stoichiometric_consistency(model)

        make_metabolite_charge_formula_csv(model, "../../data/curation/v1_3_2/v1_3_2_draft_metabolite_charge_formula.csv")

        print_all_reaction_balances(model)

        print_biocyc_ids_for_mets_without_formula(model)
 