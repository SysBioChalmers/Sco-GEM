"""
Curations for release v1.3.1

This script combines various curations on Sco-GEM v1.3.0 to produce v1.3.1.
Indicated is which issues are solved, more detailed explanation is given in the
relevant pull requests.

"""
# Import required functions, add any that are necessary for your code.
import sys
import cobra
from dotenv import find_dotenv

# Find .env + define paths. REPO_PATH refers to the repository root folder, and
# can be used to conveniently define absolute paths.
REPO_PATH = find_dotenv()
REPO_PATH = REPO_PATH[:-5]

sys.path.append(REPO_PATH + "/code")
import export

def correct_annotations(model): # Name the subroutine related to the tasks performed
    # Fixes Issue #111
    
    # Resolves incorrect annotations that were not resolved in the previous release.
    # ec-codes, preferably retain one code per reaction, multiple ec-codes indicates
    # multi-step reaction, manually curated via MetaCyc, KEGG and https://www.qmul.ac.uk/sbcs/iubmb/enzyme/ 
    AMPTASECG = model.reactions.get_by_id("AMPTASECG")
    AMPTASECG.annotation['ec-code'] = '3.4.13.18'
    ASP1DC = model.reactions.get_by_id("ASP1DC")
    ASP1DC.annotation['ec-code'] = '4.1.1.11'
    CHOLD = model.reactions.get_by_id("CHOLD")
    CHOLD.annotation['ec-code'] = '1.1.99.1'
    CYSDS = model.reactions.get_by_id("CYSDS")
    CYSDS.annotation['ec-code'] = '4.4.1.28'
    DHFS = model.reactions.get_by_id("DHFS")
    DHFS.annotation['ec-code'] = '6.3.2.12'
    DURIPP = model.reactions.get_by_id("DURIPP")
    DURIPP.annotation['ec-code'] = '2.4.2.3'
    GLYCL = model.reactions.get_by_id("GLYCL")
    GLYCL.annotation['ec-code'] = '1.4.1.27'
    GTPDPDP = model.reactions.get_by_id("GTPDPDP")
    GTPDPDP.annotation['ec-code'] = '3.6.1.40'
    GUAPRT = model.reactions.get_by_id("GUAPRT")
    GUAPRT.annotation['ec-code'] = '2.4.2.8'
    MACCOAT = model.reactions.get_by_id("MACCOAT")
    MACCOAT.annotation['ec-code'] = '2.3.1.16'
    MBCOA2 = model.reactions.get_by_id("MBCOA2")
    MBCOA2.annotation['ec-code'] = '1.3.8.1'
    NADH10b = model.reactions.get_by_id("NADH10b")
    NADH10b.annotation['ec-code'] = '1.6.5.2'
    NADS1 = model.reactions.get_by_id("NADS1")
    NADS1.annotation['ec-code'] = '6.3.1.5'
    PGM = model.reactions.get_by_id("PGM")
    PGM.annotation['ec-code'] = '5.4.2.12'
    PUNP4 = model.reactions.get_by_id("PUNP4")
    PUNP4.annotation['ec-code'] = '2.4.2.1'
    PUNP6 = model.reactions.get_by_id("PUNP6")
    PUNP6.annotation['ec-code'] = '2.4.2.1'
    SERD_L = model.reactions.get_by_id("SERD_L")
    SERD_L.annotation['ec-code'] = '4.3.1.17'
    SSALy = model.reactions.get_by_id("SSALy")
    SSALy.annotation['ec-code'] = '1.2.1.79'
    _2OXOADOX = model.reactions.get_by_id("2OXOADOX")
    _2OXOADOX.annotation['ec-code'] = ['1.2.4.2','1.8.1.4','2.3.1.61']
    AKGDH = model.reactions.get_by_id("AKGDH")
    AKGDH.annotation['ec-code'] = ['1.2.4.2','1.8.1.4','2.3.1.61']
    NSHDDS = model.reactions.get_by_id("NSHDDS")
    NSHDDS.annotation['ec-code'] = '6.3.5.12'
    PAAT = model.reactions.get_by_id("PAAT")
    PAAT.annotation['ec-code'] = '2.6.1.113'
    # KEGG ids
    CARBMODH = model.reactions.get_by_id("CARBMODH")
    CARBMODH.annotation['kegg.reaction'] = 'R07157'
    GGGABADH = model.reactions.get_by_id("GGGABADH")
    GGGABADH.annotation['kegg.reaction'] = 'R07417'
    GGGABADH.annotation['ec-code'] = '1.2.1.99'
    S7PI = model.reactions.get_by_id("S7PI")
    S7PI.annotation['kegg.reaction'] = 'R09768'
    TYRMO = model.reactions.get_by_id("TYRMO")
    TYRMO.annotation['kegg.reaction'] = 'R00031'
    UAAGLS2 = model.reactions.get_by_id("UAAGLS2")
    UAAGLS2.annotation['kegg.reaction'] = 'R02786'
    # Remove unused metabolite
    model.metabolites.remove('g1p_B_c') 
    
if __name__ == '__main__':
    # Load the latest model version, that your script aims to update
    model = export.get_earlier_model_unversioned('v1.3.0') 
    
    # Call the required subroutines or functions to curate the model
    correct_annotations(model)
    
    # Export the model using the export() function
    export.export(model, formats = ["xml", "yml"], write_requirements = 1, objective = "BIOMASS_SCO_tRNA")
    cobra.io.validate_sbml_model(REPO_PATH + "/model/Sco-GEM.xml")
