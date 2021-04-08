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
    fn = "../../data/curation/v1_3_2_curated_metabolite_charge_formula.csv"
    df = pd.read_csv(fn, header = 0)
    for _, row in df.iterrows():
        m = model.metabolites.get_by_id(row["ID"])
        m.charge = row["Charge"]
        m.formula = row["Formula"]

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
        if len(r.metabolites)>1:
            mb = r.check_mass_balance()
            if len(mb):
                print(r.id, r.reaction, mb)
                i+=1
    print(i)


if __name__ == '__main__':
    if 0:
        # Load the latest model version, that your script aims to update
        model = export.get_earlier_model_unversioned('v1.3.1') 
        
        # Call the required subroutines or functions to curate the model
        fix_stoichiometric_consistency()
        
        # Export the model using the export() function
        export.export(model, formats = ["xml", "yml"], write_requirements = 1, objective = "BIOMASS_SCO_tRNA")
        cobra.io.validate_sbml_model(REPO_PATH + "/model/Sco-GEM.xml")

    else:
        model = export.get_earlier_model_unversioned('v1.3.1')


        fix_stoichiometric_consistency(model)

        make_metabolite_charge_formula_csv(model, "../../data/curation/v1_3_2_draft_metabolite_charge_formula.csv")

        print_all_reaction_balances(model)
 