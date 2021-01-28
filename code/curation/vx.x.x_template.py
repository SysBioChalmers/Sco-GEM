"""
Curations for release vx.x.x.

This script combines various curations on Sco-GEM vx.x.y to produce vx.x.x.
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

def add_metabolite_charge(model): # Name the subroutine related to the tasks performed
    # Fixes Issue #xx
    
    # Write here the code to do the curations related to a particular issue
    # Multiple issues can be addressed in the same function or subroutine, but
    # indicate within the code where which issue is fixed.
    m = model.metabolites[0]
    m.charge = -2 
    
if __name__ == '__main__':
    # Load the latest model version, that your script aims to update
    model = export.get_earlier_model_unversioned('v1.2.1') 
    
    # Call the required subroutines or functions to curate the model
    add_metabolite_charge(model)
    
    # Export the model using the export() function
    export.export(model, formats = ["xml", "yml"], write_requirements = 1, objective = "BIOMASS_SCO_tRNA")
    cobra.io.validate_sbml_model(REPO_PATH + "/model/Sco-GEM.xml")
