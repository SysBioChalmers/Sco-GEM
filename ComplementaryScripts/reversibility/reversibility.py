from cobra.io import read_sbml_model, write_sbml_model
import pandas as pd
from collections import defaultdict

def read_metacyc_data(fn):
    """
    Read the reversibility data obtained from metacyc 
    and stored in ComplementaryData/Reversibility-based-model-MetaCyc.csv
    There is no _forward_ or _backward_ irreversible, so all _irreversible_ 
    reactions are considered _forward irreversible
    """
    df = pd.read_csv(fn)

    # Remove whitespace
    df["Reversibility in MetaCyc"] = df["Reversibility in MetaCyc"].str.strip()
    
    # make lowercase
    df["Reversibility in MetaCyc"] = df["Reversibility in MetaCyc"].str.lower()

    # Replace empty strings by nan
    df["Reversibility in MetaCyc"] = df["Reversibility in MetaCyc"].replace("", pd.np.nan)

    # replace "irreversible" with "forward"
    df["Reversibility in MetaCyc"] = df["Reversibility in MetaCyc"].replace("irreversible", "forward")
    
    df = df[df["Reversibility in MetaCyc"].notna()]

    # Renaming 
    df.rename(columns = {"Reversibility in MetaCyc": "Reversibility in DB"}, inplace = True)

    return df
    
def read_equilibrator_data(fn, threshold = 30, consider_uncertainty = True):
    df = pd.read_csv(fn)

    # Remove rows with no equilibrator value
    df = df[df["eQuilibrator - ΔrG′m [kJ/mol] - 1mM concentration"].notna()]
    
    df["dG"], df["sigma dG"] = df["eQuilibrator - ΔrG′m [kJ/mol] - 1mM concentration"].str.split("±").str

    df["dG"] = pd.to_numeric(df["dG"])
    df["sigma dG"] = pd.to_numeric(df["sigma dG"])

    direction = []
    for i, row in df.iterrows():
        dG = row["dG"]
        
        if consider_uncertainty:
            if dG < 0:
                conservative_value = dG + row["sigma dG"]
            else:
                conservative_value = dG - row["sigma dG"]
        else:
            conservative_value = dG


        if conservative_value < -threshold:
            direction.append("forward")
        elif conservative_value > threshold:
            direction.append("backward")
        else:
            direction.append("reversible")
    df["Reversibility in DB"] = direction

    return df[["Rxn", "Rxn KEGG ID", "Reversibility in model", "dG", "sigma dG", "Reversibility in DB"]]


def apply_reversibilities(model, data_df, save_fn = None, model_id = None):
    count_dict = defaultdict(int)
    for i, row in df.iterrows():
        r_id = row["Rxn"]
        r = model.reactions.get_by_id(r_id)
        db_reversibility = row["Reversibility in DB"]
        model_reversibility = check_reversibility(r)
        if model_reversibility == "blocked":
            print("Obs! Reaction {0} is blocked in model. Skipping".format(r.id))
            continue

        if db_reversibility != model_reversibility:
            if db_reversibility == "reversible":
                r.bounds = (-1000, 1000)
            elif db_reversibility == "forward":
                r.bounds = (0, 1000)
            elif db_reversibility == "backward":
                r.bounds = (-1000, 0)
            else:
                print("No direction specified for reaction {0} with dG: {1}".format(r.id, row["dG"]))
                continue
            key = "{0} to {1}".format(model_reversibility, db_reversibility)
            count_dict[key] += 1
            print("Changed direction of reaction {0} from {1} to {2}".format(r.id, model_reversibility, db_reversibility))

    n = 0
    for key, value in count_dict.items():
        print("{0}: {1}".format(key, value))
        n += value
    print("Changed the reversibility of {0} reactions in total. Data were available for {1} reactions".format(n, len(df)))

    if model_id:
        model.id = model_id

    if save_fn:
        write_sbml_model(model, save_fn)
    return model

def check_reversibility(reaction):
    """
    Return:
     - *forward* if reaction is irreversible in forward direction
     - *backward* if reaction is irreversible in the rbackward (reverse) direction
     - *reversible* if the reaction is reversible
     - *blocked* if the reactions is blocked
    """
    if (reaction.lower_bound < 0) and (reaction.upper_bound == 0):
        return "backward"
    elif (reaction.lower_bound == 0) and (reaction.upper_bound > 0):
        return "forward"
    elif (reaction.lower_bound == 0) and (reaction.upper_bound == 0):
        return "blocked"
    else:
        return "reversible"



if __name__ == '__main__':
    scoGEM = read_sbml_model("../../ModelFiles/xml/scoGEM.xml")
    if 1:
        equilibrator_data_fn = "../../ComplementaryData/curation/Reversibility-based-model-Equilibrator.csv"
        df = read_equilibrator_data(equilibrator_data_fn, threshold = 30, consider_uncertainty = True)
        eqScoGEM = apply_reversibilities(scoGEM, df, save_fn = "../../ModelFiles/xml/eqScoGEM.xml", model_id = "eqScoGEM")

    if 1:
        metacyc_data_fn = "../../ComplementaryData/curation/Reversibility-based-model-MetaCyc.csv"
        df = read_metacyc_data(metacyc_data_fn)
        apply_reversibilities(scoGEM, df, save_fn = "../../ModelFiles/xml/metacycScoGEM.xml", model_id = "metacycScoGEM")