import pytest
from memote.utils import annotate, wrapper, truncate
from cobra.io import read_sbml_model
import os.path
from pathlib import Path
import pandas as pd
import yaml
import re

REPO_ROOT_DIR = Path(__file__).parent.parent.parent
DATA_DIR = REPO_ROOT_DIR / "ComplementaryData"
SOLVER = "gurobi"
MAX_GLUCOSE_UPTAKE_RATE = -2.1
CORRESPONDING_AMMONIUM_UPTAKE = -1.85
NO_GROWTH_THRESHOLD = 0.028 # Corresponds to 1hr doubling time


@pytest.fixture(scope = "session")
def model():
    path = os.path.join(os.path.dirname(__file__), "../../ModelFiles/xml/scoGEM.xml")
    model = read_sbml_model(path)
    model.solver = SOLVER
    return model


@pytest.fixture(scope = "session")
def experiments():
    experiments_path = DATA_DIR / "experiments.yml"
    with experiments_path.open(mode = "r") as fid:
        experiments_yaml = yaml.load(fid)
    return experiments_yaml

# def growth_data():
    # DATA_DIR / 
    # df = pd.read_csv("")


@annotate(title = "Test growth in given environments", format_type = "raw")
def test_WT_growth_conditions(model, experiments):
    """
    Compare in silic growth with in vivo data for growth in different environements
    """
    # Set minimal media
    info = experiments["growth"]["experiments"]["wild_type_growth"]

    set_medium(model, info["medium"], experiments)

    growth_df = get_growth_data(info["data"])

    prediction_dict = {
        "FP": 0, # False positives
        "TP": 0, # True positives
        "FN": 0, # False negatives
        "TN": 0, # True negative
        }
    wrong_predictions = []
    for index, row in growth_df.iterrows():
        with model as M:
            carbon_exchanges = row["Carbon exchange reactions"].split(",")
            nitrogen_exchanges = row["Nitrogen exchange reactions"].split(",")
            M = set_carbon_and_nitrogen_uptake(M, carbon_exchanges+nitrogen_exchanges)
            try:
                s = M.optimize()
            except:
                s = None
            print(index, end = "\t")
            correct_prediction = classify_prediction(s, row["In vivo growth"], prediction_dict)
            if not correct_prediction:
                _print_wrong_prediction(row, s, M)
                wrong_predictions.append(str(index))

    FN_and_TN = (prediction_dict["TN"] >= 4) and (prediction_dict["FN"] == 0)
    FP_and_TP = (prediction_dict["FP"] <= 2) and (prediction_dict["TP"] >= 51)
    ann = test_WT_growth_conditions.annotation
    ann["message"] = wrapper.fill("""Wild-type growth (expected in parenthesis):\n
                                     TP: {0}(51)\t FN: {1}(0)\n 
                                     FP: {2}(2)\t  TN: {3}(4)\n
                                     The prediction is wrong for the following environments: {4}""".format(
                            *[prediction_dict[x] for x in ["TP", "FN", "FP", "TN"]], ", ".join(wrong_predictions)))

    assert FN_and_TN and FP_and_TP, ann["message"]

def _print_wrong_prediction(row, s, M):
    print("\nWrong prediction:")
    print(row)
    if s is not None:
        print("Growth rate: ", s.objective_value)
        print(M.summary())
    else:
        print("Growth rate: ", 0)


def classify_prediction(solution, in_vivo_growth, prediction_dict):
    if (solution is None) or (solution.status == "infeasible"):
        in_silico_growth = False 
    else:
        in_silico_growth = bool(solution.objective_value > NO_GROWTH_THRESHOLD)

    if in_vivo_growth:
        if in_silico_growth:
            prediction_dict["TP"] += 1
            correct_prediction = True
        else:
            prediction_dict["FN"] += 1
            correct_prediction = False

    else:
        if in_silico_growth:
            prediction_dict["FP"] += 1
            correct_prediction = False
        else:
            prediction_dict["TN"] += 1
            correct_prediction = True

    print("In vivo: {0}, in silico: {1}".format(in_vivo_growth, in_silico_growth), prediction_dict)
    return correct_prediction

def set_carbon_and_nitrogen_uptake(model, exchange_reactions):
    total_carbon_uptake = abs(MAX_GLUCOSE_UPTAKE_RATE) * 6 # Glucose has 6 carbons
    total_nitrogen_uptake = abs(CORRESPONDING_AMMONIUM_UPTAKE) * 1 # Ammonium has 1 nitrogen

    carbon_bound_exp = 0
    nitrogen_bound_exp = 0
    for r_id in exchange_reactions:
        r = model.reactions.get_by_id(r_id.strip())
        n_carbon, n_nitrogen = get_number_of_carbons_and_nitrogens_in_metabolite(list(r.metabolites)[0])
        r.lower_bound = -1000
        carbon_bound_exp += r.reverse_variable * n_carbon
        nitrogen_bound_exp += r.reverse_variable * n_nitrogen
    carbon_bound = model.problem.Constraint(carbon_bound_exp, ub = total_carbon_uptake)
    nitrogen_bound = model.problem.Constraint(nitrogen_bound_exp, ub = total_nitrogen_uptake)
    model.add_cons_vars(carbon_bound)
    model.add_cons_vars(nitrogen_bound)
    return model

def get_growth_data(fn):
    growth_path = DATA_DIR / "growth" / fn
    df = pd.read_csv(str(growth_path), sep = r";\s*", engine = "python")
    df.columns = [x.strip() for x in df.columns]
    df.set_index("N", inplace = True)
    return df

def set_medium(model, medium_name, experiments):
    """
    Set the uptake rates, 
    """
    medium_fn = experiments["medium"]["definitions"][medium_name]["filename"]
    medium_path = DATA_DIR / experiments["medium"]["path"] / medium_fn
    df = pd.read_csv(str(medium_path), sep = ",")
    df.columns = [x.strip() for x in df.columns]
    
    # Set all uptake rates to 0
    for r in model.boundary:
        r.lower_bound = 0
    
    # Set uptake rates according 
    for index, row in df.iterrows():
        model.reactions.get_by_id(row["exchange"]).lower_bound = row["uptake"]
    return model

def get_number_of_carbons_and_nitrogens_in_metabolite(metabolite):
    nitrogen_match = re.search(r"(N)(\d*)", metabolite.formula)
    if nitrogen_match is not None:
        if not nitrogen_match.group(2):
            n_nitrogen = 1
        else:
            n_nitrogen = int(nitrogen_match.group(2))
    else:
        n_nitrogen = 0

    carbon_match = re.search(r"(C)(\d*)", metabolite.formula)
    if carbon_match is not None:
        if not carbon_match.group(2):
            n_carbon = 1
        else:
            n_carbon = int(carbon_match.group(2))
    else:
        n_carbon = 0
    return n_carbon, n_nitrogen

