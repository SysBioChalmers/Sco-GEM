import pytest
from memote.utils import annotate, wrapper, truncate
from cobra.io import read_sbml_model
import os.path
from pathlib import Path
import pandas as pd
import yaml
import re

from dotenv import find_dotenv
# find .env + define paths
REPO_PATH = find_dotenv()
REPO_PATH = REPO_PATH[:-5]

DATA_DIR = REPO_PATH + "/data/"
SOLVER = "gurobi"
MAX_GLUCOSE_UPTAKE_RATE = -2.1
CORRESPONDING_AMMONIUM_UPTAKE = -1.85
NO_GROWTH_THRESHOLD = 0.028 # Corresponds to 1hr doubling time
# TRANSPOSON_NO_GROWTH_THRESHOLD = 


@pytest.fixture(scope = "session")
def model():
    path = os.path.join(os.path.dirname(__file__), REPO_PATH + "/model/Sco-GEM.xml")
    model = read_sbml_model(path)
    try:
        model.solver = SOLVER
    except:
        pass
    return model


@pytest.fixture(scope = "session")
def experiments():
    experiments_path = DATA_DIR + "experiments.yml"
    open(experiments_path, "rt")
    with open(experiments_path, "rt") as fid:
        experiments_yaml = yaml.load(fid, Loader=yaml.FullLoader)
    return experiments_yaml


@annotate(title = "Test growth for knockout-mutants from the transposon mutagenesis study by Xu et al.(2017)",
          format_type = "number")
def test_transposon_mutant_growth(model, experiments, use_tRNA_biomass = True, minimal_medium = False):
    """
    Compare in silico growth with in vivo data for the knock-out mutants from the 
    transposon mutagenesis study by xu et al. (2017)
    
    # Reference
    Xu, Zhong, et al. "Large-scale transposition
    mutagenesis of Streptomyces coelicolor identifies 
    hundreds of genes influencing antibiotic biosynthesis." 
    Applied and environmental microbiology 83.6 (2017): e02889-16.

    # Parameters
        model: 
            SBML-model
        experiments:
            yaml-file with info about data etc
        use_tRNA_biomass: bool,
            Change the objective to the biomass-function where the tRNA-complexes are used in the biomass
        minimal_medium: bool,
            Assume growth in a minimal medium with glucose and ammonium as the sole carbon and nitrogen sources.


    """

    # Set growth medium
    if minimal_medium:
        set_medium(model, "minimal_medium", experiments)
    else:
        model = set_transposon_growth_medium(model, experiments)

    if use_tRNA_biomass:
        model.reactions.BIOMASS_SCO.objective_coefficient = 0
        model.reactions.BIOMASS_SCO_tRNA.objective_coefficient = 1

    # Make X737 - mutant (which was used in the study) with the ACT-cluster removed
    model = _make_X737_mutant(model)



    info = experiments["growth"]["experiments"]["transposon_mutant_growth"]
    mutant_growth_df = get_growth_data(info["data"])

    prediction_dict = {
        "FP": 0, # False positives
        "TP": 0, # True positives
        "FN": 0, # False negatives
        "TN": 0, # True negative
        }
    wrong_predictions = []
    wt_solution = model.optimize()
    growth_threshold = wt_solution.objective_value * 0.5



    i = 0
    for index, row in mutant_growth_df.iterrows():
        with model as M:
            try:
                M = _knock_out_genes(row["Gene"], M, how = "any")
            except KeyError:
                continue
            i += 1
            try:
                s = M.optimize()
            except:
                s = None
            print(index, row["Gene"], end = "\t")
            correct_prediction = classify_prediction(s, row["In vivo growth"], prediction_dict, growth_threshold)
            if not correct_prediction:
                wrong_predictions.append("{0}: {1}".format(index, row["Gene"]))
    

    FN_and_TN = (prediction_dict["TN"] >= 33) and (prediction_dict["FN"] <= 5)
    FP_and_TP = (prediction_dict["FP"] <= 27) and (prediction_dict["TP"] >= 72)
    ann = test_transposon_mutant_growth.annotation
    ann["message"] = wrapper.fill(
        """Growth predictions for {5} of the {6} knock-out mutants from the transposon study 
           \n(expected in parenthesis):
           \nTP: {0}(72)\t FN: {1}(5)
           \nFP: {2}(27)\t  TN: {3}(33)
           \nThe prediction is wrong for the following mutants: \n{4}""".format(
            *[prediction_dict[x] for x in ["TP", "FN", "FP", "TN"]], ", ".join(wrong_predictions), i, len(mutant_growth_df)))

    ann["metric"] = int(FN_and_TN and FP_and_TP)
    ann["data"] = (prediction_dict["TN"]+prediction_dict["TP"])/sum([prediction_dict[x] for x in ["TP", "FN", "FP", "TN"]])
    assert FN_and_TN and FP_and_TP, ann["message"]


def _make_X737_mutant(model):
    for gene in model.genes:
        try:
            g_number = int(gene.id.strip("SCO"))
        except ValueError:
            continue
        else:
            if (5070 < g_number) and (g_number < 5093):
                gene.knock_out()
    return model

def set_transposon_growth_medium(model, experiments):
    """
    Because the transposon study was run in an undefined medium we assume that all sources of carbon and nitrogen is available,
    but we use an upper limit on the total uptake of carbon and nitrogen
    """    
    set_medium(model, "minimal_medium_no_carbon_no_nitrogen", experiments)
    sources = []
    for r in model.exchanges:
        n_carbon, n_nitrogen =  get_number_of_carbons_and_nitrogens_in_metabolite(list(r.metabolites)[0])
        if n_carbon + n_nitrogen > 0:
            sources.append(r.id)
    model = set_carbon_and_nitrogen_uptake(model, sources)
    return model





@annotate(title = "Test growth for knockout-mutants from the litterature in given environments", format_type = "number")
def test_single_mutant_growth(model, experiments):
    """
    Compare in silico growth with in vivo data for different knock-out mutants found in the litterature
    in the given environments
    """
    # Set minimal media
    info = experiments["growth"]["experiments"]["single_mutant_growth"]
    set_medium(model, info["medium"], experiments)
    
    # Get in vivo data for mutant growth
    mutant_growth_df = get_growth_data(info["data"])

    prediction_dict = {
        "FP": 0, # False positives
        "TP": 0, # True positives
        "FN": 0, # False negatives
        "TN": 0, # True negative
        }
    wrong_predictions = []
    for index, row in mutant_growth_df.iterrows():
        with model as M:
            carbon_exchanges = row["Carbon exchange reactions"].split(",")
            nitrogen_exchanges = row["Nitrogen exchange reactions"].split(",")
            M = set_carbon_and_nitrogen_uptake(M, carbon_exchanges+nitrogen_exchanges)
            M = _knock_out_genes(row["Genes"], M)
            try:
                s = M.optimize()
            except:
                s = None
            print(index, end = "\t")
            correct_prediction = classify_prediction(s, row["In vivo growth"], prediction_dict, NO_GROWTH_THRESHOLD)
            if not correct_prediction:
                _print_wrong_prediction(row, s, M)
                wrong_predictions.append("{0}: {1}".format(index, row["Genes"]))

    FN_and_TN = (prediction_dict["TN"] >= 9) and (prediction_dict["FN"] <= 1)
    FP_and_TP = (prediction_dict["FP"] == 0) and (prediction_dict["TP"] >= 6)
    ann = test_single_mutant_growth.annotation
    ann["message"] = wrapper.fill("""Growth predictions for mutants from the litterature 
                                     \n(expected in parenthesis):
                                     \nTP: {0}(6)\t FN: {1}(1) 
                                     \nFP: {2}(0)\t  TN: {3}(9)
                                     \nThe prediction is wrong for the following mutants: \n{4}""".format(
                            *[prediction_dict[x] for x in ["TP", "FN", "FP", "TN"]], ", ".join(wrong_predictions)))
    ann["data"] = (prediction_dict["TN"]+prediction_dict["TP"])/sum([prediction_dict[x] for x in ["TP", "FN", "FP", "TN"]])
    ann["metric"] = (prediction_dict["TN"]+prediction_dict["TP"])/sum([prediction_dict[x] for x in ["TP", "FN", "FP", "TN"]])
    assert FN_and_TN and FP_and_TP, ann["message"]

def _knock_out_genes(genes_string, model, how = "all"):
    genes = [g.strip() for g in genes_string.split(",")]
    if how == "all":
        for g_id in genes:
            g = model.genes.get_by_id(g_id)
            g.knock_out()
    else:
        i = 0
        for g_id in genes:
            try:
                g = model.genes.get_by_id(g_id)
            except KeyError:
                continue
            else:
                g.knock_out()
                i+=1
        if not i:
            raise KeyError

    return model


    

@annotate(title = "Test growth for WT in given environments", format_type = "number")
def test_WT_growth_conditions(model, experiments):
    """
    Compare in silico growth with in vivo data for growth in different environements
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
            correct_prediction = classify_prediction(s, row["In vivo growth"], prediction_dict, NO_GROWTH_THRESHOLD)
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
    ann["data"] = (prediction_dict["TN"]+prediction_dict["TP"])/sum([prediction_dict[x] for x in ["TP", "FN", "FP", "TN"]])
    ann["metric"] = int(FN_and_TN and FP_and_TP)
    assert FN_and_TN and FP_and_TP, ann["message"]

def _print_wrong_prediction(row, s, M):
    print("\nWrong prediction:")
    print(row)
    if s is not None:
        print("Growth rate: ", s.objective_value)
        print(M.summary())
    else:
        print("Growth rate: ", 0)


def classify_prediction(solution, in_vivo_growth, prediction_dict, growth_threshold):
    if (solution is None) or (solution.status == "infeasible"):
        in_silico_growth = False 
    else:
        in_silico_growth = bool(solution.objective_value > growth_threshold)

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
    growth_path = DATA_DIR + "growth/" + fn
    df = pd.read_csv(str(growth_path), sep = r";\s*", engine = "python")
    df.columns = [x.strip() for x in df.columns]
    df.set_index("N", inplace = True)
    return df

def set_medium(model, medium_name, experiments):
    """
    Set the uptake rates, 
    """
    medium_fn = experiments["medium"]["definitions"][medium_name]["filename"]
    medium_path = DATA_DIR + experiments["medium"]["path"] + "/" + medium_fn
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

# if __name__ == '__main__':
#     M = model()
#     E = experiments()
#     test_single_mutant_growth(M,E)