# -*- coding: utf-8 -*-
"""

Author: Snorre Sulheim 
Created: December 2018
Updated: 25.06.2019


This script is solving issue # 54, which is basically 3 tasks:
1) Change reaction bounds based on dG-values (change in gibbs free energy) predicted using eQuilibrator
2) Change lower bound to 0 for CPKS4a and CPKS4b to avoid unrealistic loop
3) Make most of the ATP-driven reactions irreversible, the table is given in issue #54. All ATP-driven reactions are forward unless
    a) Estimated positive dG
    b) nucleoside-diphosphate kinase, i.e. all NDPK* reactions (known to be reversible)
    This data is read from /ComplementaryData/curation/reversibility/reversibility_ATP_driven_reactions.csv,
    which contain the reversible ATP-driven reactions (before correction)
"""
from cobra.io import read_sbml_model, write_sbml_model
import pandas as pd
from collections import defaultdict
import numpy as np
import time
from matplotlib import pyplot as plt
import pickle
from itertools import combinations
import logging

KEY_TO_BOUNDS_DICT = {"forward": (0, 1000),
                      "backward": (-1000, 0),
                      "reversible": (-1000, 1000)}




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

def read_equilibrator_data2(fn, threshold = 30, consider_uncertainty = True, dG0_or_dGm = "dGm"):
    df = pd.read_csv(fn)

    if dG0_or_dGm == "dGm":
        key = "dGm_prime"
        key_std = "dGm_prime_std"
    else:
        key = "dG0_prime"
        key_std = "dG0_prime_std"

    # Remove rows with no equilibrator value
    df = df[df[key].notna()]
    
    df["dG"] = df[key]
    df["sigma dG"] = df[key_std] * 1.96 # Convert from std to 95% confidence interval

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


def create_model(model, constraining_dict, lethal_df, save_fn = None, model_id = None, skip = []):

    # Disregard all reactions in the first column of lethal_df
    lethal_reactions = list(set(list(lethal_df.loc[:, "rxn1"])+list(lethal_df.loc[:, "rxn2"])+list(lethal_df.loc[:, "rxn3"])))
    lethal_reactions += skip
    i = 0
    for r_id, bounds in constraining_dict.items():
        if r_id in lethal_reactions:
            continue
        i+=1
        r = model.reactions.get_by_id(r_id)
        logging.info("Change direction of reaction {0} from {1} to {2}".format(r_id, r.bounds, bounds))
        r.bounds = bounds
        logging.debug(model.optimize())

    logging.info("Changed the reversibility of {0} reactions in total".format(i))

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


def prep_df(model, df):
    constraining_dict = {}
    for i, row in df.iterrows():
        # Get model reaction
        r = model.reactions.get_by_id(row["Rxn"])
        db_reversibility = row["Reversibility in DB"]
        # Check the current reversibility of the reaction
        model_reversibility = check_reversibility(r)

        # Check that the reactions is not limited by 0 upper and lower bound
        if model_reversibility == "blocked":
            print("Obs! Reaction {0} is blocked in model. Skipping".format(r.id))
            continue
        if db_reversibility != model_reversibility:
            try:
                constraining_dict[r.id] = KEY_TO_BOUNDS_DICT[db_reversibility]
            except KeyError:
                print("No direction specified for reaction {0} with dG: {1}".format(r.id, row["dG"]))
                continue
    return constraining_dict

def check_smart_reversibility(model, reversibility_df, growth_rate_deviation = 0.01,
                              iterations1 = 1000, iterations2=1000, set_size = 6):
    """
    This function adds new reversibilities while ensuring growth and
    taking double and triple lethal pairs into account
    """

    # Seperate reactions into relaxing and constraining changes
    constraining_dict = prep_df(model, reversibility_df)

    test_random2(model, constraining_dict, growth_rate_deviation, iterations1 = iterations1, iterations2 = iterations2, set_size = set_size)

def check_smart_reversibility(model, constraining_dict, growth_rate_deviation = 0.01, set_size = 6, iterations1 = 1000, iterations2 = 1000, plot = True, always_lethal_ratio = 0.9):
    # initialize variables    
    r_id_list = list(constraining_dict.keys())
    N = len(r_id_list)
    possible_idxs = np.arange(N)

    set_array = np.zeros((iterations1+iterations2, N))
    solution_arr = np.zeros(iterations1+iterations2)
    start = time.time()

    growth_array = np.zeros(iterations1+iterations2)
    base_growth_rate = model.optimize().objective_value

    def run(model, idxs, i):
        for r_id in [r_id_list[x] for x in idxs]:
            model.reactions.get_by_id(r_id).bounds = constraining_dict[r_id]
        s = model.optimize().objective_value
        growth_array[i] = s
        dev = s - base_growth_rate
        if abs(dev) < growth_rate_deviation:
            solution_arr[i] = 0
        else:
            solution_arr[i] = np.sign(dev)                


    # Normal growth rate
    # s = model.optimize()
    for i in range(iterations1):
        idxs = np.random.choice(possible_idxs, set_size, False)
        set_array[i, idxs] = 1
        with model:
            run(model, idxs, i)

    # Remove always lethal from set
    lethal = set_array[solution_arr != 0, :]
    cols = set_array.sum(axis=0)*always_lethal_ratio <= lethal.sum(axis=0)
    possible_idxs = possible_idxs[~cols]
    print("Finished first iterations")

    # Continue run
    for i in range(iterations1, iterations2):
        idxs = np.random.choice(possible_idxs, set_size, False)
        set_array[i, idxs] = 1
        with model:
            run(model, idxs, i)
    print("Ran {0}+{1} iterations in {2:.0f} seconds".format(iterations1, iterations2, time.time()-start))

    with open("temp_{0}.pkl".format(time.strftime("%d%H%M")), "wb") as f:
        pickle.dump([set_array, solution_arr, r_id_list, constraining_dict, growth_array, always_lethal_ratio, growth_rate_deviation], f)

def analyse_random2(model, pkl_filename, save_fn,  plot = True):
    with open(pkl_filename, "rb") as f:
        set_array, solution_arr, r_id_list, constraining_dict, growth_array, always_lethal_ratio, growth_rate_deviation = pickle.load(f)

    lethal = set_array[solution_arr != 0, :]
    solution_lethal = solution_arr[solution_arr != 0]
    cols = set_array.sum(axis=0)*always_lethal_ratio <= lethal.sum(axis=0)
    rows = (set_array[:, cols] == 1).any(axis=1)
    set_array_except = set_array[~rows, :]
    solution_except = solution_arr[~rows]
    lethal_except = set_array_except[solution_except != 0, :]
    lethal_pairs, lethal_leftover, solution_leftover = analyze_lethal_pairs(set_array_except, solution_except, lethal_except)
    new_lethals = find_lethals_in_leftovers(model, lethal_leftover, solution_leftover, r_id_list, growth_rate_deviation, constraining_dict)

    # Store as csv
    df = pd.DataFrame(columns = ["rxn1","rxn1_lb","rxn1_ub","rxn2","rxn2_lb","rxn2_ub","rxn3","rxn3_lb","rxn3_ub","model growth"])
    
    for idx, j in enumerate(np.where(cols)[0]):
        r = r_id_list[j]
        df.loc[idx, "rxn1"] = r
        df.loc[idx, "rxn1_lb"] = constraining_dict[r][0]
        df.loc[idx, "rxn1_ub"] = constraining_dict[r][1]
    
    for pair in lethal_pairs:
        idx += 1
        r1 = r_id_list[pair[0]]
        r2 = r_id_list[pair[1]]
        
        df.loc[idx, "rxn1"] = r1
        df.loc[idx, "rxn1_lb"] = constraining_dict[r1][0]
        df.loc[idx, "rxn1_ub"] = constraining_dict[r1][1]

        df.loc[idx, "rxn2"] = r2
        df.loc[idx, "rxn2_lb"] = constraining_dict[r2][0]
        df.loc[idx, "rxn2_ub"] = constraining_dict[r2][1]
    
    for pair in new_lethals:
        idx += 1
        r1 = pair[0]
        r2 = pair[1]
        
        df.loc[idx, "rxn1"] = r1
        df.loc[idx, "rxn1_lb"] = constraining_dict[r1][0]
        df.loc[idx, "rxn1_ub"] = constraining_dict[r1][1]

        df.loc[idx, "rxn2"] = r2
        df.loc[idx, "rxn2_lb"] = constraining_dict[r2][0]
        df.loc[idx, "rxn2_ub"] = constraining_dict[r2][1]
        
        if len(pair) > 2:
            r3 = pair[2]
            df.loc[idx, "rxn3"] = r3
            df.loc[idx, "rxn3_lb"] = constraining_dict[r3][0]
            df.loc[idx, "rxn3_ub"] = constraining_dict[r3][1]

        if len(pair) > 3:
            print("Four-combination: ", pair)

    df.to_csv(save_fn, index = False)

    if plot:
        x = np.arange(set_array.shape[1])
        fig, ax = plt.subplots(1)

        ratio_lethal = lethal.sum(axis=0) / set_array.sum(axis = 0)
        ratio_pos = set_array[solution_arr == 1, :].sum(axis=0) / set_array.sum(axis = 0)
        ratio_neg = set_array[solution_arr == -1, :].sum(axis=0) / set_array.sum(axis = 0)
        
        ax.bar(x, height = ratio_lethal, color = "b", width = 0.4,  align = "edge", label = "Lethal")
        ax.bar(x+.4, height = ratio_pos, color = "g", width = 0.4,  align = "edge", label = "Too high growth")
        ax.bar(x+.4, height = ratio_neg, bottom = ratio_pos, color = "r", width = 0.4,  align = "edge", label = "Too low growth")
        # plt.xticks(x, r_id_list, rotation = 90)
        plt.show()

        fig, ax = plt.subplots(1)
        ax.bar(x, height = set_array.sum(axis=0), color = "b", width = 0.4,  align = "edge", label = "Lethal")
        plt.show()


   

def find_lethals_in_leftovers(model, lethal_leftover, solution_leftover, r_id_list, growth_rate_deviation, constraining_dict):
    s = model.optimize()
    base_growth_rate = s.objective_value
    lethals = []
    new_solutions = []
    print(lethal_leftover.shape, solution_leftover.shape, solution_leftover)
    
    def _test_model(model):
        s = model.optimize()
        dev = s.objective_value - base_growth_rate
        if abs(dev) < growth_rate_deviation:
            return False
        else:
            if np.sign(dev) != solution_leftover[row]:
                print(j, "Solution, but not a solution for {0}, {1}, {2}".format(pair, dev, solution_leftover[row]))
            else:
                print(j, "Solution for {0}".format(pair))

            lethals.append(pair)
            new_solutions.append(np.sign(dev))
            return True
    

    for row in range(lethal_leftover.shape[0]):
        arr = lethal_leftover[row, :]
        row_reaction_ids = [r_id_list[i] for i, x in enumerate(arr) if x]
        print(row, row_reaction_ids)

        # Test singles

        solution = False
        for n in [1, 2, 3, 4]:
            for j, pair in enumerate(combinations(row_reaction_ids, n)):
                with model:
                    for r_id in pair:
                        model.reactions.get_by_id(r_id).bounds = constraining_dict[r_id]

                    solution = _test_model(model)
                if solution:
                    break
            if solution:
                break
        if not solution:                
            print("No solution")

    lethals = list(set(lethals))

    return lethals



def analyze_lethal_pairs(set_array_except, solution_except, lethal_except):

    temp_arr = set_array_except.copy()
    temp_solution = solution_except.copy()
    lethal_temp = lethal_except.copy()
    order_n_lethals = np.argsort(lethal_temp.sum(axis=0))[::-1]
    lethal_pairs = []
    i = 0
    j = 1

    while True:
        idx_i = order_n_lethals[i]
        idx_j = order_n_lethals[j]
        n_i = lethal_except.sum(axis=0)[i]
        if n_i == 1:
            break
        # Find rows where both of these reactions are selected
        rows = (temp_arr[:, [idx_i, idx_j]] == 1).all(axis = 1)
        # Get the solution array for these rows
        rows_solution = temp_solution[rows]
        # If the solution is != 0 for all these rows we considere these two a lethal pair
        if (not (rows_solution == 0).any()) and (sum(rows)>2):
            print(lethal_temp.sum(axis=0)[idx_i], lethal_temp.sum(axis=0)[idx_j])
            print("n rows: ", sum(rows))
            # print(i,j, lethal_pairs, lethal_temp.shape, temp_arr.shape)
            lethal_pairs.append((idx_i,idx_j))
            temp_arr, temp_solution = strip_rows(temp_arr, temp_solution, rows)
            lethal_temp = temp_arr[temp_solution != 0, :]
            order_n_lethals = np.argsort(lethal_temp.sum(axis=0))[::-1]
            i = 0
            j = 1
        else:
            if j == len(order_n_lethals)-1:
                if i == len(order_n_lethals)-1:
                    break
                else:
                    i+=1
            else:
                j+=1
    if lethal_temp.shape[0]>0:
        print(lethal_temp.sum(axis=0))
    lethal_solution_temp = temp_solution[temp_solution != 0]
    return lethal_pairs, lethal_temp, lethal_solution_temp

def strip_rows(arr, solution, rows):
    return arr[~rows, :], solution[~rows]

def fill_reversibility_csv(fn, model):
    """
    An interesting observation is that there are only a few different
    growth values for this table.
    """
    df = pd.read_csv(fn, sep = ",", header = 0)
    for i in df.index:
        row  = df.loc[i,:]
        with model:
            for key in ["rxn1", "rxn2", "rxn3"]:
                r_id = df.loc[i, key]
                if pd.notna(r_id):
                    bounds = (df.loc[i, "{0}_lb".format(key)], df.loc[i, "{0}_ub".format(key)])
                    model.reactions.get_by_id(r_id).bounds = bounds
            s = model.optimize()
            df.loc[i, "model growth"] = np.round(s.objective_value, 4)
    df.to_csv(fn, index = False)
    print(df)

def change_bounds_according_to_eQuilibrator(model, equilibrator_data_fn, eq_lethals_fn):
    """
    Change bounds according to eQuilibrator dG - values, with abs(dGm_prime) < 30 as threshold
     for reversible reactions.

    The values are pre-calculated using the functions (3,4 and 5), and
    the total pipeline is then:
    1) read_equilibrator_data2
    2) prep_df
    3) check_smart_reversibility
    4) analyse_random2
    5) fill_reversibility
    6) create_model

    """
    df = read_equilibrator_data2(equilibrator_data_fn, threshold = 30, consider_uncertainty = False)
    lethal_df = pd.read_csv(eq_lethals_fn, sep = ",", header = 0)
    
    # discard reactions that have similar bounds
    reversibility_dict = prep_df(model, df)

    # List of reaction discarded from bounds being changed    
    skip = ["PROD2"] # Cause false positive prediction with proline https://www.ncbi.nlm.nih.gov/pubmed/22201764 
    skip += ["ARGSS", "OCT"] # Cause false positive prediction with arginine. The latter gas dG -28. The first one sets up loop
    skip += ["URIK1", "URIK2", "UPPRT"] # Cause false positive prediction with SCO5626 mutant.
    model = create_model(model, reversibility_dict, lethal_df, skip = skip)
    revert_backward_reactions(model)

def revert_backward_reactions(model):
    for reaction in model.reactions:
        if (reaction.upper_bound == 0) and (reaction.lower_bound < 0):
            old_bounds = reaction.bounds
            old_string = reaction.reaction
            reverse_reaction_dict = {}
            for m, coeff in reaction.metabolites.items():
                reverse_reaction_dict[m] = -2*coeff
            reaction.add_metabolites(reverse_reaction_dict)            
            reaction.bounds = (0, -old_bounds[0])

            print("{0}: Changed bounds from {1} to {2}".format(reaction.id, old_bounds, reaction.bounds))
            print("{0}: Changed reaction direction, i.e. from {1} to {2}".format(reaction.id, old_string, reaction.reaction))
            

def change_lower_bound_on_CPKS_reactions(model):
    model.reactions.CPKS4a.lower_bound = 0
    model.reactions.CPKS4b.lower_bound = 0
    logging.info("Changed lower bound of CPKS4a and CPKS4b to 0 to avoid loop")

def change_bounds_on_ATP_driven_reactions(model, ATP_driven_reactions_fn = "../../ComplementaryData/curation/reversibility/reversibility_ATP_driven_reactions.csv"):
    df = pd.read_csv(ATP_driven_reactions_fn, index_col = 0, sep = ";")
    for r_id, row in df.iterrows():
        r = model.reactions.get_by_id(r_id)
        bounds = KEY_TO_BOUNDS_DICT[str(row["Direction"]).lower()]
        if bounds[0] != r.lower_bound:
            logging.info("Change lower bound of {0} from {1} to {2}".format(r_id, r.lower_bound, r.bounds[0]))
            r.lower_bound = bounds[0]

    


if __name__ == '__main__':
    scoGEM = read_sbml_model("../../ModelFiles/xml/scoGEM.xml")
    if 0:
        # Pipeline for metacyc data
        metacyc_data_fn = "../../ComplementaryData/curation/reversibility/Reversibility-based-model-MetaCyc.csv"
        df = read_metacyc_data(metacyc_data_fn)
        constraining_dict = prep_df(scoGEM, df)
        # check_smart_reversibility(scoGEM, constraining_dict, iterations1 = 1000, iterations2 = 20000, set_size = 10)
        metacyc_lethals_fn = "../../ComplementaryData/curation/reversibility/metacyc_reversibility_lethals.csv"
        analyse_random2(scoGEM, "temp_071028.pkl", metacyc_lethals_fn)
        fill_reversibility_csv(metacyc_lethals_fn, scoGEM)
        lethal_df = pd.read_csv(metacyc_lethals_fn, sep = ",", header = 0)
        save_fn = "../../ModelFiles/xml/metacycScoGEM.xml"
        create_model(scoGEM, constraining_dict, lethal_df, save_fn = save_fn, model_id = "metacycScoGEM")
        
   
    if 0:
        # Pipeline for eQuilibrator
        equilibrator_data_fn = "../../ComplementaryData/curation/reversibility/eQuilibrator_reversibility.csv"
        df = read_equilibrator_data2(equilibrator_data_fn, threshold = 30, consider_uncertainty = False)
        constraining_dict = prep_df(scoGEM, df)
        print(len(df.index))
        print(len(constraining_dict))
        # check_smart_reversibility(scoGEM, constraining_dict, iterations1 = 1000, iterations2 = 20000, set_size = 10)
        eq_lethals_fn =  "../../ComplementaryData/curation/reversibility/eQuilibrator_reversibility_lethals.csv"
        # analyse_random2(scoGEM, "temp_111313.pkl", eq_lethals_fn)
        lethal_df = pd.read_csv(eq_lethals_fn, sep = ",", header = 0)
        fill_reversibility_csv(eq_lethals_fn, scoGEM)
        save_fn = "../../ModelFiles/xml/eqScoGEM.xml"
        skip = ["PROD2"] # Cause false positive prediction with proline https://www.ncbi.nlm.nih.gov/pubmed/22201764 
        skip += ["ARGSS", "OCT"] # Cause false positive prediction with arginine. The latter gas dG -28. The first one sets up loop
        skip += ["URIK1", "URIK2", "UPPRT"] # Cause false positive prediction with SCO5626 mutant.
        # skip LEUTA, VALTA, ILETA?
        create_model(scoGEM, constraining_dict, lethal_df, save_fn = save_fn, model_id = "eqScoGEM", skip = skip)

    if 0:
        # Print dG for a few reactions
        equilibrator_data_fn = "../../ComplementaryData/curation/reversibility/eQuilibrator_reversibility.csv"
        df = read_equilibrator_data2(equilibrator_data_fn, threshold = 30, consider_uncertainty = False)
        ALL_ATP_DRIVEN_REACTIONS = [
        "ADNCYC","GALKr","PGK","PYK","CYTK2","PRPPS","CDPMEK","NDPK9","CBIAT","UGMDDS","DBTS","DPCOAK","PTPATi","ADNK1","SHKK","PACCOAL","CBLAT","ACS","D5KGK","MCCC","DHAK","ACS2","UMPK","CBPS","NMNAT",
        "ALALAC","TMDK1","PYDXK","ASPK","BUTCOAL","NDPK5","UAMAS","PPDK","PYDXNK","RMK","NDPK6","PPCK","NDPK7","PFK_2","SUCOAS","PYDAMK","PPK2r","NDPK8","NDPK1","PMPK","DGK1","ADK1","XYLK2","ACTCOALIG",
        "DADK","ADSK","HEX4","NDPK2","UAMAGS","GLU5K","PNTK","NDPK3","KHK2","PRAGSr","SUCBZL","UAAGDS","NDPK4","XYLK","ACKr","DTMPK","CTPS1","MMCOAL","FACOAL160","FMNAT","HPPK2","HMPK1","GLGC","DURIK1",
        "CYTDK1","LTHRK","ACGAMK","MACPAL","GK1i","DDGLK","RBFK","GLCOAS","MALTHIK","CYTK1"]
        for r in ATP_DRIVEN_REACTIONS:
            print(r)
            print(df.loc[df["Rxn"] == r, "dG"].values)
    if 1:
        change_bounds_on_ATP_driven_reactions(scoGEM)
