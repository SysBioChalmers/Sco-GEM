import cobra
from cobra.flux_analysis.sampling import OptGPSampler, ACHRSampler
import pandas as pd
from pathlib import Path
import argparse

REPO_DIR = Path(__file__).parent.parent.parent
XML_FILES_DIR = REPO_DIR / "ModelFiles" / "xml"
RESULTS_FOLDER = REPO_DIR / "ComplementaryData" / "ecmodel" / "random_sampling"

def get_model(strain, timepoint):
    fn = XML_FILES_DIR / "ec{0}_{1}.xml".format(strain, timepoint)
    try:
        model = cobra.io.read_sbml_model(str(fn))
    except IOError:
        print("Could not read file: {0}. \n The strain must be either M145 or M1152".format(str(fn)))
        exit()
    return model

def run_random_sampling(strain, timepoint, iterations = 10, processes = 1):
    print("Initating random sampling for {0} at {1}h for {2} iterations divided on {3} processes".format(strain, timepoint, iterations, processes))
    model = get_model(strain, timepoint)
    model.solver = "gurobi"
    print("Model imported. Preparing model by removing blocked reactions")
    prep_model(model)
    print("The growth is fixed to: ", model.reactions.BIOMASS_SCO_tRNA.bounds)
    # fba_solution = model.optimize()
    # print(model.summary())

    # optgp = ACHRSampler(model)
    optgp = OptGPSampler(model, processes = processes)
    print("Start sampling: {0} iterations on {1} cores".format(iterations, processes))
    s = optgp.sample(iterations)
    print(s.head())

    check_or_make_folder(RESULTS_FOLDER)
    save_name = "randomsampling_{0}_{1}h_{2}_iterations_{3}_processes.csv".format(strain, timepoint, iterations, processes)
    s.to_csv(str(RESULTS_FOLDER / save_name))

def check_or_make_folder(folder_path):
    if not folder_path.is_dir():
        folder_path.mkdir()

def prep_model(model):
    model.reactions.BIOMASS_SCO_tRNA.objective_coefficient = 1
    # model.reactions.ATPM.lower_bound = 2.55
    set_ATPM_lower_bound(model)
    remove_blocked_reactions(model)
    finite_upper_bound(model)

def set_ATPM_lower_bound(model):
    with model as temp:
        temp.reactions.BIOMASS_SCO_tRNA.objective_coefficient = 0
        temp.reactions.ATPM.objective_coefficient = 1
        temp.optimize()
        max_atpm = temp.reactions.ATPM.flux
    model.reactions.ATPM.lower_bound = 0.99 * max_atpm
        

def test_model(strain, timepoint):
    model = get_model(strain, timepoint)
    model.solver = "gurobi"
    print("Model imported. Preparing model by removing blocked reactions")
    prep_model(model)
    s = model.optimize()
    if s.status == 'infeasible':
        print("Status infeasible for model {0} at {1} hours".format(strain, timepoint))
    else:
        print("Status OK for model {0} at {1} hours".format(strain, timepoint))


def finite_upper_bound(model):
    n = 0
    for r in model.reactions:
        if r.upper_bound > 2e3:
            r.upper_bound = 1e3
            n+=1
    print("Changed upper bound on{0} reactions".format(n))

def remove_blocked_reactions(model):
    n_start = len(model.reactions)
    
    r_list = []
    for r in model.reactions:
        if r.bounds == (0,0):
            r_list.append(r)
    model.remove_reactions(r_list, True)
    

    while True:
        r_list = []
        for m in model.metabolites:
            if len(m.reactions) == 1:
                r = list(m.reactions)[0]
                r_list.append(r.id)
                r.delete()
        print(len(r_list))
        if not len(r_list):
            break
    n_1 = len(model.reactions)
    print("Removed {0} reactions".format(n_start-n_1))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Sample the solution space of an enzyme-constrained model")
    parser.add_argument("-s", "--strain", help = "Which strain (either M145 or M1152)", default = "M145", type = str)
    parser.add_argument("-t", "--time", help = "Which timepoint", type = int)
    parser.add_argument("-n", "--n_iter", help = "Number of iterations", type = int)
    parser.add_argument("-p", "--n_proc", help = "Number of processes", type = int)
    args = parser.parse_args()
    run_random_sampling(args.strain, args.time, args.n_iter, args.n_proc)
    




