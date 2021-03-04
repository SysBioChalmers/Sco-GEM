import pytest
from memote.utils import annotate, wrapper, truncate
from cobra.io import read_sbml_model
import os.path

from dotenv import find_dotenv
# find .env + define paths
REPO_PATH = find_dotenv()
REPO_PATH = REPO_PATH[:-5]

@pytest.fixture(scope = "session")
def model(model = None):
    if model is None:
        path = os.path.join(os.path.dirname(__file__), REPO_PATH + "/model/Sco-GEM.xml")
        print("Loading model {0}".format(path))
        return read_sbml_model(path)
    else:
        return model

# @annotate(title="Some human-readable descriptive title for the report", format_type="raw")
# def test_read_model(model):
#     """
#     Checks that the model can be read using cobrapy. Should be extended to test different versions
#     """
#     ann = test_read_model.annotation
#     ann["message"] = wrapper.fill("The model can be read by cobrapy succsessfuly")
#     assert bool(model), ann["message"]

@annotate(title="Test if all reactions have been given a name", format_type="percent")
def test_reaction_names(model):
    """
    Test reaction names
    """
    ann = test_reaction_names.annotation
    
    n_r = 0
    for r in model.reactions:
        if (not isinstance(r.name, str)) or (not len(r.name)):
            n_r += 1
    ann["metric"] = 1 - (n_r/len(model.reactions))
    ann["data"] = len(model.reactions) - n_r
    ann["message"] = wrapper.fill("{0} of {1} reactions don't have names".format(n_r, len(model.reactions)))
    assert n_r == 0, ann["message"]        

@annotate(title="Test if all metabolites have been given a name", format_type="percent")
def test_metabolite_names(model):
    """
    Test metabolite names
    """
    ann = test_metabolite_names.annotation
    
    n_m = 0
    for m in model.metabolites:
        if (not isinstance(m.name, str)) or (not len(m.name)):
            n_m += 1

    ann["metric"] = 1 - (n_m/len(model.metabolites))
    ann["data"] = len(model.metabolites) - n_m
    ann["message"] = "{0} of {1} metabolites don't have names".format(n_m, len(model.metabolites))
    assert n_m == 0, ann["message"]


@annotate(title="Test that the growth rate is around 0.075", format_type="number")
def test_growth_rate(model):
    """
    Test metabolite names
    """
    ann = test_growth_rate.annotation
    solution = model.optimize()
    
    ann["message"] = wrapper.fill(
        """The growth rate is {0}
        """.format(solution.objective_value))
    ann["data"] = solution.objective_value 
    ann["metric"] = int(0.071 <= solution.objective_value < 0.078)

    assert 0.071 <= solution.objective_value < 0.078, ann["message"]


@annotate(title="Test germicidinA production", format_type="number")
def test_germicidinA_production(model):
    """
    Test that the model can produce germicidin A
    """
    with model:
        model.reactions.DM_germicidinA_c.objective_coefficient = 1
        model.reactions.BIOMASS_SCO.objective_coefficient = 0
        solution = model.optimize()
    ann = test_germicidinA_production.annotation
    ann["message"] = "The production of germicidinA is {0}".format(solution.objective_value)
    ann["data"] = solution.objective_value 
    ann["metric"] = int(solution.objective_value > 1e-3)
    assert solution.objective_value > 1e-3, ann["message"]

@annotate(title="Test germicidinB production", format_type="number")
def test_germicidinB_production(model):
    """
    Test that the model can produce germicidin B
    """
    with model:
        model.reactions.DM_germicidinB_c.objective_coefficient = 1
        model.reactions.BIOMASS_SCO.objective_coefficient = 0
        solution = model.optimize()
    ann = test_germicidinB_production.annotation
    ann["message"] = "The production of germicidinB is {0}".format(solution.objective_value)
    ann["data"] = solution.objective_value 
    ann["metric"] = int(solution.objective_value > 1e-3)
    assert solution.objective_value > 1e-3, ann["message"]

@annotate(title="Test germicidinC production", format_type="number")
def test_germicidinC_production(model):
    """
    Test that the model can produce germicidin C, and report how much
    """
    with model:
        model.reactions.DM_germicidinC_c.objective_coefficient = 1
        model.reactions.BIOMASS_SCO.objective_coefficient = 0
        solution = model.optimize()
    ann = test_germicidinC_production.annotation
    ann["data"] = solution.objective_value
    ann["metric"] = int(solution.objective_value > 1e-3)
    ann["message"] = "The production of germicidinC is {0}".format(solution.objective_value)
    assert solution.objective_value > 1e-3, ann["message"]

@annotate(title="Test ACT production", format_type="number")
def test_ACT_production(model):
    """
    Test ACT production
    """
    with model:
        model.reactions.EX_act_e.objective_coefficient = 1
        model.reactions.BIOMASS_SCO.objective_coefficient = 0
        solution = model.optimize()
    ann = test_ACT_production.annotation
    ann["data"] = solution.objective_value
    ann["metric"] = int(solution.objective_value > 1e-3)
    ann["message"] = "The production of ACT is {0}".format(solution.objective_value)
    assert solution.objective_value > 1e-3, ann["message"]


@annotate(title="Test RED production", format_type="number")
def test_RED_production(model):
    """
    Test RED production
    """
    with model:
        model.reactions.DM_RED_c.objective_coefficient = 1
        model.reactions.BIOMASS_SCO.objective_coefficient = 0
        solution = model.optimize()
    ann = test_RED_production.annotation
    ann["data"] = solution.objective_value
    ann["metric"] = int(solution.objective_value > 1e-3)
    ann["message"] = "The production of RED is {0}".format(solution.objective_value)
    assert solution.objective_value > 1e-3, ann["message"]


@annotate(title="Test CDA production", format_type="number")
def test_CDA_production(model):
    """
    Test CDA production
    """
    with model:
        model.reactions.EX_CDA_e.objective_coefficient = 1
        model.reactions.BIOMASS_SCO.objective_coefficient = 0
        solution = model.optimize()
    ann = test_CDA_production.annotation
    ann["data"] = solution.objective_value
    ann["metric"] = int(solution.objective_value > 1e-3)
    ann["message"] = "The production of CDA is {0}".format(solution.objective_value)
    assert solution.objective_value > 1e-3, ann["message"]


