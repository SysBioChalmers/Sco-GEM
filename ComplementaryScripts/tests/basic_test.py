import pytest
from memote.utils import annotate, wrapper, truncate
from cobra.io import read_sbml_model
import os.path




@pytest.fixture(scope = "session")
def model():
    path = os.path.join(os.path.dirname(__file__), "../../ModelFiles/xml/scoGEM.xml")
    return read_sbml_model(path)

@annotate(title="Some human-readable descriptive title for the report", format_type="raw")
def test_read_model(model):
    """
    Checks that the model can be read using cobrapy. Should be extended to test different versions
    """
    ann = test_read_model.annotation
    ann["message"] = wrapper.fill("The model can be read by cobrapy succsessfuly")
    assert bool(model), ann["message"]

@annotate(title="Test if all reactions have been given a name", format_type="raw")
def test_reaction_names(model):
    """
    Test reaction names
    """
    ann = test_reaction_names.annotation
    
    n_r = 0
    for r in model.reactions:
        if not isinstance(r.name, str):
            n_r += 1
    ann["message"] = wrapper.fill("{0} of {1} reactions don't have names".format(n_r, len(model.reactions)))
    assert n_r == 0, ann["message"]        

@annotate(title="Test if all metabolites have been given a name", format_type="raw")
def test_metabolite_names(model):
    """
    Test metabolite names
    """
    ann = test_metabolite_names.annotation
    
    n_m = 0
    for m in model.metabolites:
        if not isinstance(m.name, str):
            n_m += 1
    ann["message"] = "{0} of {1} metabolites don't have names".format(n_m, len(model.metabolites))
    assert n_m == 0, ann["message"]


@annotate(title="Test that the growth rate is around 0.075", format_type="raw")
def test_growth_rate(model):
    """
    Test metabolite names
    """
    ann = test_growth_rate.annotation
    solution = model.optimize()
    
    ann["message"] = wrapper.fill(
        """The growth rate is {0}
        """.format(solution.objective_value))
    assert 0.073 <= solution.objective_value < 0.077, ann["message"]


@annotate(title="Test that the model can produce germicidinA", format_type="raw")
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
    assert solution.objective_value > 1e-3, ann["message"]

@annotate(title="Test that the model can produce germicidinB", format_type="raw")
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
    assert solution.objective_value > 1e-3, ann["message"]

@annotate(title="Test that the model can produce germicidinC", format_type="raw")
def test_germicidinC_production(model):
    """
    Test that the model can produce germicidin C
    """
    with model:
        model.reactions.DM_germicidinC_c.objective_coefficient = 1
        model.reactions.BIOMASS_SCO.objective_coefficient = 0
        solution = model.optimize()
    ann = test_germicidinC_production.annotation
    ann["message"] = "The production of germicidinC is {0}".format(solution.objective_value)
    assert solution.objective_value > 1e-3, ann["message"]

@annotate(title="Test that the model can produce ACT", format_type="raw")
def test_ACT_production(model):
    """
    Test that the model can produce ACT
    """
    with model:
        model.reactions.DM_ACT_c.objective_coefficient = 1
        model.reactions.BIOMASS_SCO.objective_coefficient = 0
        solution = model.optimize()
    ann = test_ACT_production.annotation
    ann["message"] = "The production of ACT is {0}".format(solution.objective_value)
    assert solution.objective_value > 1e-3, ann["message"]


# def test_your_custom_case(read_only_model):
# """
# Docstring that briefly outlines the test function.

# A more elaborate explanation of why this test is important, how it works,
# and the assumptions/ theory behind it. This can be more than one line.
# """
# ann = test_your_custom_case.annotation
# ann["data"] = list(your_support_module.specific_model_quality(read_only_model))
# ann["metric"] = len(ann["data"]) / len(read_only_model.total_model_quality)
# ann["message"] = wrapper.fill(
#     """A concise message that displays and explains the test results.
#     For instance, if data is a list of items the amount: {} and
#     percentage ({:.2%}) values can be recorded here, as well as an
#     excerpt of the list itself: {}""".format(
#     len(ann["data"]), ann['metric'], truncate(ann['data'])
#     ))
# )
# assert len(ann["data"]) == 0, ann["message"]