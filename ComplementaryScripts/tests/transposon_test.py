import pytest
from memote.utils import annotate, wrapper, truncate
from cobra.io import read_sbml_model
import os.path




@pytest.fixture(scope = "session")
def model():
    path = os.path.join(os.path.dirname(__file__), "../../ModelFiles/xml/scoGEM.xml")
    return read_sbml_model(path)

