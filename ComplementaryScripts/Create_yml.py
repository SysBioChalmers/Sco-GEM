# -*- coding: utf-8 -*-
"""
Created on Tue May 14 16:06:50 2019

@author: kumelj
"""


import cobra
import logging
from cobra.io import read_sbml_model
solver='gurobi';
scoGEM=cobra.io.read_sbml_model("C:/Users/kumelj/Downloads/OptKnock_050218/scoGEM_model_fix_transporters_140519/scoGEM.xml");
cobra.io.save_yaml_model(scoGEM, "C:/Users/kumelj/Downloads/OptKnock_050218/scoGEM_model_fix_transporters_140519/scoGEM.xml");
