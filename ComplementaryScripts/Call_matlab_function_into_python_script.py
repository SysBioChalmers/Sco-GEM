# -*- coding: utf-8 -*-
"""
Created on Sun May 05 18:20:15 2019

@author: kumelj
"""

import matlab.engine
eng=matlab.engine.start_matlab();
scoGEM=eng.automatically_adding_transporters_git;
