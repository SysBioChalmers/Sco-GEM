# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 10:45:54 2018

@author: kumelj
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  3 15:28:15 2017

@author: noore
"""
import cobra
import csv
import sys
import os
import inspect
import pandas as pd
import numpy as np
import scipy as sp
import SBtabTools
import settings

from cobra.io import read_sbml_model
from six import iteritems
from cobra import Reaction, Metabolite
from cobra.flux_analysis.parsimonious import optimize_minimal_flux
from cobra.manipulation.modify import convert_to_irreversible
from stfba import settings
from cobra.flux_analysis.loopless import construct_loopless_model

sys.path.append(os.path.expanduser('~/git/SBtab/python'))
sys.path.append(os.path.expanduser('https://github.com/SysBioChalmers/sco-GEM/tree/develNew/ModelFiles/xml/scoGEM.xml'))
cobra_model=model_scoGEM;

class stFBA(object):
    def __init__(self, cobra_model, config_fname='C:/Users/kumelj/Downloads/OptKnock_050218/ELADNOOR/data/stfba_config.tsv',
                 solver=settings.LP_SOLVER):
        """
        construct a semi-thermodynamic model

        This adds MILP constraints to prevent internal flux cycling and
        energy generating cycles, as done in XXXXXX.

        This must be solved with an MILP capable solver.

        """

        sbtab = SBtabTools.openMultipleSBtab(config_fname)

        formation_energy = sbtab[0].toDataFrame()
        formation_energy.set_index('Compound', inplace=True)
        formation_energy['Value'] = formation_energy['Value'].apply(float)
        formation_energy.rename(columns={'Value': 'formation_energy'},
                                inplace=True)
        concentrations = sbtab[1].toDataFrame()
        concentrations.set_index('Compound', inplace=True)
        concentrations['Concentration:Min'] = \
            np.log(concentrations['Concentration:Min'].apply(float))
        concentrations['Concentration:Max'] = \
            np.log(concentrations['Concentration:Max'].apply(float))
        df = pd.merge(formation_energy, concentrations,
                      left_index=True, right_index=True, how='outer')
        R = 8.314e-3  # kJ/mol/K
        T = 298  # K

        self.model = cobra_model.copy()
        self.solver = solver

        convert_to_irreversible(self.model)
        max_ub = max(self.model.reactions.list_attr("upper_bound"))
        kappa = 10000  # used for the 'kappa' parameter

        # a dict for storing S^T
        thermo_stoic = {metabolite.id: {}
                        for metabolite in self.model.metabolites}
        # Slice operator is so that we don't get newly added metabolites
        self.original_metabolites = self.model.metabolites[:]
        self.original_objective = self.model.objective

        for reaction in self.model.reactions[:]:
            # Boundary reactions are not subjected to these constraints
            if len(reaction._metabolites) == 1:
                continue
            # populate the S^T dict
            bound_id = "thermo_bound_" + reaction.id
            for met, stoic in iteritems(reaction._metabolites):
                thermo_stoic[met.id][bound_id] = stoic

            # encode the st-FBA constrainsts:
            #   1) v <= z*v_max
            #   2) S^T * dfG <= K - K*z
            # first, we create the boolean variable (z), in the form of a
            # special reactions that can have only two flux values (0 or 1)
            reaction_ind = Reaction(reaction.id + "_indicator")
            reaction_ind.variable_kind = "integer"
            reaction_ind.upper_bound = 1

            # For the first constraint we create a "fake" metabolite
            # in the indicator reaction it has a stoichiometric coeff = v_max
            # in the original reaction it has a stoichiometric coeff = -1
            # The mass-balance constraint is set to be >= 0
            # Mathematically, this translates to:
            #   -1 * v + v_max * z >= 0
            reaction_ub = Metabolite(reaction.id + "_ind_ub")
            reaction_ub._constraint_sense = "G"
            reaction.add_metabolites({reaction_ub: -1})
            reaction_ind.add_metabolites({reaction_ub: max_ub})

            # For the second constraint we create anotehr "fake" metabolite and
            # add it also to the indicator reaction with coefficint K
            # The mass-balance constraint is set to be <= K.
            # Mathematically, this translates to:
            #   K * z <= K
            # Later, we will add the reaction Gibbs energy to the left side.
            reaction_bound = Metabolite(bound_id)
            reaction_bound._constraint_sense = "L"
            reaction_bound._bound = kappa
            reaction_ind.add_metabolites({reaction_bound: kappa})
            self.model.add_reaction(reaction_ind)

        # Now we create the thermodynamic variables (dfG), each will be a new
        # reaction with the corresponding metabolite name.
        for metabolite in self.original_metabolites:
            metabolite_var = Reaction("thermo_var_" + metabolite.id)
            if metabolite.id in df.index:
                dfG0 = df.loc[metabolite.id, 'formation_energy']
                lb = dfG0 + R*T*df.loc[metabolite.id, 'Concentration:Min'] - 1e-5
                ub = dfG0 + R*T*df.loc[metabolite.id, 'Concentration:Max'] + 1e-5
                if lb > ub:
                    raise ValueError('''lower concentration bound is higher
                                        than upper bound''')
                metabolite_var.bounds = (lb, ub)
            else:
                # if there are no specific bounds for this metabolite, use a very
                # wide range for its formation energy
                metabolite_var.bounds = (-kappa, kappa)
            self.model.add_reaction(metabolite_var)

            # Finally, we add the stoichiometric coefficients of each reaction that
            # creates or uses this metabolite (i.e. the transposed stoichiometric
            # matrix).
            for k, v in iteritems(thermo_stoic[metabolite.id]):
                m = self.model.metabolites.get_by_id(k)
                metabolite_var.add_metabolites({m: v})

    def optimize(self):
        return self.model.optimize(solver=self.solver)

    @staticmethod
    def find_egc(cobra_model, already_irreversible=False):
        """
            try to locate EGCs by blocking all transport reactions and
            maximizing the objective of creating ATP from ADP + Pi
        """
        model = cobra_model.copy()
        obj = Reaction("EGC_tester")
        obj.add_metabolites({model.metabolites.atp_c: -1,
                             model.metabolites.h2o_c: -1,
                             model.metabolites.adp_c: 1,
                             model.metabolites.pi_c: 1})
        model.add_reaction(obj)
        model.objective = "EGC_tester"

        # disable exchange reactions and ATM maintenance
        for r in model.reactions:
            if r.id[0:3] == 'EX_':
                r.lower_bound = 0
                r.upper_bound = 0
        model.reactions.ATPM.lower_bound = 0
        model.reactions.ATPM.upper_bound = 0

        # protons are sometimes not balanced well, so we ignore them
        model.reactions.EX_h_e.lower_bound = -1000
        model.reactions.EX_h_e.upper_bound = 1000

        FBA_sol = model.optimize()
        if FBA_sol.status == 'optimal' and FBA_sol.f > 1e-9:
            FBA_sol = \
                optimize_minimal_flux(model,
                                      already_irreversible=already_irreversible)
            stFBA.print_solution(FBA_sol)
            return FBA_sol
        else:
            print('No EGCs found')
            return None

    @staticmethod
    def print_solution(sol):
        if sol is None:
            print("Linear problem is not feasible")
        elif sol.status == 'optimal':
            print("maximal ATP futile production: ", sol.f)
            for rid, flux in sol.x_dict.iteritems():
                if abs(flux) > 1e-9:
                    print('%20s: %6.2f' % (rid, flux))
        else:
            print("Linear problem is %s" % sol.status)

stFBA.find_egc(cobra_model, already_irreversible=False);
model=model_scoGEM;
fba_sol = model.optimize(solver=settings.LP_SOLVER)
print("FBA max yield: %.3f" % fba_sol.f)

#print "ll-FBA max yield: ",
loopless_model = construct_loopless_model(model)
llfba_sol = loopless_model.optimize(solver=settings.LP_SOLVER)
print("looplessFBA max yield: %.3f" % llfba_sol.f)



    
