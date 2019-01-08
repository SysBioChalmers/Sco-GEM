"""
Author: Snorre Sulheim
Date: 27.11.2018
"""

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from scipy.misc import derivative
import matplotlib
from matplotlib import pyplot as plt
from collections import defaultdict, OrderedDict
from pathlib import Path
from cobra.io import read_sbml_model

plt.style.use("seaborn-white")
# matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams["axes.labelsize"] = 12
matplotlib.rcParams["errorbar.capsize"] = 2
pd.set_option("display.max_columns", 10)


MOLAR_MASS_DICT = {"Glucose": 180.156, #g/mol
                   "Glutamic acid": 147.13,  #g/mol
                   "Undecylprodigiosin 2": 393.575, # g/mol
                   "TBP":   634.546, # g/mol (mol weight for Actinorhodin)
                   "Germicidin-A": 196.246, # g/mol
                   "Germicidin-B": 182.219, #g/mol

                    }

CARBON_DICT = {
    "Glucose": 6, 
    "Glutamic acid": 5,  
    "Undecylprodigiosin 2": 25, 
    "Germicidin-A": 10, 
    "Germicidin-B": 10,
    "CO2": 1, 
}

MODEL_REACTION_DICT = {
    "Glucose": "EX_glc__D_e", 
    "Glutamic acid": "EX_glu__L_e",  
    "Undecylprodigiosin 2": "DM_RED_c", 
    "Germicidin-A": "DM_germicidinA_c", 
    "Germicidin-B": "DM_germicidinB_c",
    "CO2": "EX_co2_e", 
}
RED_M145_OD = [0.1, 0.0, 0.0, 0.0, 0.1, 0.4, 0.7, 1.4, 2.1, 3.7]

REPO_MAIN_FOLDER = Path(__file__).resolve().parent.parent
GROWTH_DATA_FOLDER = REPO_MAIN_FOLDER / "ComplementaryData"/ "growth"
scoGEM_FN = REPO_MAIN_FOLDER / "ModelFiles" / "xml" / "scoGEM.xml"

class RateEstimator(object):
    def __init__(self, strain_name, online_fn, offline_fn, offline_std_fn = None, proteome_timepoints = None):
        self.online_df = read_online(online_fn)
        self.offline_df = read_offline(offline_fn)
        self.offline_std_df = read_offline(offline_std_fn)
        self.phases = dict()
        self.proteome_timepoints = proteome_timepoints
        self.strain_name = strain_name

        # Results
        self.exponential_fit_CO2 = None
        self.linear_fit_CDW = OrderedDict()
        self.growth_rates = OrderedDict()
        self.substrate_fits = defaultdict(OrderedDict)
        self.substrate_rates = defaultdict(OrderedDict)
        self.max_growth_rate = None
        self.estimated_cdw = OrderedDict()
        self.rate_df = None

    def set_phase(self, name, start, stop):
        self.phases[name] = (start, stop)

    def fit_exponential_phase_to_CO2(self, phase_name = None, phase_limits = None):
        if phase_name:
            phase_limits = self.phases[phase_name]
        elif phase_limits:
            phase_name = "{0}-{1}".format(*phase_limits)
        else:
            print("No time limits given")

        df = self.online_df[(phase_limits[0] <= self.online_df["TAI"]) & (self.online_df["TAI"] <= phase_limits[1])]
        df = df[df["CO2"].notna()]

        time_arr = np.array(df["TAI"])
        CO2_arr = np.array(df["CO2"])

        popt, _ = curve_fit(exp_fun, time_arr, CO2_arr)
        fit_time = np.arange(phase_limits[0], phase_limits[1]+0.5, 0.5)
        fit = exp_fun(fit_time, *popt)
        self.exponential_fit_CO2 = {"time": fit_time, "fit": fit, "popt": popt}
        print("Fitted CO2: x0*e^ut: x0:{0:.2f}, u:{1:.2f}".format(*popt))
        self.max_growth_rate = popt[1]


    def plot_exponential_fit_CO2(self, save_name = None):
        if not self.exponential_fit_CO2:
            print("Run fit_exponential_phase_to_CO2 first")
            return False

        time = self.exponential_fit_CO2["time"]
        fit  = self.exponential_fit_CO2["fit"]
        popt = self.exponential_fit_CO2["popt"]

        fig, ax = plt.subplots(1)
        ax.plot(self.online_df["TAI"], self.online_df["CO2"], ".", label = "CO2 measurements")
        ax.plot(time, fit, lw = 2, label = "Fit: {0:.2f}*e^{1:.2f}t".format(*popt))
        ax.legend()
        ax.set_xlabel("Time after inoculation [h]")
        ax.set_ylabel("CO2 [mmol/L/h]")
        ax.set_title("Growth rate in exponential phase {0}".format(self.strain_name))
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        if save_name:
            fig.savefig(save_name, format = "svg")
        plt.show()


    def fit_piecewise_linear_CDW(self, phase1, phase2, phase3):
        tmin = self.phases[phase1][0]
        x0 = self.phases[phase1][1]
        x1 = self.phases[phase2][1]
        tmax = self.phases[phase3][1]
        
        df = self.offline_df[(tmin <= self.offline_df["TAI"]) & (self.offline_df["TAI"] <= tmax)]
        df = df[df["CDW"].notna()]
        popt, _ = curve_fit(piecewise_linear_2, np.array(df["TAI"]), np.array(df["CDW"]), bounds = ([33, 3, 0,0], [50, 8, 0.7, 0.05]), method = "dogbox")
        time_arr = np.arange(tmin, tmax, 0.1)
        fit = piecewise_linear_2(time_arr, *popt)
        plt.plot(df["TAI"], df["CDW"], ".")
        plt.plot(time_arr, fit)
        plt.show()

    def get_offline_std(self, key, phase_limits = None, mean_df = None):
        if not phase_limits:
            phase_limits = (0, 1e4)
        try:
            df = self.offline_std_df[self.offline_std_df[key].notna()]
        except TypeError:
            print("Standard deviation data not provided")
            return np.ones(len(mean_df)), False
        except KeyError:
            print("Wrong key for standard deviation: {0}".format(key))
            return np.ones(len(mean_df)), False
        else:
            df = df[(phase_limits[0] <= df["TAI"]) & (df["TAI"] <= phase_limits[1])]
            return np.array(df[key]), True


    def fit_linear_CDW_for_phases(self, phase_names):
        for phase_name in phase_names:
            self.fit_linear_CDW(phase_name)

    def fit_linear_CDW(self, phase_name = None, phase_limits = None):
        if phase_name:
            phase_limits = self.phases[phase_name]
        elif phase_limits:
            phase_name = "{0}-{1}".format(*phase_limits)
        else:
            print("No time limits given")

        cdw_df = self.offline_df[self.offline_df["CDW"].notna()]
        cdw_df = cdw_df[(phase_limits[0] <= cdw_df["TAI"]) & (cdw_df["TAI"] <= phase_limits[1])]

        # Standard deviation
        std_values, absolute_sigma = self.get_offline_std("CDW", phase_limits, cdw_df)

        popt, pcov = curve_fit(lin_fun, cdw_df["TAI"], cdw_df["CDW"], sigma = std_values, absolute_sigma = absolute_sigma)
        time_arr = np.array(cdw_df["TAI"])
        fit = lin_fun(time_arr, *popt)

        self.linear_fit_CDW[phase_name] = {"time": time_arr, "fit": fit, "popt": popt, 
                                           "pcov": pcov}
        print("{2}: fitted CDW: ax + b: a:{0:.2f}, b:{1:.2f}".format(*popt, phase_name))

    def plot_CDW_fit(self):
        fig, ax = plt.subplots(1)

        cdw_df = self.offline_df[self.offline_df["CDW"].notna()]
        std_values, has_std = self.get_offline_std("CDW", mean_df = cdw_df)

        # Plot measurements with or without standard deviation
        if has_std:
            ax.errorbar(cdw_df["TAI"], cdw_df["CDW"], yerr = std_values, fmt = "k.", ms = 10, label = "CDW measurements",
                        barsabove = True, ecolor = "k", elinewidth = 1)
        else:
            ax.plot(cdw_df["TAI"], cdw_df["CDW"], '.', ms = 12, label = "CDW measurements")
        

        # Plot fit
        for phase_name, dic in self.linear_fit_CDW.items():
            ax.plot(dic["time"], dic["fit"], lw = 2, label = "Fit {0}".format(phase_name))
        
        proteome_times = []
        at_proteome_times = []
        for t, x in self.estimated_cdw.items():
            proteome_times.append(t)
            at_proteome_times.append(x)

        ax.plot(proteome_times, at_proteome_times, "rx", label = "Proteome timepoints")
        ax.legend()
        ax.set_xlabel("Time after inoculation [h]")
        ax.set_ylabel("CDW [g/L]")
        ax.set_title(self.strain_name)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.show()

    def fit_linear_CDW_and_predict_growth_rate_for_phases(self, timepoint_phasename_list):
        for (timepoints, phase_name) in timepoint_phasename_list:
            self.fit_linear_CDW(phase_name)
            self.predict_growth_rate(timepoints, phase_name)

    def predict_growth_rate(self, timepoints, phase_name):
        """
        dX/dt = u*X
        u = (1/X) * dX/dt
        """
        dic = self.linear_fit_CDW[phase_name]
        X = lin_fun(np.array(timepoints), *dic["popt"])
        for t, x in zip(timepoints, X):
            self.growth_rates[t] = dic["popt"][0] / x
            self.estimated_cdw[t] = x

    def fit_substrates_for_phases(self, substrate_list, phase_list):
        for substrate in substrate_list:
            for phase in phase_list:
                self.fit_substrate(substrate, phase)

    def fit_substrate(self, substrate_name, phase_name = None, phase_limits = None):
        if phase_name:
            phase_limits = self.phases[phase_name]
        elif phase_limits:
            phase_name = "{0}-{1}".format(*phase_limits)
        else:
            print("No time limits given")


        # Get timepoints between the limits
        df = self.offline_df[(phase_limits[0] <= self.offline_df["TAI"]) & (self.offline_df["TAI"] <= phase_limits[1])]
        df = df[df[substrate_name].notna()]

        # Standard deviation
        std_values, absolute_sigma = self.get_offline_std(substrate_name, phase_limits, df)


        popt, pcov = curve_fit(lin_fun, df["TAI"], df[substrate_name], sigma = std_values, absolute_sigma = absolute_sigma)

        time_arr = np.arange(phase_limits[0], phase_limits[1]+1, 1)
        fit = lin_fun(time_arr, *popt)

        phases_dict = self.substrate_fits[substrate_name]
        phases_dict[phase_name] = {"time": time_arr, "fit": fit, "popt": popt, "type": "linear", "pcov": pcov}
        # print(substrate_name, popt, df[substrate_name])

    def predict_uptake_rates_for_phases(self, substrate_list,  timepoint_phasename_list):
        for substrate in substrate_list:
            for (timepoints, phase_name) in timepoint_phasename_list:
                self.predict_uptake_rates(substrate, phase_name, timepoints)


    def predict_uptake_rates(self, substrate, phase_name, timepoints, cdw_phase_name = None):
        # print(self.substrate_fits)
        phases_dict = self.substrate_fits[substrate]
        dic = phases_dict[phase_name]
        if dic["type"] == "linear":
            ds_dt = np.ones(len(timepoints)) * dic["popt"][0]
             # Save timepoints for rate estimation
            dic["rate_times"] = timepoints
            dic["at_rate_times"] = lin_fun(np.array(timepoints), *dic["popt"])
        
        elif dic["type"] == "sigmoid":
            dt = dic["time"][1] - dic["time"][0]
            ds_dt = np.zeros(len(timepoints))
            for i, t in enumerate(timepoints):
                ds_dt[i] = derivative(fsigmoid, x0 = t, dx = dt, args = (dic["popt"]))
            # ds_dt = np.gradient(dic["fit"], dt)
            # Save timepoints for rate estimation
            dic["rate_times"] = timepoints
            dic["at_rate_times"] = fsigmoid(np.array(timepoints), *dic["popt"])

        else: 
            print("Fit not implemented")
            return None

        if cdw_phase_name:
            cdw_dic = self.linear_fit_CDW[cdw_phase_name]
        else:
            cdw_dic = self.linear_fit_CDW[phase_name]
        
        cdw_values = lin_fun(np.array(timepoints), *cdw_dic["popt"])
        
        rate_dic = self.substrate_rates[substrate]
        for i, (t, x) in enumerate(zip(timepoints, cdw_values)):
            rate_dic[t] = 1e3 * ds_dt[i] / x / MOLAR_MASS_DICT[substrate]
            # print("# ", substrate,  ds_dt, x)
        # print(self.substrate_rates)

    def plot_RED(self):
        substrate = "Undecylprodigiosin 2"
        phases_dict = self.substrate_fits[substrate]
        
        fig, ax = plt.subplots(1)
        ax2 = ax.twinx()
        ax2._get_lines.prop_cycler = ax._get_lines.prop_cycler

        df = self.offline_df[self.offline_df[substrate].notna()]
        
        # Get standard deviation values
        std_values, has_std = self.get_offline_std(substrate, phase_limits = None, mean_df = df)

        # Plot measurements with or without standard deviation
        if has_std:
            ax.errorbar(df["TAI"], df[substrate], yerr = std_values, fmt = "k.", ms = 10, label = "{0} measurements".format(substrate),
                        barsabove = True, ecolor = "k", elinewidth = 1)
        
        labels = []
        l = ax.plot(df["TAI"], df[substrate], 'k.', ms = 10, label = "{0} measurements".format(substrate))        
        labels += l

        proteome_times = []
        proteome_rates = []
        for phase_name, dic in phases_dict.items():
            # print(phase_name, dic)
            l = ax.plot(dic["time"], dic["fit"], label = "Fit {0} {1}".format(substrate, phase_name))
            proteome_times += list(dic["rate_times"])
            proteome_rates += list(dic["at_rate_times"])
            labels += l

        l = ax.plot(proteome_times, proteome_rates, "r+", label = "Proteome timepoints")
        labels += l

        l = ax2.plot(df["TAI"], RED_M145_OD, marker = ".", alpha = 0.8, lw = 0, ms = 10, label = "RED OD measurements")
        labels += l
        ax2.set_ylabel("OD [spec. units]")

        ax.legend(labels, [l.get_label() for l in labels])
        ax.set_xlabel("Time after inoculation [h]")
        ax.set_ylabel("{0} [g/L]".format(substrate))
        # ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_title("{0} {1} concentration".format(self.strain_name, substrate))
        plt.show()


    def plot_uptake_fit(self, substrate):
        phases_dict = self.substrate_fits[substrate]

        fig, ax = plt.subplots(1)
        
        df = self.offline_df[self.offline_df[substrate].notna()]
        
        # Get standard deviation values
        std_values, has_std = self.get_offline_std(substrate, phase_limits = None, mean_df = df)

        # Plot measurements with or without standard deviation
        print(len(df[substrate]), len(std_values))
        if has_std:
            ax.errorbar(df["TAI"], df[substrate], yerr = std_values, fmt = "k.", ms = 10, label = "{0} measurements".format(substrate),
                        barsabove = True, ecolor = "k", elinewidth = 1)
        else:
            ax.plot(df["TAI"], df[substrate], 'k.', ms = 12, label = "{0} measurements".format(substrate))        

        proteome_times = []
        proteome_rates = []
        for phase_name, dic in phases_dict.items():
            # print(phase_name, dic)
            ax.plot(dic["time"], dic["fit"], label = "Fit {0} {1}".format(substrate, phase_name))
            proteome_times += list(dic["rate_times"])
            proteome_rates += list(dic["at_rate_times"])
        ax.plot(proteome_times, proteome_rates, "r+", label = "Proteome timepoints")
        ax.legend()
        ax.set_xlabel("Time after inoculation [h]")
        ax.set_ylabel("{0} [g/L]".format(substrate))
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_title("{0} {1} concentration".format(self.strain_name, substrate))
        plt.show()

    def rates_as_df(self, save_name = None, column_order = None, float_format = "%.3f"):
        index = sorted(list(self.growth_rates.keys()))
        self.rate_df = pd.DataFrame(index = index)
        self.rate_df["Growth rate"] = pd.Series(self.growth_rates)
        self.rate_df["Estimated CDW"] = pd.Series(self.estimated_cdw)
        CO2 = self.get_CO2_at_proteome_timepoints()
        self.rate_df["CO2"] = CO2
        print(self.rate_df)
        if self.rate_df["Growth rate"].max() > self.max_growth_rate:
            print("Values above max rate:\n", self.rate_df["Growth rate"] > self.max_growth_rate)
            
            self.rate_df.loc[self.rate_df["Growth rate"] > self.max_growth_rate, "Growth rate"] = self.max_growth_rate
            print("Replaced by ", self.max_growth_rate)

        # uptake rates
        for substrate, sub_dict in self.substrate_rates.items():
            self.rate_df[substrate] = pd.Series(sub_dict)

        if save_name:
            self.rate_df.to_csv(save_name, sep = ",", index_label = "TAI", columns = column_order, float_format = float_format)

    def fit_sigmoidial_to_substrate(self, substrate_name, phase_name, n = 3, bounds = None):
        phase_limits = self.phases[phase_name]
        df = self.offline_df[(phase_limits[0] <= self.offline_df["TAI"]) & (self.offline_df["TAI"] <= phase_limits[1])]
        df = df[df[substrate_name].notna()]

        popt, _ = curve_fit(fsigmoid, df["TAI"], df[substrate_name], bounds = bounds)

        time_arr = np.arange(df["TAI"].min(), df["TAI"].max(), 0.1)
        fit = fsigmoid(time_arr, *popt)

        phases_dict = self.substrate_fits[substrate_name]
        phases_dict[phase_name] = {"time": time_arr, "fit": fit, "popt": popt, "type": "sigmoid"}

        # if 0:
        #     plt.plot(df["TAI"], df[substrate_name], ".")
        #     plt.plot(time, fit)
        #     plt.show()
    def plot_CDW_and_CO2(self):
        fig, ax = plt.subplots(1)
        l2 = ax.plot(self.online_df["TAI"], self.online_df["CO2"], label = "CO2 measurements")
        ax2 = ax.twinx()
        ax2._get_lines.prop_cycler = ax._get_lines.prop_cycler


        cdw_df = self.offline_df[self.offline_df["CDW"].notna()]
        std_values, has_std = self.get_offline_std("CDW", mean_df = cdw_df)
        # has_std = False
        # Plot measurements with or without standard deviation
        if has_std:
            err_container = ax2.errorbar(cdw_df["TAI"], cdw_df["CDW"], yerr = std_values, color = "orange", fmt = ".", ms = 1, label = "CDW measurements",
                        barsabove = True, elinewidth = 1)
        l1 = ax2.plot(cdw_df["TAI"], cdw_df["CDW"], color = "orange", marker = ".", lw = 0, ms = 12, label = "CDW measurements")#c = "#b4d8ee"
        

        labels = l1 + l2
        ax.legend(labels, [l.get_label() for l in labels], loc = 0)
        ax.set_xlabel("Time after inoculation [h]")
        ax.set_ylabel("CO2 [mmol/L/h]")
        ax2.set_ylabel("CDW [g/L]")
        # ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_title("CDW and CO2 measurements {0}".format(self.strain_name))
        plt.show()


    def get_CO2_at_proteome_timepoints(self):

        CO2_values = []
        for i, (t, row) in enumerate(self.rate_df.iterrows()):
            idx = (self.online_df["TAI"] - t).abs().idxmin()
            CO2_values.append(self.online_df.loc[idx, "CO2"] / row["Estimated CDW"])
        return CO2_values
    
    def model_rates(self, model_fn, carbon, save_fn = None):
        model = read_sbml_model(model_fn)
        
        new_df = pd.DataFrame(columns = ["Estimated growth rate", "FBA growth rate", "Carbon balance"], dtype="float64")

        for i, (idx, row) in enumerate(self.rate_df.iterrows()):
            new_df.loc[idx, "Estimated growth rate"] = row["Growth rate"]
            
            tmp_balance = []
            with model:
                for c in carbon:
                    if np.isnan(row[c]):
                        tmp_balance.append(0)
                    else:
                        # Carbon balance
                        tmp_balance.append(row[c]*CARBON_DICT[c])
                        
                        # FBA growth
                        model.reactions.get_by_id(MODEL_REACTION_DICT[c]).lower_bound = row[c]
                s = model.optimize()
                if s.status == "infeasible":
                    new_df.loc[idx, "FBA growth rate"] = np.nan
                else:
                    new_df.loc[idx, "FBA growth rate"] = s.objective_value

            new_df.loc[idx, "Carbon balance"] = -sum(tmp_balance)
        
        print(new_df.dtypes)
        if save_fn:
            new_df.to_csv(str(save_fn), sep = ",", index_label = "TAI", float_format = "%.4f")
        return new_df





def replace_str_by_0(column, replace_by = 0):
    l = []
    for x in column:
        if pd.isna(x):
            y = x
        elif isinstance(x, float):
            y = x
        else:
            try:
                y = float(x.replace(",", "."))
            except:
                y = replace_by
        l.append(y)
    return l

def read_online(fn):
    df = pd.read_csv(fn,  sep = ";", header = 0, skiprows = [0,1,2,3,5], decimal = ".")
    return df


def read_offline(fn):
    if fn is None:
        return None

    df = pd.read_csv(fn,  sep = ",", header = 0, skiprows = [0,1,2,3,5], decimal = ".")
    # Replace " < 0.5 " by 0
    df["Undecylprodigiosin 2"] = replace_str_by_0(df["Undecylprodigiosin 2"])
    df["Germicidin-A"] = replace_str_by_0(df["Germicidin-A"])
    df["Germicidin-B"] = replace_str_by_0(df["Germicidin-B"])
    

    df[["Undecylprodigiosin 2", "Germicidin-A", "Germicidin-B"]] *= 1e-6
    return df

def piecewise_linear_3(x, x0, x1, y0, a1, a2, a3):
    condlist = [x <= x0, x <= x1, x1 < x]
    funclist = [lambda x: a1*x + y0, lambda x: a2*x + (a1*x0 + y0 - a2*x0), lambda x: a3*x + (a2*(x1-x0) + a1*x0 + y0) - a3*x1, lambda x: 0]
    print(len(condlist), len(funclist))
    return np.piecewise(x, condlist, funclist)

def piecewise_linear_2(x, x0, y0, a1, a2):
    condlist = [x <= x0]
    funclist = [lambda x: a1*x + y0 - a1*x0, lambda x: a2*x + y0 - a2*x0]
    return np.piecewise(x, condlist, funclist)

def fsigmoid(x, a, b, c):
    return c *  1.0 / (1.0 + np.exp(-a*(x-b)))

def lin_fun(x, a, b):
    return a*x + b

def exp_fun(x, a, b):
    return a*np.exp(b*x)


if __name__ == '__main__':
    column_order = ["Estimated CDW", "Growth rate", "Glucose", "Glutamic acid", "Undecylprodigiosin 2", "Germicidin-A", "Germicidin-B", "CO2"]
    if 1:
        # M145
        folder = "Average"
        online_M145_fn = GROWTH_DATA_FOLDER / "M145" / folder / "M145_online_data.csv"
        offline_M145_fn = GROWTH_DATA_FOLDER / "M145" / folder / "M145_offline_data.csv"
        offline_M145_std_fn = GROWTH_DATA_FOLDER / "M145" / folder / "M145_offline_data_std.csv"
        proteomic_timepoints_M145 = [21, 29, 33, 37, 41, 45, 49, 53, 57]
        p1 = proteomic_timepoints_M145[:3] 
        p2 = proteomic_timepoints_M145[3:6]
        p3 = proteomic_timepoints_M145[6:] 

        # RE = RateEstimator("M145", online_M145_fn, offline_M145_fn, proteome_timepoints = proteomic_timepoints_M145)
        RE = RateEstimator("M145", online_M145_fn, offline_M145_fn, offline_M145_std_fn, proteomic_timepoints_M145)

        RE.set_phase("Exponential phase", 14, 21)
        RE.set_phase("Linear phase 1", 21, 37)
        RE.set_phase("Linear phase 2", 37, 49)
        RE.set_phase("Linear phase 3", 49, 66)
        RE.set_phase("Red linear phase", 41, 54)
        RE.set_phase("Germicidin phase", 38, 66)


        RE.fit_exponential_phase_to_CO2("Exponential phase")
        # RE.plot_exponential_fit_CO2()
        # RE.plot_CDW_and_CO2()
    
        
        RE.fit_linear_CDW_and_predict_growth_rate_for_phases([(p1,  "Linear phase 1"),
                                                              (p2,  "Linear phase 2"), 
                                                              (p3,  "Linear phase 3")])
        # RE.plot_CDW_fit()

        RE.fit_substrates_for_phases(["Glucose", "Glutamic acid"], ["Linear phase 1", "Linear phase 2", "Linear phase 3"])
        RE.predict_uptake_rates_for_phases(["Glucose", "Glutamic acid"], [(p1, "Linear phase 1"), (p2, "Linear phase 2"), (p3, "Linear phase 3")])

        # RE.plot_uptake_fit("Glucose")
        # RE.plot_uptake_fit("Glutamic acid")

        # RE.fit_sigmoidial_to_substrate("Undecylprodigiosin 2", "Full time serie", 3, bounds = ([0, 20, 2], [1, 60, 3]))
        RE.fit_substrate("Undecylprodigiosin 2", "Red linear phase")
        RE.predict_uptake_rates("Undecylprodigiosin 2", "Red linear phase", p2, "Linear phase 2")
        RE.predict_uptake_rates("Undecylprodigiosin 2", "Red linear phase", p3, "Linear phase 3")

        
        RE.fit_substrates_for_phases(["Germicidin-A", "Germicidin-B"], ["Germicidin phase"])
        RE.predict_uptake_rates("Germicidin-A", "Germicidin phase", p2, "Linear phase 2")
        RE.predict_uptake_rates("Germicidin-A", "Germicidin phase", p3, "Linear phase 3")
        RE.predict_uptake_rates("Germicidin-B", "Germicidin phase", p2, "Linear phase 2")
        RE.predict_uptake_rates("Germicidin-B", "Germicidin phase", p3, "Linear phase 3")
        
        # RE.plot_uptake_fit("Undecylprodigiosin 2")
        # RE.plot_RED()
        RE.rates_as_df(save_name = GROWTH_DATA_FOLDER / "M145" / folder /  "M145_estimated_rates.csv", column_order = column_order, float_format = "%.6f")
        
        carbon = ["Glutamic acid", "Glucose", "Germicidin-A", "Germicidin-B", "CO2", "Undecylprodigiosin 2"]
        RE.model_rates(str(scoGEM_FN), carbon, save_fn =  GROWTH_DATA_FOLDER / "M145" / folder /  "M145_carbon_balance.csv")
    if 1:

        # M1152
        folder = "Average"
        online_M1152_fn = GROWTH_DATA_FOLDER / "M1152" / folder / "M1152_online_data.csv"
        offline_M1152_fn = GROWTH_DATA_FOLDER / "M1152" / folder / "M1152_offline_data.csv"
        offline_M1152_std_fn = GROWTH_DATA_FOLDER / "M1152" / folder / "M1152_offline_data_std.csv"


        proteomic_timepoints_M1152 = [33, 41, 45, 49, 53, 57, 61, 65]


        RE = RateEstimator("M1152", online_M1152_fn, offline_M1152_fn, offline_M1152_std_fn, proteomic_timepoints_M1152)
        RE.set_phase("Exponential phase", 10, 29)
        RE.set_phase("Linear phase 1", 29, 50)
        RE.set_phase("Linear phase 1b", 27, 41)
        RE.set_phase("Linear phase 2", 49, 54)
        RE.set_phase("Linear phase 3", 54, 75)
        RE.set_phase("Linear phase 2b", 41, 57)
        RE.set_phase("Germicidin phase", 40, 66)


        p1 = proteomic_timepoints_M1152[:4]
        p1b = proteomic_timepoints_M1152[:2]
        p2 = proteomic_timepoints_M1152[4:5]
        p2b = proteomic_timepoints_M1152[2:5]
        p3 = proteomic_timepoints_M1152[5:]

        RE.fit_exponential_phase_to_CO2("Exponential phase")
        # RE.plot_exponential_fit_CO2()
        # RE.plot_CDW_and_CO2()

        
        RE.fit_linear_CDW_and_predict_growth_rate_for_phases([(p1,  "Linear phase 1"),
                                                              (p2,  "Linear phase 2"), 
                                                              (p3,  "Linear phase 3")])
        # RE.plot_CDW_fit()


        # RE.fit_substrates_for_phases(["Glucose", "Glutamic acid"], ["Linear phase 1", "Linear phase 2b", "Linear phase 3"])
        RE.fit_substrates_for_phases(["Glucose", "Glutamic acid"], ["Linear phase 1b", "Linear phase 2b", "Linear phase 3"])
        
        RE.predict_uptake_rates("Glucose", "Linear phase 1b", p1b, cdw_phase_name = "Linear phase 1")
        RE.predict_uptake_rates("Glutamic acid", "Linear phase 1b", p1b, cdw_phase_name = "Linear phase 1")

        RE.predict_uptake_rates("Glucose", "Linear phase 2b", p2b, cdw_phase_name = "Linear phase 2")
        RE.predict_uptake_rates("Glutamic acid", "Linear phase 2b", p2b, cdw_phase_name = "Linear phase 2")

        RE.predict_uptake_rates_for_phases(["Glucose", "Glutamic acid"], [(p3, "Linear phase 3")])
        # RE.predict_uptake_rates_for_phases(["Glucose", "Glutamic acid"], [(p1, "Linear phase 1b"), (p3, "Linear phase 3")])

        # RE.plot_uptake_fit("Glucose")
        # RE.plot_uptake_fit("Glutamic acid")


        RE.fit_substrates_for_phases(["Germicidin-A", "Germicidin-B"], ["Germicidin phase"])
        RE.predict_uptake_rates("Germicidin-A", "Germicidin phase", p2, "Linear phase 2")
        RE.predict_uptake_rates("Germicidin-A", "Germicidin phase", p3, "Linear phase 3")
        RE.predict_uptake_rates("Germicidin-B", "Germicidin phase", p2, "Linear phase 2")
        RE.predict_uptake_rates("Germicidin-B", "Germicidin phase", p3, "Linear phase 3")
        
        # RE.plot_uptake_fit("Germicidin-A")
        # RE.plot_uptake_fit("Germicidin-B")

        

        RE.rates_as_df(save_name = GROWTH_DATA_FOLDER / "M1152" / folder / "M1152_estimated_rates.csv", column_order = column_order, float_format = "%.7f")
        carbon = ["Glutamic acid", "Glucose", "Germicidin-A", "Germicidin-B", "CO2"]
        RE.model_rates(str(scoGEM_FN), carbon, save_fn = GROWTH_DATA_FOLDER / "M1152" / folder / "M1152_carbon_balance.csv")


