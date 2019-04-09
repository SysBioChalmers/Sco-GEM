# -*- coding: utf-8 -*-
"""
Author: Snorre Sulheim
Date: 09.04.2019

# Description
Create plots of how the flux changes during the timecourse and the corresponding genes

"""

import cobra
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import seaborn as sns
import re


pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_columns', 1000)

M145_FERMENTORS = ["F516", "F517", "F518"]
M1152_FERMENTORS = ["F519", "F521", "F522"]

# Delayed reactions in M1152
DELAYED_M1152 = ["MCOATA", "TAGS140"]

# Different reactions in M1152
DIFFERENT_M1152 = ["PDH", "G6PDH2r", "PGL", "CYO2a", "ACONTa", "GLUDxi"]

def convert_to_timepoint(strain, hour):
    if strain == "M145":
        timepoint_dict = {21:1, 29:2, 33:3, 37:4, 41:5, 45:6, 49:7, 53:8, 57:9}
    else:
        timepoint_dict = {33:1, 41:2, 45:3, 49:4, 53:5, 57:6, 61:7, 65:8}
    return timepoint_dict[int(hour)]

def get_all_flux(folder, reactions):
    folder = Path(folder)
    df_dict = {}
    for i, file in enumerate(folder.glob("*.csv")):
        _, strain, hours = file.stem.split("_")
        print(strain, hours)
        df = pd.read_csv(folder / file, sep = "\t", header = None, index_col = 0)
        df = df.loc[df.index.str.contains("|".join(reactions)), :]
        df = revert_GECKO_reactions(df)
        # print(df.head())
        df_dict[(strain, hours)] = df.loc[reactions, :]
        
    for r in reactions:
        df_list = []
        for (strain, hours), df in df_dict.items():
            new_df = pd.DataFrame()
            new_df[r] = df.T[r]
            new_df["Strain"] = strain
            new_df["Hours"]  = hours
            new_df["Timepoint"] = convert_to_timepoint(strain, hours)
            df_list.append(new_df)
        r_df = pd.concat(df_list, axis = 0)
        fn = folder / "reactions"/ "{0}.csv".format(r)
        r_df.to_csv(fn, index = False)

def plot_all_flux(folder, reaction):
    df = pd.read_csv(Path(folder) / "reactions" / "{0}.csv".format(reaction))
    print(df)
    sns.violinplot(x = "Timepoint", y = reaction, hue = "Strain", data = df, inner = None, split = True)
    plt.show()

def get_gene_data(gene_expression_csv):
    # The csv file contains normalized rna-seq data for each of the three biological replicates for each strain
    df_gene_fermentors = pd.read_csv(gene_expression_csv, sep = "\t")
    df_gene_fermentors.set_index("Identifier", inplace = True)
    # print(df_gene_fermentors.T)

    df_gene = df_gene_fermentors.T

    # Set strain based on index
    df_gene.loc[df_gene.index.str.contains("|".join(M145_FERMENTORS)), "Strain"] = "M145"    
    df_gene.loc[df_gene.index.str.contains("|".join(M1152_FERMENTORS)), "Strain"] = "M1152"    
    

    # Set timepoint
    df_gene.loc[:, "Hours"] = df_gene.index.str[5:7]
    return df_gene

def get_single_flux(mean_flux_fn):
    df_flux = pd.read_csv(mean_flux_fn, index_col = 0)

    #Remove the gene annotation and name
    df_flux = df_flux.iloc[:, 2:].T
    # df_flux.loc[df_flux.index.str.contains("M145"), "Strain"] = "M145"
    # df_flux.loc[df_flux.index.str.contains("M1152"), "Strain"] = "M1152"
    # df_flux[["Strain", "Hours"]] = 
    df_flux.index = df_flux.index.str.split("_", expand = True)
    df_flux.index.set_names(["Strain", "Hours"], inplace = True)
    df_flux.reset_index(inplace = True)
    return df_flux 


def plot_gene_and_single_flux(mean_flux_fn, gene_expression_csv):
    df_flux = get_single_flux(mean_flux_fn)
    df_gene = get_gene_data(gene_expression_csv)
    
def plot_reactions(mean_flux_fn, reaction_ids):
    df_flux = get_single_flux(mean_flux_fn)
    # df_flux = df_flux.loc[:, reaction_ids]
    #
    print(df_flux.head())
    for r in reaction_ids:
        sns.lineplot(x = "Hours", y = r, hue = "Strain", data = df_flux)
        plt.show()



def revert_GECKO_reactions(full_flux_frame):
    """
    Sum all the different parallell reactions created by the GECKO-method
    - irreversible to reversible (REV)
    - arm_reactions?
    - _Nox reactions
    """
    all_reaction_ids = list(full_flux_frame.index)
    all_columns = list(full_flux_frame.columns)
    full_flux_frame.loc[:, "keep"] = False
    full_flux_frame.loc[:, "base id"] = list(full_flux_frame.index)
    # get arm reactions
    arm_ids = []
    non_arm_ids  = []

    for r_id in all_reaction_ids:
        if r_id[:4] == "arm_":
            full_flux_frame.loc[r_id, "keep"] = True
            if r_id[-4:] == "_REV":
                full_flux_frame.loc[r_id, "base id"] = r_id[4:-4]
                full_flux_frame.loc[r_id, all_columns] *= -1
            else:
                full_flux_frame.loc[r_id, "base id"] = r_id[4:]
            arm_ids.append(r_id[4:])
        else:
            non_arm_ids.append(r_id)


    for r_id in non_arm_ids:
        no_stripped_id = re.sub(r"No\d+\Z", "", r_id)
        if not no_stripped_id in arm_ids:
            full_flux_frame.loc[r_id, "keep"] = True
        
            if no_stripped_id[-4:] == "_REV":
                full_flux_frame.loc[r_id, "base id"] = no_stripped_id[:-4]
                full_flux_frame.loc[r_id, all_columns] *= -1
            else:
                full_flux_frame.loc[r_id, "base id"] = no_stripped_id
        
    full_flux_frame.set_index("base id", inplace = True)
    flux_frame = full_flux_frame.loc[full_flux_frame["keep"], all_columns]
    flux_frame = flux_frame.groupby(flux_frame.index, axis = 0).sum()
    # print(flux_frame)
    # print(flux_frame_sum)
    return flux_frame


if __name__ == '__main__':
    mean_flux_co2_scaled = "C:/Users/snorres/Google Drive/scoGEM community model/Supporting information/Model/ec-RandSampComb_AllSamples_co2.csv"
    all_random_samples_folder = "C:/Users/snorres/OneDrive - SINTEF/SINTEF projects/INBioPharm/scoGEM/random sampling/Eduard random sampling"
    gene_expression_csv = "C:/Users/snorres/OneDrive - SINTEF/SINTEF projects/INBioPharm/scoGEM/manuscript/m145-m1152-normalized_allfermenters.tsv"

    # plot_gene_and_single_flux(mean_flux_co2_scaled, gene_expression_csv)
    # plot_reactions(mean_flux_co2_scaled, reaction_ids = DIFFERENT_M1152)

    if 1:
        # get_all_flux(all_random_samples_folder, DIFFERENT_M1152+DELAYED_M1152)
        for r in DIFFERENT_M1152+DELAYED_M1152:
            plot_all_flux(all_random_samples_folder, r)