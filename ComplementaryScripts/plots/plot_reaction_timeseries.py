# -*- coding: utf-8 -*-
"""
Author: Snorre Sulheim
Date: 09.04.2019

# Description
Create plots of how the flux changes during the timecourse and the corresponding genes

## Plot 1
Plot the sampled flux for a given reaction for all timesteps for both strains
function: plot_all_flux()

##



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

# Genes in pho-regulon
DELAYED_GENES = ["SCO1196", "SCO1968", "SCO2286", "SCO1906"]
NOT_RESPONDING_GENES = ["SCO1565", "SCO4229"]
EXPORTED_GDPD_GENES = ["SCO1565", "SCO1968", "SCO7550"]
CYTOSOLIC_GDPD_GENES = ["SCO1090", "SCO1419", "SCO3976", "SCO5661"]
TAT_SECRETED_GENES = ["SCO0736", "SCO1172", "SCO1196", "SCO1356", "SCO1432", "SCO1565", "SCO1590", "SCO1639", "SCO1906", "SCO2068", "SCO2286", "SCO2758", "SCO2780", "SCO2786", "SCO3484", "SCO3790", "SCO4672", "SCO6052", "SCO6198", "SCO6272", "SCO6281", "SCO6457", "SCO6580", "SCO6594", "SCO6691", "SCO7631", "SCO7677"]
CYO2a_GENES = ["SCO2150", "SCO2149", "SCO7236", "SCO2148", "SCO7120"]
ATP_SYNTHASE_GENES = ["SCO5368", "SCO5370", "SCO5373", "SCO5374"]
INOSITOL_DEHYDROGENASE_GENES = ["SCO6255", "SCO6984", "SCO7254"]

def convert_sample_number_to_timepoint(row):
    if row["Strain"] == "M145":
        timepoint_dict = {"P6":1, "P14":2, "P18":3, "P22":4, "P26":5, "P30":6, "P34":7, "P38":8, "P42":9}
    else:
        timepoint_dict = {"P18":1, "P26":2, "P30":3, "P34":4, "P38":5, "P42":6, "P46":7, "P50":8}
    return timepoint_dict[row["Sample number"]]

def convert_to_timepoint2(row):
    if row["Strain"] == "M145":
        timepoint_dict = {21:1, 29:2, 33:3, 37:4, 41:5, 45:6, 49:7, 53:8, 57:9}
    else:
        timepoint_dict = {33:1, 41:2, 45:3, 49:4, 53:5, 57:6, 61:7, 65:8, 69:9}
    return timepoint_dict[int(row["Hours"])]

def convert_to_timepoint(strain, hour):
    if strain == "M145":
        timepoint_dict = {21:1, 29:2, 33:3, 37:4, 41:5, 45:6, 49:7, 53:8, 57:9}
    else:
        timepoint_dict = {33:1, 41:2, 45:3, 49:4, 53:5, 57:6, 61:7, 65:8, 69:9}
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
            new_df["Time point"] = convert_to_timepoint(strain, hours)
            df_list.append(new_df)
        r_df = pd.concat(df_list, axis = 0)
        fn = folder / "reactions"/ "{0}.csv".format(r)
        r_df.to_csv(fn, index = False)

def plot_all_flux(folder, reaction):
    df = pd.read_csv(Path(folder) / "reactions" / "{0}.csv".format(reaction))
    print(df)
    sns.violinplot(x = "Time point", y = reaction, hue = "Strain", data = df, inner = "stick", split = False, cut = 0, scale = "count")
    # sns.stripplot(x = "Time point", y = reaction, hue = "Strain", data = df, inner = None, split = False)
    # sns.swarmplot(x = "Time point", y = reaction, hue = "Strain", data = df)
    # sns.boxplot(x = "Time point", y = reaction, hue = "Strain", data = df)

    plt.show()

def get_gene_data(gene_expression_csv, gene_ids):
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

    df_selected_genes = df_gene.loc[:, gene_ids+["Strain", "Hours"]]
    df_selected_genes["Time point"] = df_selected_genes.apply(convert_to_timepoint2, axis = 1)

    return df_selected_genes

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

def plot_genes_transcriptome_proteome(proteome_norm_tsv, gene_expression_csv, gene_ids):
    """
    Plot the timecourse for a list of genes, one gene ni each plot for both
    transcriptome and proteome data
    """
    df_transcr = get_gene_data(gene_expression_csv, gene_ids)
    df_prot = get_proteome(proteome_norm_tsv, gene_ids)


    for i, gene in enumerate(gene_ids):
        # ax = axes[i]
        fig, ax = plt.subplots(1, figsize  = (8, 6))
        ax = sns.lineplot(x = "Time point", y = gene, hue = "Strain", data = df_transcr, legend = False)#, ax = ax)
        ax2 = ax.twinx()
        ax2 = sns.lineplot(x = "Time point", y = gene, hue = "Strain", data = df_prot, err_style = "bars", ax = ax2,
                           legend = False, alpha = 0.8)
        for i in range(len(ax2.lines)):
            # print(ax2.lines[i].__dict__)
            # if ax2.lines[i]._linestyle is not None:
            ax2.lines[i].set_linestyle("--")
        lines = ax.lines + ax2.lines[::2]
        labels = ["M145 transcriptome", "M1152 transcriptome", "M145 proteome", "M1152 proteome"]
        ax.legend(lines, labels, bbox_to_anchor=(1,0), loc="lower right", 
                bbox_transform=fig.transFigure, ncol=1)

        ax.set_ylabel("Normalized transcriptome data")
        ax2.set_ylabel("Normalized proteome data")
        ax.set_title(gene)
        plt.subplots_adjust(right = 0.85, bottom = 0.2)
        plt.show()

def plot_genes(gene_expression_csv, gene_ids):
    df_gene = get_gene_data(gene_expression_csv, gene_ids)
    
    df_melt = df_selected_genes.melt(["Time point", "Strain"], value_vars = gene_ids, var_name = "Gene", value_name = "Normalized count")
    print(df_melt)
    sns.lineplot(x = "Time point", y = "Normalized count", hue = "Gene", style = "Strain", data = df_melt, err_style = "bars")
    plt.show()




def get_proteome(proteome_norm_tsv, gene_ids):
    df = pd.read_csv(proteome_norm_tsv, sep = "\t", index_col = 0).T
    df_selected_genes = df.loc[:, gene_ids] # Select only a few genes
    df_selected_genes.index = df_selected_genes.index.str.split("_", expand = True) # Split index into strain, f-number and sample time
    df_selected_genes.index.set_names(["Strain", "F number", "Sample number"], inplace = True) # Set index names
    df_selected_genes.reset_index(inplace = True) # Convert multiindex frame to df with three extra colums
    print(df_selected_genes)
    df_selected_genes["Time point"] = df_selected_genes.apply(convert_sample_number_to_timepoint, axis = 1)
    return df_selected_genes

def plot_proteome(proteome_norm_tsv, gene_ids):
    df_prot = get_proteome(proteome_norm_tsv, gene_ids)
    sns.lineplot(data = df_prot)
    plt.show()
    df_melt = df_prot.melt(["Time point", "Strain"], value_vars = gene_ids, var_name = "Gene", value_name = "Normalized counts")
    sns.lineplot(x = "Time point", y = "Normalized counts", hue = "Gene", style = "Strain", data = df_melt, err_style = "bars")
    plt.show()

if __name__ == '__main__':
    mean_flux_co2_scaled = "C:/Users/snorres/Google Drive/scoGEM community model/Supporting information/Model/ec-RandSampComb_AllSamples_co2.csv"
    all_random_samples_folder = "C:/Users/snorres/OneDrive - SINTEF/SINTEF projects/INBioPharm/scoGEM/random sampling/Eduard random sampling"
    gene_expression_csv = "C:/Users/snorres/OneDrive - SINTEF/SINTEF projects/INBioPharm/scoGEM/manuscript/m145-m1152-normalized_allfermenters.tsv"
    proteome_norm_tsv = "C:/Users/snorres/Google Drive/scoGEM community model/Supporting information/Correlation_ProtRNA_Ed/proteomeNorm.tsv"



    # plot_gene_and_single_flux(mean_flux_co2_scaled, gene_expression_csv)
    # plot_reactions(mean_flux_co2_scaled, reaction_ids = DIFFERENT_M1152)

    if 0:
        # Plot no 1
        get_all_flux(all_random_samples_folder, DIFFERENT_M1152+DELAYED_M1152)
        for r in DIFFERENT_M1152+DELAYED_M1152:
            plot_all_flux(all_random_samples_folder, r)

    if 0:
        # plot no 2
        gene_ids = ["SCO1196", "SCO1968", "SCO2286"]
        plot_genes(gene_expression_csv, gene_ids)

    if 0:
        #  Plot proteome (no 3)
        gene_ids = ["SCO1196", "SCO1968", "SCO2286"]
        plot_proteome(proteome_norm_tsv, gene_ids)

    if 0:
        # Plot proteome and transcriptome in same plot
        
        # plot_delayed_genes_transcriptome_proteome(proteome_norm_tsv, gene_expression_csv, CYO2a_GENES)
        # plot_delayed_genes_transcriptome_proteome(proteome_norm_tsv, gene_expression_csv, ATP_SYNTHASE_GENES)
        # plot_delayed_genes_transcriptome_proteome(proteome_norm_tsv, gene_expression_csv, DELAYED_GENES)
        plot_genes_transcriptome_proteome(proteome_norm_tsv, gene_expression_csv, INOSITOL_DEHYDROGENASE_GENES)
    if 1:
        get_all_flux(all_random_samples_folder, ["INS2D"])
        plot_all_flux(all_random_samples_folder, "INS2D")
        