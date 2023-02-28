#!/usr/bin/python3

import os
import sys
import numpy as np
from pprint import pprint
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt


#############################################################
### See the following web sites for plt.rcParams          ###
### https://qiita.com/qsnsr123/items/325d21621cfe9e553c17 ###
### https://beiznotes.org/matplot-cmap-list/              ###
#############################################################
plt.rcParams["font.family"]       = "sans-serif"
plt.rcParams["font.size"]         = 14
plt.rcParams["xtick.direction"]   = "in"
plt.rcParams["ytick.direction"]   = "in"
plt.rcParams["xtick.major.width"] = 1.0
plt.rcParams["ytick.major.width"] = 1.0
plt.rcParams["axes.linewidth"]    = 1.0
plt.rcParams["image.cmap"]        = "viridis"


# hold data to be plotted
class Data:
    def __init__(self, names):
        for indx in range(len(names)):
            self.__dict__[names[indx]] = None

    def set_values(self, name, dist, hf_ene, casci_ene, vqe_ene):
        self.__dict__[name] = (dist, hf_ene, casci_ene, vqe_ene)

    def get_plot_values1(self, name):
        return self.__dict__[name][0], self.__dict__[name][1]  # dist vs. hf_ene

    def get_plot_values2(self, name):
        return self.__dict__[name][0], self.__dict__[name][2]  # dist vs. casci_ene

    def get_plot_values3(self, name):
        return self.__dict__[name][0], self.__dict__[name][3]  # dist vs. vqe_ene


file_names = sys.argv[1:]  # Exclude ".py" file
data = Data(file_names)


np.set_printoptions(formatter={"float": "{:.6e}".format})
fig, axes = plt.subplots(1, 1, tight_layout=True, figsize=(8, 6))  # (horizon, vertical)


for indx, file_name in enumerate(file_names):
    print(f"{file_name} is now read...")
    dist      = []
    hf_ene    = []
    casci_ene = []
    vqe_ene   = []

    with open(file_name) as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if i > 0:
                dist.append(float(line.split()[0]))
                hf_ene.append(float(line.split()[1]))
                casci_ene.append(float(line.split()[2]))
                vqe_ene.append(float(line.split()[3]))

    data.set_values(file_name, np.array(dist), np.array(hf_ene), np.array(casci_ene), np.array(vqe_ene))

    # plot data of file_name
    axes.plot(*data.get_plot_values1(file_name), label=file_name +" HF")
    axes.plot(*data.get_plot_values2(file_name), label=file_name +" CASCI")
    axes.plot(*data.get_plot_values3(file_name), label=file_name +" VQE")


ymin = None
ymax = None
axes.set_ylim(ymin, ymax)

xmin = None
xmax = None
axes.set_xlim(xmin, xmax)


#######################################################################################
### linestyle: "solid" == "-", "dashed" == "--", "dashdot" == "-.", "dotted" == ":" ###
### See the following web sites for grid setting                                    ###
### https://python.atelierkobato.com/linestyle/                                     ###
### https://www.pythoncharts.com/matplotlib/customizing-grid-matplotlib/            ###
#######################################################################################
axes.grid(visible=True, which="major", color="lightgrey", linestyle="-",  linewidth=0.8)
axes.grid(visible=True, which="minor", color="lightgrey", linestyle="-.", linewidth=0.8)


####################################################################################################################
### See https://www.pythoncharts.com/matplotlib/customizing-grid-matplotlib/ for minor grid (Not log scale ver.) ###
####################################################################################################################
axes.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
axes.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))


########################################################################################################################
### See https://www.delftstack.com/ja/howto/matplotlib/python-matplotlib-plot-superscript/ for superscript/subscript ###
########################################################################################################################
axes.set_xlabel(r"Distance (Angstrom)",  fontsize=16)
axes.set_ylabel(r"Energy (Eh)",          fontsize=16)

#################################################################################################
### See https://qiita.com/matsui-k20xx/items/291400ed56a39ed63462 for legend position setting ###
#################################################################################################
axes.legend(loc="best", borderaxespad=1, fontsize=12, edgecolor="black")
axes.set_title("PES for BeH2", fontsize=16)
plt.savefig("PES.jpg", dpi=300)
plt.show()
