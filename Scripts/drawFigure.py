#!/usr/bin/python3
import csv
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Plot bar charts for different types of data
def plotBarChart(type, extension, showFig=True):
  results = pd.read_csv("results_" + type + ".csv")

  if "cmp_ms" in type or type == "stellar" or "cmp_cp" in type or type == "cmp_testing":
    selected = results
    hueStr = "Data Structure"
    # comparisons between TTK triangulation 
    if "cmp_cp" in type or "cmp_ms" in type:
      palette1 = sns.color_palette(muted1)
      hueOrder1 = ["Implicit", "Explicit", "TTK"]
      palette2 = sns.color_palette(muted2)
      hueOrder2 = ["TTK", "Explicit", "Implicit"]
  else:
    print("The type is not supported!\n")
    sys.exit(-1)
  # print(selected)
  
  sns.set(style="whitegrid", font="Times New Roman", font_scale=1.8)
  fig, sub = plt.subplots()
  if extension == "memory":
    if type == "cmp_estimated":
      ax = sns.barplot(ax=sub, x="Dataset", y="Estimated Storage", hue=hueStr, hue_order=hueOrder1, data=selected, palette=palette1, ci=None)
      h, l = sub.get_legend_handles_labels()
      sub.legend(h, ["TTK", "$k_L$ Stellar tree", "$k_S$ Stellar tree", "Implicit", "Explicit"])
      sub.set_ylabel("Estimated Storage (MB)")
    else:
      ax = sns.barplot(ax=sub, x="Dataset", y="Memory Usage", hue=hueStr, hue_order=hueOrder1, data=selected, palette=palette1, ci=None)
      sub.set_ylabel("Memory (MB)")
    plt.xticks(rotation=45)
    sub.set_xlabel("")
    if type == "cmp_testing":
      h, l = sub.get_legend_handles_labels()
      # Implicit, kL, kS, Explicit, TTK
      sub.legend(h, ["Implicit", "$k_L$ Stellar tree", "$k_S$ Stellar tree", "Explicit", "TTK"])
    sub.autoscale_view()
  elif extension == "time":
    ax = sns.barplot(ax=sub, x="Dataset", y="Total Time", hue=hueStr, hue_order=hueOrder2, data=selected, palette=palette2, ci=None)
    plt.xticks(rotation=45)
    sub.set_xlabel("")
    # sub.set_ylabel("Time (s)")
    # for visualizing speedup data
    if "speedup" in type:
      sub.set_ylabel("Speedup")
    if type == "cmp_topocluster":
      ax = sns.barplot(ax=sub, x="Dataset", y="Preprocessing", hue=hueStr, hue_order=hueOrder2, data=selected, palette=sns.color_palette(["#00ff00"]), ci=None)
      h, l = sub.get_legend_handles_labels()
      sub.legend(h, ["Explicit","Implicit", "Preprocessing"])
      sub.set_xlabel("")
      sub.set_ylabel("Time (s)")
    elif type == "cmp_estimated":
      h, l = sub.get_legend_handles_labels()
      sub.legend(h, ["TTK", "Explicit", "$k_S$ Stellar tree", "$k_L$ Stellar tree", "Implicit"])
    elif type == "cmp_testing":
      h, l = sub.get_legend_handles_labels()
      sub.legend(h, ["TTK", "Explicit", "$k_S$ Stellar tree", "$k_L$ Stellar tree", "Implicit"])
    sub.autoscale_view()

  sub.legend(fontsize=18)
  # For visualizing speedup data
  if "speedup" in type:
    ax.set(ylim=(0, 5))

  fig.set_figheight(6)
  fig.set_figwidth(8)
  fig.tight_layout(pad=0)

  if showFig:
    plt.show()
  
  fig.savefig(type + "_" + extension + ".png")

# Main entry of the program
if __name__ == "__main__":
  # folderNames = ["Explicit_Serial_Threshold", "Explicit_Serial_Combined/all", "Implicit_Serial_Threshold", "Implicit_Serial_Combined/all"]
  datasets = ["Redsea", "Engine", "Cat", "Sphere", "Foot", "Shapes", "Hole", "Stent"]
  # newDatasets = ["CAT", "SPHERE", "SHAPES", "HOLE"]
  thresholds = ['20', '100', '500', '1000', '5000']
  cacheSizes = ['5', '10', '20', '50', '100']
  colorScheme = ['#abdda4', '#fdae61', '#ffffbf', '#d7191c', '#2b83ba'] # TTK, kL, kS, Implicit, Explicit

  # color palettes
  muted1 = [colorScheme[3], colorScheme[4], colorScheme[0]] # Implicit; Explicit; TTK
  muted2 = [colorScheme[0], colorScheme[4], colorScheme[3]] # TTK, Explicit, Implicit

  # parameters for plots
  # plt.rcParams['text.usetex'] = True
  plt.rcParams["font.family"] = 'Times New Roman'
  plt.rcParams["font.size"] = 16

  plotBarChart("cmp_cp", "memory")
  plotBarChart("cmp_cp", "time")

  print("Finished!\n")
