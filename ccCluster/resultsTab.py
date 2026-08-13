from __future__ import print_function

__author__ = "Gianluca Santoni & Tiankun Zhou"
__copyright__ = "Copyright 2015-2026"
__credits__ = ["Gianluca Santoni, Alexander Popov"]
__license__ = ""
__version__ = "1.0"
__maintainer__ = "Tiankun Zhou"
__email__ = "tiankun.zhou@esrf.fr"
__status__ = "Beta"


from PyQt5 import QtGui, QtCore, QtWidgets
from PyQt5.QtWidgets import QWidget, QApplication
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from pathlib import Path
from scipy.cluster import hierarchy
from .clustering import extractXSCALEStat
import collections
import operator
import stat

# implement the default mpl key bindings
import sys
import os
import subprocess
import matplotlib.patches as mpatches

#create results tab to plot the SCALE.LP statistics
class SinglePlotTab(QtWidgets.QWidget):
    def __init__(self, ProcessedDir:str):
        super().__init__()

        self.ProcessedDir= ProcessedDir
        #setup widget
        self.tabLayout= QtWidgets.QVBoxLayout(self)
        self.Title=QtWidgets.QLabel(self)
        ResultName = Path(self.ProcessedDir).name.split('_')
        self.Title.setText(f"Threshold: {ResultName[2]} and Group number: {ResultName[3]}")
        self.Title.setStyleSheet("color: black; font-weight: bold; font-size: 14px;")

        #set up the plot
        self.statsPlot = Figure()
        self.Ax = self.statsPlot.add_subplot(111)
        self.statsCanvas = FigureCanvas(self.statsPlot)
        self.statsBar= NavigationToolbar(self.statsCanvas, self)

        #Buttons to plot stats
        self.buttonBar = QtWidgets.QWidget()
        self.barLayout= QtWidgets.QHBoxLayout(self.buttonBar)

        self.ccVsR= QtWidgets.QPushButton()
        self.ccVsR.setText("CC. vs Res")
        self.ccVsR.clicked.connect(lambda:self.plotStats(0, 4, "CC. vs Res"))

        self.compVsR= QtWidgets.QPushButton()
        self.compVsR.setText("comp vs Res")
        self.compVsR.clicked.connect(lambda:self.plotStats(0, 10, "comp vs Res"))

        self.RobsVsR= QtWidgets.QPushButton()
        self.RobsVsR.setText("Robs vs Res")
        self.RobsVsR.clicked.connect(lambda:self.plotStats(0, 5, "Robs vs Res"))

        self.IsigmaVsR= QtWidgets.QPushButton()
        self.IsigmaVsR.setText("<I/\u03C3I> vs Res")
        self.IsigmaVsR.clicked.connect(lambda:self.plotStats(0, 8, "<I/\u03C3I> vs Res"))

        self.SanoVsR= QtWidgets.QPushButton()
        self.SanoVsR.setText("Sig. Ano. vs Res")
        self.SanoVsR.clicked.connect(lambda:self.plotStats(0, 12, "Sig. Ano. vs Res"))

        self.barLayout.addWidget(self.ccVsR)
        self.barLayout.addWidget(self.compVsR)
        self.barLayout.addWidget(self.RobsVsR)
        self.barLayout.addWidget(self.SanoVsR)
        self.barLayout.addWidget(self.IsigmaVsR)

        self.tabLayout.addWidget(self.Title, 1, alignment=QtCore.Qt.AlignCenter)
        self.tabLayout.addWidget(self.buttonBar, 1)
        self.tabLayout.addWidget(self.statsBar, 1)
        self.tabLayout.addWidget(self.statsCanvas, 17)
        
    def plotStats(self, res, value, title):
        self.Ax.clear()
        plotDataX= []
        plotDataY= []
        plotList, _ = extractXSCALEStat(f"{self.ProcessedDir}/XSCALE.LP")

        #Check if plotList is empty
        if not plotList:
            self.Ax.text(0.5, 0.5, "No data available to plot", 
                         horizontalalignment='center', verticalalignment='center')
            self.statsCanvas.draw()
            return

        LowestRes = float(plotList[0][res]) + 0.25
        HighestRes = float(plotList[-1][res]) - 0.25 if float(plotList[-1][res]) > 0.25 else 0
        for line in plotList:
            plotDataX.append(float(line[res])) 
            plotDataY.append(float(line[value].strip('*').strip('%')))
        self.Ax.plot(plotDataX, plotDataY, 'r-^')
        self.Ax.set_xlim(LowestRes, HighestRes)
        self.Ax.set_title(title)
        self.statsCanvas.draw()


#plot multple statistics for multiple processed directories
class MultiPlotTab(QtWidgets.QWidget):
    def __init__(self, ProcessedDirs:list):
        super().__init__()
        self.ProcessedDirs= ProcessedDirs
        #set up widget
        self.tabLayout= QtWidgets.QVBoxLayout(self)
        self.Title=QtWidgets.QLabel(self)
        self.Title.setText(f"Multi statistics for {len(self.ProcessedDirs)} processed directories")
        self.Title.setStyleSheet("color: black; font-weight: bold; font-size: 14px;")

        #set up the plot
        self.MultiStatusPlots = Figure()
        self.Ax = self.MultiStatusPlots.add_subplot(111)
        self.MultiStatsCanvas = FigureCanvas(self.MultiStatusPlots)
        self.MultiplestatsBar= NavigationToolbar(self.MultiStatsCanvas, self)

        #Buttons to plot stats
        self.buttonBar = QtWidgets.QWidget()
        self.barLayout= QtWidgets.QHBoxLayout(self.buttonBar)

        self.ccVsR= QtWidgets.QPushButton()
        self.ccVsR.setText("CC. vs Res")
        self.ccVsR.clicked.connect(lambda:self.MultiPlotStats(0, 4, "CC. vs Res Multi"))

        self.compVsR= QtWidgets.QPushButton()
        self.compVsR.setText("comp vs Res")
        self.compVsR.clicked.connect(lambda:self.MultiPlotStats(0, 10, "comp vs Res Multi"))

        self.RobsVsR= QtWidgets.QPushButton()
        self.RobsVsR.setText("Robs vs Res")
        self.RobsVsR.clicked.connect(lambda:self.MultiPlotStats(0, 5, "Robs vs Res Multi"))

        self.IsigmaVsR= QtWidgets.QPushButton()
        self.IsigmaVsR.setText("<I/\u03C3I> vs Res")
        self.IsigmaVsR.clicked.connect(lambda:self.MultiPlotStats(0, 8, "<I/\u03C3I> vs Res Multi"))

        self.SanoVsR= QtWidgets.QPushButton()
        self.SanoVsR.setText("Sig. Ano. vs Res")
        self.SanoVsR.clicked.connect(lambda:self.MultiPlotStats(0, 12,"Sig. Ano. vs Res Multi"))

        self.barLayout.addWidget(self.ccVsR)
        self.barLayout.addWidget(self.compVsR)
        self.barLayout.addWidget(self.RobsVsR)
        self.barLayout.addWidget(self.SanoVsR)
        self.barLayout.addWidget(self.IsigmaVsR)

        self.tabLayout.addWidget(self.Title, 1, alignment=QtCore.Qt.AlignCenter)
        self.tabLayout.addWidget(self.buttonBar, 1)
        self.tabLayout.addWidget(self.MultiplestatsBar, 1)
        self.tabLayout.addWidget(self.MultiStatsCanvas, 17)


    #plot the statistics for the selected button
    def MultiPlotStats(self, res, value, title):
        self.Ax.clear()

        #store legends
        legend_handles = []

        #setup plot for each processed directory
        for ProcessedDir in self.ProcessedDirs:
            plotList, _ = extractXSCALEStat(f"{ProcessedDir}/XSCALE.LP")

            #Check if plotList is empty
            if not plotList:
                self.Ax.text(0.5, 0.5, f"No data available to plot for {Path(ProcessedDir).name}", 
                             horizontalalignment='center', verticalalignment='center')
                continue

            plotDataX= []
            plotDataY= []
            LowestResList = []
            HighestResList = []
            
            for line in plotList:
                #check whether x and v values are valid numbers
                try:
                    LowestResList.append(float(line[res]) + 0.25)
                    HighestResList.append(float(line[res]) - 0.25 if float(line[res]) > 0.5 else 0)
                    x_value = float(line[res]) if line[res] else 0
                    y_value = float(line[value].strip('*').strip('%')) if line[value] else 0
                    plotDataX.append(x_value)
                    plotDataY.append(y_value)
                except (ValueError, TypeError) as e:
                    print(f"Error converting data: {e} - line[{value}] = {line[value]}")
                    continue

            #add figure legend with the name of the processed directory
            line, = self.Ax.plot(plotDataX, plotDataY, '-^', label=Path(ProcessedDir).name)

            legend_handles.append(line)

        self.Ax.set_xlim(max(LowestResList), min(HighestResList))
        self.Ax.set_title(title)

        #setup the legend
        if legend_handles:
            self.Ax.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(1.02, 1),\
                           borderaxespad=0, fontsize="small", title_fontsize="medium")
            self.MultiStatusPlots.tight_layout()
        else:
            self.MultiStatusPlots.tight_layout()

        #plot the figure
        self.MultiStatsCanvas.draw()



#create pre-plot dendrogram tab to plot the dendrogram with given Threshold
class PrePlotDendrogram(QtWidgets.QWidget):
    def __init__(self, ccClusterfile:QtWidgets.QLineEdit, Threshold:QtWidgets.QLineEdit, setupCC_method):
        super().__init__()
        self.ccClusterfile_widget = ccClusterfile
        self.threshold_widget = Threshold
        self.setupCC_method = setupCC_method

        #setup widget
        self.tabLayout = QtWidgets.QVBoxLayout(self)
        self.Title=QtWidgets.QLabel(self)
        self.Title.setText(f"Dendrogram with threshold: {self.threshold_widget.text().strip() if self.threshold_widget else ''}")
        self.Title.setStyleSheet("color: black; font-weight: bold; font-size: 14px;")
    
        #set up the plot
        self.dendroPlot = Figure()
        self.Ax = self.dendroPlot.add_subplot(111)
        self.dendroCanvas = FigureCanvas(self.dendroPlot)
        self.dendrostatsBar= NavigationToolbar(self.dendroCanvas, self)

        #plot button
        self.plotButton = QtWidgets.QPushButton("Plot Dendrogram")
        self.plotButton.setFixedSize(200, 30)
        self.plotButton.clicked.connect(self.on_plot_clicked)

        #add widgets to layout
        self.tabLayout.addWidget(self.Title, 1, alignment=QtCore.Qt.AlignCenter)
        self.tabLayout.addWidget(self.plotButton, 1, alignment=QtCore.Qt.AlignCenter)
        self.tabLayout.addWidget(self.dendrostatsBar, 1)
        self.tabLayout.addWidget(self.dendroCanvas, 17)


    #check whether the threshold is valid float and plot the dendrogram
    def on_plot_clicked(self):
        #get current values from the widgets
        ccClusterfile = self.ccClusterfile_widget.text() if self.ccClusterfile_widget else None
        threshold_text = self.threshold_widget.text().strip() if self.threshold_widget else None
        
        if threshold_text:
            try:
                threshold = round(float(threshold_text), 2)
                self.PlotDendrogram(ccClusterfile, threshold)
            except ValueError:
                #Show error for invalid threshold
                self.Ax.clear()
                self.Ax.text(0.5, 0.5, "Invalid threshold value", 
                            horizontalalignment='center', verticalalignment='center')
                self.dendroCanvas.draw()
        else:
            self.PlotDendrogram(ccClusterfile, None)


    #prepare Dendrogram with current threshold
    def PlotDendrogram(self, ccClusterfile:str=None, threshold:float=None):
        #Clean plot
        self.Ax.clear()
        
        # Check parameters
        if ccClusterfile is None:
            self.Ax.text(0.5, 0.5, "No ccCluster file in Work dir", \
                         horizontalalignment='center', verticalalignment='center')
            self.dendroCanvas.draw()
            return None
        
        if threshold is None:
            self.Ax.text(0.5, 0.5, "No threshold provided", \
                         horizontalalignment='center', verticalalignment='center')
            self.dendroCanvas.draw()
            return None

        #set up CC
        _, Tree, etiquets, etiquetsWithNum, _ = self.setupCC_method(ccClusterfile)

        #check whether the dendrogram tree is valid
        if Tree is None:
            self.Ax.text(0.5, 0.5, "Error: Could not generate dendrogram",
                        horizontalalignment='center', verticalalignment='center')
            self.dendroCanvas.draw()
            return None

        X = hierarchy.dendrogram(Tree, color_threshold=threshold, ax=self.Ax)

        #try to match the color to cluster number using count
        #As they contains same amount od datasets
        #get the FlatC with the same threshold
        FlatC = hierarchy.fcluster(Tree, threshold, criterion='distance')
        cluster_counts = collections.Counter(FlatC)
        color_counts = collections.Counter(X['leaves_color_list'])

        #Remove unwanted color
        for remove_color in ['C0', 'k', 'grey']:
            if remove_color in color_counts:
                del color_counts[remove_color]

        #Remove cluster with only 1 dataset
        filtered_clusters = {cluster: count for cluster, count in cluster_counts.items() if count > 1}

        #Sort by count
        sorted_clusters = sorted(filtered_clusters.items(), key=lambda x: x[1])
        sorted_colors = sorted(color_counts.items(), key=lambda x: x[1])
        print(f"sorted cluster: {sorted_clusters}")
        print(f"sorted color: {sorted_colors}")

        #Match clusters to colors by count
        cluster_to_color = {}
        for (cluster, cluster_count), (color, color_count) in zip(sorted_clusters, sorted_colors):
            if cluster_count == color_count:
                cluster_to_color[cluster] = color
            else:
                print(f"Warning: Count mismatch! Cluster {cluster} has {cluster_count}, Color {color} has {color_count}")

        #Create legend using matched clusters
        legend_handles = []
        for cluster in sorted(cluster_to_color.keys()):
            color = cluster_to_color[cluster]
            legend_handles.append(mpatches.Patch(color=color, label=f"$\\mathbf{{Cluster\\ {cluster}}}$ :\n{cluster_counts[cluster]} datasets"))

        #Add legend
        self.Ax.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(1.02, 1),
                    borderaxespad=0, fontsize="small", handleheight=3, handlelength=2,
                    title_fontsize="medium")

        self.dendroPlot.tight_layout()
        self.Title.setText(f"Dendrogram with threshold: {threshold}")

        self.dendroCanvas.draw()
