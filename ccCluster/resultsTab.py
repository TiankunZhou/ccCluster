from __future__ import print_function

__author__ = "Gianluca Santoni"
__copyright__ = "Copyright 2015-2019"
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
from .clustering import Clustering

# implement the default mpl key bindings


import sys
import os
import subprocess
import matplotlib.patches as mpatches

#Read XSCALE.LP and extract information to plot the statistics
def extractXSCALEStat(XSCALEFile):
    plotList = []
    plotText = ""
    with open(XSCALEFile, 'r') as LogFile:
        for line in LogFile:
            if line.strip().startswith('LIMIT'):
                break
        for line in LogFile:
            break
        for line in LogFile:
            if line.strip().startswith('total'):
                break
            plotList.append(line.split())

        #align the columns, use the max length of each column to determine the width
        if plotList:
            col_widths = [max(len(row[col_idx]) for row in plotList) + 3 for col_idx in range(len(plotList[0]))]
            for line in plotList:
                aligned_line = "".join(f"{token:>{col_widths[i]}}" for i, token in enumerate(line))
                plotText += f"{aligned_line}\n"
    
    return plotList, plotText


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
        self.ccVsR.clicked.connect(lambda:self.plotStats(0, 4,"CC. vs Res" ))

        self.compVsR= QtWidgets.QPushButton()
        self.compVsR.setText("comp vs Res")
        self.compVsR.clicked.connect(lambda:self.plotStats(0, 10,"comp vs Res" ))

        self.RobsVsR= QtWidgets.QPushButton()
        self.RobsVsR.setText("Robs vs Res")
        self.RobsVsR.clicked.connect(lambda:self.plotStats(0, 5, "Robs vs Res" ))

        self.IsigmaVsR= QtWidgets.QPushButton()
        self.IsigmaVsR.setText("I/sigmaI vs Res")
        self.IsigmaVsR.clicked.connect(lambda:self.plotStats(0, 8,"I/sigmaI vs Res" ))

        self.SanoVsR= QtWidgets.QPushButton()
        self.SanoVsR.setText("Sig. Ano. vs Res")
        self.SanoVsR.clicked.connect(lambda:self.plotStats(0, 12,"Sig. Ano. vs Res" ))

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
        for line in plotList:
            plotDataX.append(float(line[res])) 
            plotDataY.append(float(line[value].strip('*').strip('%')))
        self.Ax.plot(plotDataX, plotDataY, 'r-^')
        self.Ax.set_xlim(10,0)
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
        self.ccVsR.clicked.connect(lambda:self.MultiPlotStats(0, 4,"CC. vs Res Multi"))

        self.compVsR= QtWidgets.QPushButton()
        self.compVsR.setText("comp vs Res")
        self.compVsR.clicked.connect(lambda:self.MultiPlotStats(0, 10,"comp vs Res Multi"))

        self.RobsVsR= QtWidgets.QPushButton()
        self.RobsVsR.setText("Robs vs Res")
        self.RobsVsR.clicked.connect(lambda:self.MultiPlotStats(0, 5, "Robs vs Res Multi"))

        self.IsigmaVsR= QtWidgets.QPushButton()
        self.IsigmaVsR.setText("I/sigmaI vs Res")
        self.IsigmaVsR.clicked.connect(lambda:self.MultiPlotStats(0, 8,"I/sigmaI vs Res Multi"))

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
            plotDataX= []
            plotDataY= []
            
            for line in plotList:
                #check whether x and v values are valid numbers
                try:
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

        self.Ax.set_xlim(10,0)
        self.Ax.set_title(title)

        #setup the legend
        if legend_handles:
            self.Ax.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(1.02, 1),\
                           borderaxespad=0, fontsize="small", title="Directories", title_fontsize="medium")
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
        _, Tree, _, _ = self.setupCC_method(ccClusterfile)

        #check whether the dendrogram tree is valid
        if Tree is None:
            self.Ax.text(0.5, 0.5, "Error: Could not generate dendrogram",
                        horizontalalignment='center', verticalalignment='center')
            self.dendroCanvas.draw()
            return None

        X = hierarchy.dendrogram(Tree, color_threshold=threshold, ax=self.Ax)

        #Show figure legend about what color is what cluster
        legend_handles = [mpatches.Patch(color=c, label=f"Cluster {i+1}") \
                            for i, c in enumerate(dict.fromkeys(X['color_list'])) \
                            if c not in ['C0', 'k', 'grey']]
        self.Ax.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(1.02, 1), \
                       borderaxespad=0, fontsize="small", title="Clusters", title_fontsize="medium")

        self.dendroPlot.tight_layout()
        self.Title.setText(f"Dendrogram with threshold: {threshold}")

        #plot the figure
        self.dendroCanvas.draw()
