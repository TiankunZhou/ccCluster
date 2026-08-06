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
from Pathlib import Path

# implement the default mpl key bindings


import sys
import os
import subprocess

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
        QtWidgets.QWidget.__init__(self)
        self.ProcessedDir= ProcessedDir
        #setup widget
        self.tabLayout= QtWidgets.QVBoxLayout(self)
        self.Title=QtWidgets.QLabel(self)
        ResultName = Path(self.ProcessedDir).name.split('_')
        self.Title.setText(f"Threshold: {ResultName[2]} and Group number: {ResultName[3]}")

        #set up the plot
        self.statsPlot, self.Ax= plt.subplots()
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

        self.tabLayout.addWidget(self.Title)
        self.tabLayout.addWidget(self.buttonBar)
        self.tabLayout.addWidget(self.statsBar)
        self.tabLayout.addWidget(self.statsCanvas)
        
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

        #set up the plot
        self.MultiStatusPlots, self.Ax= plt.subplots()
        self.MultiStatsCanvas = FigureCanvas(self.MultiStatusPlots)
        self.statsBar= NavigationToolbar(self.MultiStatsCanvas, self)

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

        self.tabLayout.addWidget(self.Title)
        self.tabLayout.addWidget(self.buttonBar)
        self.tabLayout.addWidget(self.statsBar)
        self.tabLayout.addWidget(self.MultiStatsCanvas)


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
            self.MultiStatusPlots.tight_layout(rect=[0, 0, 0.82, 1])
        else:
            self.MultiStatusPlots.tight_layout()

        #plot the figure
        self.MultiStatsCanvas.draw()
