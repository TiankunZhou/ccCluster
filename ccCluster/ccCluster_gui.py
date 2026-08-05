#! /usr/bin/env python3
from __future__ import print_function, absolute_import

__author__ = "Gianluca Santoni"
__copyright__ = "Copyright 2015-2019"
__credits__ = ["Gianluca Santoni, Alexander Popov"]
__license__ = ""
__version__ = "1.0"
__maintainer__ = "Gianluca Santoni"
__email__ = "gianluca.santoni@esrf.fr"
__status__ = "Beta"


#implement the default mpl key bindings
from PyQt5 import QtGui, QtCore, QtWidgets
from PyQt5.QtWidgets import QWidget, QApplication
import matplotlib.pyplot as plt
import sys
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.patches as mpatches

#import CalcClass
from scipy.cluster import hierarchy
import collections
import operator
from time import sleep
import os
from .resultsTab import resultsTab
from .summary import resultsSummary
from .clustering import Clustering
from .ccCalc import ccList

#Insert parse  to change the file path from command line
import argparse

#Deal with wild card
import glob
from pathlib import Path
import textwrap


#Startup message
print(r"""ccCluster - HCA for protein crystallography 
G. Santoni and A. Popov, 2015
              v .   ._, |_  .,
           `-._\/  .  \ /    |/_
               \\  _\, y | \//
         _\_.___\\, \\/ -.\||
           `7-,--.`._||  / / ,
           /'     `-. `./ / |/_.'
                     |    |//
                     |_    /
                     |-   |
                     |   =|
                     |    |
--------------------/ ,  . \--------._
""")


#Set color for printing
class colors:
    BOLD = '\033[1m'
    GREEN = '\033[32m'
    BLUE = '\033[34m'
    RED = '\033[31m'
    ENDC = '\033[m'



#here is the class that creates the tab for ccCal, 
#to generate ccClusterLog.txt in the target folder if it is not exist
class tab_ccCal(QWidget):
    def __init__(self, realTimeUpdates:MainWindow, **kwargs):
        #pass self to parent
        super().__init__()

        #get real time update
        self.realTimeUpdate = realTimeUpdates

        #set initial work dir as pwd and pass it to shared place
        self.WorkDir = os.getcwd()
        self.realTimeUpdate.updateWorkDir(self.WorkDir)

        #vertical layout
        self.ccCal_layout = QtWidgets.QVBoxLayout(self)

        #order is important, two lines below must above addwidget
        #set up general widget
        self.general_area()
        #set up buttons widget
        self.ccCal_area()

        #Add widget to the layout
        self.ccCal_layout.addWidget(self.general_widget, 1)
        self.ccCal_layout.addWidget(self.ccCal_widget, 4)


    #define general area
    def general_area(self):
        #define wedget
        self.general_widget = QtWidgets.QWidget()

        #show and select work dir
        self.WorkDir_entry = QtWidgets.QLineEdit()
        self.WorkDir_entry.setText(self.WorkDir)  
        #self.WorkDir_entry.setPlaceholderText("Please select a work dir")

        #update the work dir in real time when the text is changed
        self.WorkDir_entry.textChanged.connect(self.realTimeUpdat.updateWorkDir)
        #do not need label for work dir, as it is self-explanatory
        #self.WorkDir_entry.textChanged.connect(lambda: self.realTimeUpdat.updateWorkDir(self.WorkDir_entry.text()))
        #check ccClusterlog.txt when change the work dir:
        self.WorkDir_entry.textChanged.connect(self.check_ccCalLogStatus)
        #update the result list in ccCluster tab when work dir changed
        self.WorkDir_entry.textChanged.connect(self.realTimeUpdate.tab_ccCluster.CheckAndShowResult)

        #select workdir button
        self.ChooseWorkDir = QtWidgets.QPushButton("Select work dir")
        self.ChooseWorkDir.clicked.connect(self.select_workDir)

        #create layout and add the content (button etc)
        layout = QtWidgets.QVBoxLayout(self.general_widget)

        #workdir layout
        WorkDir_layout = QHBoxLayout()
        WorkDir_layout.addWidget(self.WorkDir_entry, 4)
        WorkDir_layout.addWidget(self.ChooseWorkDir, 1)

        #put in main layout
        layout.addLayout(WorkDir_layout)


    #define ccCal area
    def ccCal_area(self):
        #define widget
        self.ccCal_widget = QtWidgets.QWidget()

        #check whether ccClusterlog.txt exist and show it in the log bar
        self.ccCallog_status = QtWidgets.QLabel()
        self.check_ccCalLogStatus()

        #select paths of HKL files for ccCal
        self.HKLPaths_text = QtWidgets.QTextEdit()
        self.HKLPaths_text.setPlaceholderText(textwrap.dedent(f"""
                                        Please put absolute path of HKL files here
                                        e.g. 
                                        /gpfs/easy/data/id30a3/inhouse/tiankun/mesh_n_collect/test6/codgas_ccCluster/0.HKL
                                        /gpfs/easy/data/id30a3/inhouse/tiankun/mesh_n_collect/test6/codgas_ccCluster/1.HKL
                                        or use wild card
                                        /gpfs/easy/data/id30a3/inhouse/tiankun/mesh_n_collect/test7/HKL_files/HKL*/*.HKL
                                        Using wild card is highly recommended
                                        """))
        self.insert_HKLPaths = QtWidgets.QPushButton("Add a HKL file")
        self.insert_HKLPaths.clicked.connect(self.add_HKLPath)

        #button to run ccCla
        self.run_ccCal = QtWidgets.QPushButton("Run ccCal")
        self.run_ccCal.clicked.connect(self.submit_ccCal)

        #create main layout
        layout = QtWidgets.QVBoxLayout(self.ccCal_widget)

        #put hkltext and button together
        hklLayout = QtWidgets.QHBoxLayout()
        hklLayout.addWidget(self.HKLPaths_text, 4)     
        hklLayout.addWidget(self.insert_HKLPaths, 1)

        #other things layout
        layout.addWidget(self.ccCallog_status)
        #have hkllayout inside main layout
        layout.addLayout(hklLayout)
        layout.addWidget(self.run_ccCal)


    #get HKL file list for ccCal job
    def getHKLList (self):
        abs_file_list = []
        input_path = [line.strip() for line in self.HKLPaths_text.toPlainText().splitlines() if line.strip()]
        for path in input_path:
            hkl_files = glob.glob(path)
            if HKL_files:
                for hkl_file in hkl_files:
                    if os.path.isfile(hkl_file) and Path(hkl_file).suffix.lower() == ".hkl":
                        print(f"{colors.BLUE}Adding {hkl_file} to the HKL merge list")
                        abs_path_list.append(os.path.abspath(hkl_file))
                    else:
                        print(f"{colors.RED}No HKL file: {hkl_file}, please check{colors.ENDC}")
            else:
                print(f"{colors.RED}No files or folder in {path}, please check{colors.ENDC}")
        
        return abs_file_list


    #submit ccCal jobs
    def submit_ccCal(self):
        HKL_list = self.getHKLList()
        if os.path.isdir(self.realTimeUpdate.shareWorkDir) and if HKL_list:
            ccList(HKL_list, self.realTimeUpdate.shareWorkDir)
        elif not of.path.isdir(self.realTimeUpdate.shareWorkDir):
            print(f"{color.RED}Working dir does not exist, please check: {self.realTimeUpdate.shareWorkDir}{colors.ENDC}")
        elif not HKL_list:
            print(f"{colors.RED}No HKL file list, please check HKL paths{colors.ENDC}")
        else:
            print(f"{colors.RED}Unknow problem with ccCal input HKLs or working dir, please check{colors.ENDC}")
        
        #check ccClusterLog.txt after generation
        self.check_ccCalLogStatus()


    #select output folder
    def select_WorDir(self):
        WorkFolder = QtWidgets.QFileDialog.getExistingDirectory(self, "Select work dir")

        if WorkFolder:
            self.WorkDir_entry.setText(WorkFolder)


    #add HKL path 
    def add_HKLPath(self):
        file_path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select HKL File", "", "HKL Files (*.HKL *.hkl);;All Files (*)")

        if file_path:
            self.HKLPaths_text.append(os.path.abspath(file_path))


    #check whether ccClusterlog.txt exists in the work dir
    def check_ccCalLogStatus(self):
        ccClusterLog = os.path.join(self.realTimeUpdate.shareWorkDir, "ccClusterLog.txt")

        #update the ccClusterLog.txt status
        if os.path.isfile(ccClusterLog):
            self.ccCallog_status.setText(f"ccClusterLog.txt found: {ccClusterLog}")
            self.ccCallog_status.setStyleSheet("color: green; font-weight: bold")
        else:
            self.ccCallog_status.setText(f"ccClusterLog.txt not found in {self.realTimeUpdate.shareWorkDir}. Please generate one or select a diffenent path")
            self.ccCallog_status.setStyleSheet("color: red; font-weight: bold")



#define Class for the tab for ccClustering
#tabs for result and summary are generated 
#through different modules
class tab_ccCluster(QWidget):
    def __init__(self, realTimeUpdates:MainWindow, **kwargs):
        super().__init__()

        #get real time update
        self.realTimeUpdate = realTimeUpdates

        #read processed result in WorkDir:
        self.MergeResult = []
        self.CheckAndShowResult()

        #vertical layout
        self.ccCal_layout = QtWidgets.QVBoxLayout(self)

        #set up buttons
        self.ccClusterSetup_area()
        #set up plot area
        self.plotDendroAndStatistic_area()

        #other parameters

        #Add widget to the layout
        self.ccCal_layout.addWidget(self.ccClusterSetup_widget, 1)
        self.ccCal_layout.addWidget(self.plotDendroAndStatistic_widget, 4)


    #setups for ccCluster
    def ccClusterSetup_area(self):
        #Define widget
        self.ccClusterSetup_widget = QtWidgets.QWidget()

        #status bar and button to run ccCluster
        self.ccClusterLogPath_text = QtWidgets.QLineEdit()
        if os.path.isfile(f"{self.realTimeUpdate.shareWorkDir}/ccClusterLog.txt"):
            self.ccClusterLogPath_text.setText(f"{self.realTimeUpdate.shareWorkDir}/ccClusterLog.txt")
        else:
            self.ccClusterLogPath_text.setPlaceholderText(f"No ccClusterLog.txt in working Dir, please change working dir or generate one")
        self.insert_ccClusterLogPath = QtWidgets.QPushButton("select log file")
        self.insert_ccClusterLogPath.clicked.connect(self.select_ccClusterLogPath)

        #button to show auto threshold result and add it to the Threshold line
        #also show the largest group the default threshold
        self.ShowThreshold = QtWidgets.QLineEdit()
        self.ShowThreshold.setPlaceholderText("Put threshod or click Auto Threshold")
        self.AutoThreshols = QtWidgets.QPushButton("Auto Threshold")
        self.AutoThreshols.clicked.connect(self.getAutoThreshold)

        #button to show the largest group number with current threshold
        self.ShowLargestGroup = QtWidgets.QLineEdit()
        self.ShowLargestGroup.setPlaceholderText("Show largest group with current threshold")
        self.CheckLargestGroup = QtWidgets.QPushButton("Check Largest Group")
        self.CheckLargestGroup.clicked.connect(self.getLargestGroup)

        #anomalous flag
        self.anomBox = QtWidgets.QCheckBox("Anomalous data")
        self.anomBox.setChecked(False)

        #status bar to show information and button to run ccCluster job
        self.ccClusterStatusBar = QtWidgets.Qlabel()
        self.update_ccClusterStatusBar("ready to work")
        self.RunccCluster = QtWidgets.QPushButton("Run ccCluster")
        self.RunccCluster.clicked.connect(self.submit_ccCluster)

        #create main layout
        layout = QtWidgets.QVBoxLayout(self.ccClusterSetup_widget)

        #put ccClusterLogPath and button together
        ccClusterLogLayout = QtWidgets.QHBoxLayout()
        ccClusterLogLayout.addWidget(self.ccClusterLogPath_text, 4)     
        ccClusterLogLayout.addWidget(self.insert_ccClusterLogPath, 1)

        #put AutoThreshold and show largest group together
        AutoThresholsLayout = QtWidgets.QHBoxLayout()
        AutoThresholsLayout.addWidget(self.ShowThreshold, 2)
        AutoThresholsLayout.addWidget(self.AutoThreshols, 1)
        AutoThresholsLayout.addWidget(self.ShowLargestGroup, 2)
        AutoThresholsLayout.addWidget(self.CheckLargestGroup, 1)
        AutoThresholsLayout.addWidget(self.anomBox, 1)

        #put status bar and merge button together
        StatusBarLayout = QtWidgets.QHBoxLayout()
        StatusBarLayout.addWidget(self.ccClusterStatusBar, 4)
        StatusBarLayout.addWidget(self.RunccCluster, 1)
        #pack everything in the layout
        layout.addLayout(ccClusterLogLayout)
        layout.addLayout(AutoThresholsLayout)
        layout.addLayout(StatusBarLayout)
        layout.addWidget(self.run_ccCal)


    #plot the Dendrogram and SCALE.LP statistics in different sub-tabs
    def plotDendroAndStatistic_area(self):
        #define widget
        self.plotDendroAndStatistic_widget = QtWidgets.QWidget()

        #Create layout
        layout = QtWidgets.QVBoxLayout(self.plotDendroAndStatistic_widget)

        #Create tabs for Dendrogram and statistics
        self.DendroAndStatsPlot()


    #Add ccCluster log path is needed
    def select_ccClusterLogPath(self):
        file_path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select ccCluster log File", "", "ccCluster Log Files (*.txt);;All Files (*)")

        if file_path:
            self.WorkDir_entry.setText(os.path.abspath(file_path))


    #Auto elect ccCluster log file if exists in work Dir
    def auto_select_ccClusterLogPath(self):
        ccClusterLog = os.path.join(self.workDir.text(), "ccClusterLog.txt")

        #update the ccClusterLog.txt status
        if os.path.isfile(ccClusterLog):
            self.ccCallog_status.setText(f"ccClusterLog.txt found: {ccClusterLog}")
            self.ccCallog_status.setStyleSheet("color: green; font-weight: bold")
        else:
            self.ccCallog_status.setText(f"ccClusterLog.txt not found in {self.workDir.text()}. Please generate one or select a diffenent path")
            self.ccCallog_status.setStyleSheet("color: red; font-weight: bold")


    #update ccCluster bar:
    def update_ccClusterStatusBar(self, status:str):
        self.ccClusterStatusBar.setText(f"{status}")
        #self.ccCallog_status.setStyleSheet("color: green; font-weight: bold")

    
    #auto define thershold
    def getAutoThreshold(self):
        CC, _, _, statue_text = self.realTimeUpdate.setupCC(self.ccClusterLogPath_text.text())
        if CC is None:
            self.update_ccClusterStatusBar(statue_text)
        else:
            Threshold = CC.thrEstimation()
            self.ShowThreshold.set(Threshold)
            GroupNum, largestGroup, totalHKL = CC.checkMultiplicity(Threshold)
            self.update_ccClusterStatusBar(f"Auto Threshold is {Threshold}, the largest cluster number is {GroupNum} with {largestGroup}/{totalHKL} files")


    #Show largest cluster
    def getLargestGroup(self):
        CC, _, _, status_text= self.realTimeUpdate.setupCC(self.ccClusterLogPath_text.text())
        if CC is None:
            self.update_ccClusterStatusBar(status_text)
        else:
            if self.ShowThreshold:
                GroupNum, largestGroup, totalHKL = CC.checkMultiplicity(self.ShowThreshold) # need change
                self.update_ccClusterStatusBar(f"Auto Threshold is {self.ShowThreshold}, the largest cluster number is {GroupNum} with {largestGroup}/{totalHKL} files")
            else:
                self.update_ccClusterStatusBar(f"Please input a threshold value")


    #submit ccCluster job
    def submit_ccCluster(self):
        CC, Tree, etiquets, status_text = self.realTimeUpdate.setupCC(self.ccClusterLogPath_text.text())
        threshold = self.ShowThreshold.text()
        if CC is None:
            self.update_ccClusterStatusBar(status_text)
        else:
            #check file type
            fileType = CC.inputType()
            #check anomalous flag
            if self.anomBox.isChecked():
                anomlous = "ano"
            else:
                anomlous = "no_ano"

            if fileType=="HKL":
                #prepare and run XSCALE job
                if not self.reference_HKL:
                    xscale_checker, xscale_path = CC.prepareXSCALE(anomlous, threshold)
                elif os.path.isfile(self.reference_HKL):
                    xscale_checker, xscale_path = CC.prepareXSCALE(anomlous, threshold, refHKL=self.reference_HKL)
                else:
                    xscale_checker, xscale_path = CC.prepareXSCALE(anomlous, threshold)
                if xscale_checker == True:
                    self.update_ccClusterStatusBar(f"Running XSCALE job in {xscale_path}")
                    CC.scaleAndMerge(anomlous, threshold, xscale_path)
                    #get jason from XSCALE
                    CC.flatClusterPrinter(threshold, etiquets, anomlous, xscale_path)

                #prepare and run Pointless
                pointless_checker, pointless_path = CC.preparePointless(anomlous, threshold)
                if pointless_checker == True:
                    self.update_ccClusterStatusBar(f"Running Pointless job in {pointless_path}")
                    CC.pointlessRun(anomlous, threshold, pointless_path)
                    #prepare and run aimless
                    self.update_ccClusterStatusBar(f"Running Aimless job in {xscale_path}, please be patient")
                    CC.aimlessRun(anomlous, threshold, pointless_path)

            #CC.passOInfoToGA(threshold, etiquets, anomlous)
            elif fileType=="mtz":
                #prepare and run Pointless
                pointless_checker, pointless_path = CC.preparePointless(anomlous, threshold)
                if pointless_checker == True:
                    self.update_ccClusterStatusBar(f"Running Pointless job in {pointless_path}")
                    CC.pointlessRun(anomlous, threshold, pointless_path)
                    #prepare and run aimless
                    self.update_ccClusterStatusBar(f"Running Aimless job in {xscale_path}, please be patient")
                    CC.aimlessRun(anomlous, threshold, pointless_path)
                    CC.flatClusterPrinter(threshold, etiquets, anomlous, pointless_path)
            else:
                self.update_ccClusterStatusBar(f"Unknown input file format, please check ccCluster log file: {self.ccClusterLogPath_text.text()}")
                print(f"Unknown input file format, please check ccCluster log file: {self.ccClusterLogPath_text.text()}")


    #check and show results will also be used by other tabs
    def CheckAndShowResult(self):
        abs_FolderPaths = Path(self.realTimeUpdate.shareWorkDir)
        for folder_path in abs_FolderPaths.glob("cc_Cluster_*"):
            if folder_path.is_dir():
                if (folder_path/"XSCALE.LP").is_file():
                    folder_name = folder_path.name
                    if folder_name not in self.MergeResult:
                        self.MergeResult.append(folder_name)
                else:
                    print(f"{colors.RED}No XSCALE.LP in folder: {folder_path}{colors.ENDC}")
            else:
                print(f"{colors.RED}Folder path does not exist: {folder_path}{colors.ENDC}")

        self.MergeResult.sort()


    #create tabs in the Result tab to show dendrogram and statistics for each merged result
    def DendroAndStatsPlot(self):
        #Show Dendrogram with current threshold with a button
        def ShowDendrogramThreshold(threshold):
            fig, ax = plt.subplots()

            CC, Tree, _, _ = self.realTimeUpdate.setupCC(self.ccClusterLogPath_text.text())
            X = hierarchy.dendrogram(Tree, color_threshold=threshold, ax=ax)

            #Show figure legend about what color is what cluster
            legend_handles = [mpatches.Patch(color=c, label=f"Cluster {i+1}") \
                                for i, c in enumerate(dict.fromkeys(X['color_list'])) \
                                if c not in ['C0', 'k', 'grey']]
            plt.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(1.02, 1), \
                        borderaxespad=0, fontsize="small", title="Clusters", title_fontsize="medium")

            fig.tight_layout()

            return FigureCanvas(fig)


        #set up the Dendrogram plot tab
        def DendrogramPlotProcessed(DendragramPath):
            DendroproicessedWedget = QtWidgets.QWidget()
            Dendrolayout = QtWidgets.QVBoxLayout(DendroproicessedWedget)

            #plot the Dendrogram from the pong file
            self.ImageBox = QtWidgets.QLabel(self)
            self.ImageBox.setPixmap(f"{DendragramPath}")
            self.ImageBox.setScaledContents(True)

            Dendrolayout.addWidget(self.ImageBox)

            return DendroproicessedWedget
        
        #grab the mergnig statistcs form XSCALE.LP
        def XSCALEStat(XSCALEFile):
            XSCALEstat = QtWidgets.QWidget()
            XSCALEstatlayout = QtWidgets.QVBoxLayout(XSCALEstat)

            #create a read-only text edit to show the XSCALE.LP statistics
            XSCALEText = QtWidgets.QTextEdit()
            XSCALEText.setReadOnly(True)

            #Get the content
            _, plotText = self.extractXSCALEStat(XSCALEFile)
            XSCALEText.setText(plotText)

            #put in to layout
            XSCALEstatlayout.addWidget(XSCALEText)

            #return the widget
            return XSCALEstat
            
        #set up the widget and the tabs
        self.plotDendroAndStatistic_widget = QtWidgets.QWidget()

        #setup the tab for plotting Dendrogram with given threshold
        DendroGramTreshTab = QtWidgets.QWidget()
        DendroGramTreshLayout = QtWidgets.QVBoxLayout(DendroGramTreshTab)
        DendroGramTreshPlot = ShowDendrogramThreshold(float(self.ShowThreshold.text())) if self.ShowThreshold.text() else None

        #set up the tab for Dendrogram and statistics
        if not self.MergeResult:
            print(f"{colors.RED}No merged result found in {self.realTimeUpdate.shareWorkDir}, please check{colors.ENDC}")
            return

        #set up the tab for Dendrogram and statistics
        for result_folder in self.MergeResult:
            resultTab = QtWidgets.QWidget()
            resultTabLayout = QtWidgets.QVBoxLayout(resultTab)


    #Extrace statistics from XSCALE.LP
    def extractXSCALEStat(self, XSCALEFile):
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

            #align the columns
            if plotList:
                col_widths = [max(len(row[col_idx]) for row in plotList) + 3 for col_idx in range(len(plotList[0]))]
                for line in plotList:
                    aligned_line = "".join(f"{token:>{col_widths[i]}}" for i, token in enumerate(line))
                    plotText += f"{aligned_line}\n"

        return plotList, plotText
                

"""
    #old not use
        ###########
        #Buttons Widget:
        #to contain all the buttons and inputs
        #to run clusterings
        ###########
        #Set the main window
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.centralWidgetLayout= QtWidgets.QHBoxLayout(self.centralwidget)
        #set button window

        self.verticalLayoutWidget = QtWidgets.QWidget(self.centralwidget)
        self.verticalLayoutWidget.setMaximumSize(200, 191)
        self.verticalLayoutWidget.setObjectName(_fromUtf8("verticalLayoutWidget"))
        self.verticalLayout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        #self.verticalLayout.setMargin(0)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.PlotButton = QtWidgets.QPushButton(self.verticalLayoutWidget)
        self.PlotButton.setObjectName(_fromUtf8("PlotButton"))
        self.PlotButton.clicked.connect(self.createDendrogram)

        self.verticalLayout.addWidget(self.PlotButton)
        self.lineEdit = QtWidgets.QLineEdit(self.verticalLayoutWidget)
        self.lineEdit.setObjectName(_fromUtf8("lineEdit"))
        self.verticalLayout.addWidget(self.lineEdit)
        self.checkBox = QtWidgets.QCheckBox(self.verticalLayoutWidget)
        self.checkBox.setChecked(True)
        self.checkBox.setObjectName(_fromUtf8("checkBox"))
        self.verticalLayout.addWidget(self.checkBox)
        self.anomBox = QtWidgets.QCheckBox(self.verticalLayoutWidget)
        self.anomBox.setChecked(True)
        self.checkBox.setObjectName(_fromUtf8("anomBox"))
        self.verticalLayout.addWidget(self.anomBox)
        self.processButton = QtWidgets.QPushButton(self.verticalLayoutWidget)
        self.processButton.setObjectName(_fromUtf8("processButton"))
        self.processButton.clicked.connect(self.processClusters)
        self.verticalLayout.addWidget(self.processButton)
        self.summaryButton= QtWidgets.QPushButton(self.verticalLayoutWidget)
        self.summaryButton.setObjectName(_fromUtf8("summaryButton"))
        self.summaryButton.clicked.connect(self.showSummary)
        self.verticalLayout.addWidget(self.summaryButton)
        
        self.centralWidgetLayout.addWidget(self.verticalLayoutWidget)

        ###########
        #Tab Widget:
        #to show Dendrogram and clustering results
        #Each Run of XSCALE will create a new tab for results
        ###########
        self.verticalLayoutWidget_2 = QtWidgets.QWidget(self.centralwidget)
        self.verticalLayoutWidget_2.setMinimumSize(1020, 410)
        self.verticalLayoutWidget_2.setObjectName(_fromUtf8("verticalLayoutWidget_2"))
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.verticalLayoutWidget_2)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.tabWidget = QtWidgets.QTabWidget(self.verticalLayoutWidget_2)
        self.tabWidget.setTabPosition(QtWidgets.QTabWidget.North)
        self.tabWidget.setObjectName(_fromUtf8("tabWidget"))

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.Dendrogram = plt.figure()

        self.TreeCanvas = FigureCanvas(self.Dendrogram)
        cid = self.Dendrogram.canvas.mpl_connect('button_press_event', self.getThreshold)
        self.TreeBar= NavigationToolbar(self.TreeCanvas, self)
        self.createDendrogram()

        self.plotTab = QtWidgets.QWidget()
        self.plotTabLayout=QtWidgets.QVBoxLayout(self.plotTab)
        self.plotTabLayout.addWidget(self.TreeBar)
        self.plotTabLayout.addWidget(self.TreeCanvas)
        self.plotTab.setObjectName(_fromUtf8("plotTab"))
        self.tabWidget.addTab(self.plotTab, _fromUtf8("Dendrogram"))
        self.verticalLayout_2.addWidget(self.tabWidget)
        self.centralWidgetLayout.addWidget(self.verticalLayoutWidget_2)

        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        #self.menubar.setGeometry(QtCore.QRect(0, 0, 1215, 21))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        self.statusbar.showMessage('Ready!')
        MainWindow.setStatusBar(self.statusbar)
"""
        ##########
        #Show previous results
        # to be removed from the setupUi class
        ###########
        """
        self.alreadyDone= []
        if os.path.isfile(os.getcwd()+'/.cc_cluster.log'):
            with open(os.getcwd()+'/.cc_cluster.log') as log:
                for line in log:
                    L = line.split(',')
                    try :
                        self.tabWidget.addTab(resultsTab(float(L[1]), L[2].strip(), L[3].strip(), L[4]),L[1]+L[2]+L[3].strip())
                        self.alreadyDone.append([L[1], L[2].strip(), L[3].strip()])

                    except:
                        self.tabWidget.addTab(resultsTab(float(L[1]), L[2].strip(), L[3].strip(), 'unk'),L[1]+L[2]+L[3].strip())
                        self.alreadyDone.append([L[1], L[2].strip(), L[3].strip()])

        #L[1] = threshold, L[2]=number, L[3] = anomFlag
        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)



    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "Cluster and Merge", None))
        self.PlotButton.setText(_translate("MainWindow", "Plot Dendrogram", None))
        self.checkBox.setText(_translate("MainWindow", "Merge only biggest cluster", None))
        self.anomBox.setText(_translate("MainWindow", "Anomalous data", None))
        self.processButton.setText(_translate("MainWindow", "Merge clusters", None))
        self.summaryButton.setText(_translate("MainWindow", "Summary", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.TreeCanvas), _translate("MainWindow", "Dendrogram", None))
        #self.tabWidget.setTabText(self.tabWidget.indexOf(self.clusterTab), _translate("MainWindow", "Clustering result", None))

    """


    def showSummary(self):
        self.sum = resultsSummary()
        self.sum.show()


    def onChanged(self, text):
        self.threshold = float(text)



class tab_plotStats():
    def __init__(self, realTimeUpdates:MainWindow, **kwargs):
        super().__init__()

        #get real time update from self of MainWindow
        self.realTimeUpdate = realTimeUpdates

        #get self from ccCluster tab to get the result list
        self.ccClusterTab = self.realTimeUpdate.Tab_ccCluster
        #setup buttons


#put the tabs to gether in one GUI
class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, **kwargs):
        #make one or several linEdit or textEdit accessable to all tabs in realtime
        self.shareWorkDir = ""

        ###pass the self of MainWindow to the tabs as argument
        ###to enalbe realtime update on lineEdit
        self.Tab_ccCal = tab_ccCal(self)
        self.Tab_ccCluster = tab_ccCluster(self)
        self.Tab_plotStats = tab_plotStats(self)

        #add tabs
        self.tabWidget.addTab(self.Tab_ccCal, "ccCal tab")
        self.tabWidget.addTab(self.Tab_ccCluster, "ccCluster tab")
        self.tabWidget.addTab(self.Tab_compare, "Compare merging result tab")

        #Set up main window
        MainWindow.setObjectName("MainWindow")
        MainWindow.setWindowTitle("Cluster and merge")
        MainWindow.resize(1600, 1200)
    
    #update workDir, put it here to avoid duplication of same function
    #use it in other class: self.realTimeUpdate.updateWorkDit(WorkDir)
    def updateWorkDir(self, text):
        if not text.strip():
            print(f"{colors.RED}Work dir path is empty.{colors.ENDC}")
            return #End function here (earlier)

        abs_dir = os.path.abspath(text)
        if os.path.isdir(abs_dir):
            self.shareWorkDir = abs_dir
        else:
            print(f"{colors.RED}Update work dir failed, input dir does not exist: {abs_dir}{colors.ENDC}")
    

    #set up clustering scripts:
    def setupCC(self, correlationFile:str):
        if not correlationFile:
            text = f"No ccClusterlog file provided"
            return None, None, None, text

        if not os.path.isfile(correlationFile):
            text = f"ccClusterlog file does not exist: {correlationFile}"
            return None, None, None, text

        if not self.shareWorkDir or not os.path.isdir(self.shareWorkDir):
            text = f"Work directory not set or invalid: {self.shareWorkDir}"
            return None, None, None, text

        CC = Clustering(correlationFile, self.shareWorkDir)
        Tree = CC.avgTree() #needed, for set up self.Tree in clustering.py
        etiquets = CC.createLabels()
        return CC, Tree, etiquets, text



#Main part of the program
#with the different options, we can chose 
# to process through the shell,
#count the multiplicity of the highest cluster
#run the interface
def main():
    app = QApplication(sys.argv)
    ex = MainWindow()
    ex.show()
    sys.exit(app.exec_())      


#run the GUI
if __name__== '__main__':
    main()
