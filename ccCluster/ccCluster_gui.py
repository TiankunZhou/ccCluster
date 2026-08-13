#! /usr/bin/env python3
from __future__ import print_function, absolute_import

__author__ = "Gianluca Santoni & Tiankun Zhou"
__copyright__ = "Copyright 2015-2026"
__credits__ = ["Gianluca Santoni, Alexander Popov"]
__license__ = ""
__version__ = "1.0"
__maintainer__ = "Tiankun Zhou"
__email__ = "tiankun.zhou@esrf.fr"
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
from .resultsTab import SinglePlotTab, MultiPlotTab, PrePlotDendrogram
from .summary import resultsSummary
from .clustering import Clustering, extractXSCALEStat
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

        #vertical layout
        self.ccCal_layout = QtWidgets.QVBoxLayout(self)

        #order is important, two lines below must above addwidget
        #set up general widget
        self.general_area()
        #set up buttons widget
        self.ccCal_area()

        #separation line:
        Septline = QtWidgets.QFrame()
        Septline.setFrameShape(QtWidgets.QFrame.HLine)
        Septline.setStyleSheet("background-color: #888888; max-height: 3px; border: none; margin: 5px 0;")

        #Add widget to the layout
        self.ccCal_layout.addWidget(self.general_widget, 1)
        self.ccCal_layout.addWidget(Septline)
        self.ccCal_layout.addWidget(self.ccCal_widget, 19)


    #define general area
    def general_area(self):
        #define widget
        self.general_widget = QtWidgets.QWidget()

        #show and select work dir
        self.WorkDir_entry = QtWidgets.QLineEdit()
        self.WorkDir_entry.setText(self.realTimeUpdate.WorkDir)  
        #self.WorkDir_entry.setPlaceholderText("Please select a work dir")

        #update the work dir in real time when the text is changed
        self.WorkDir_entry.textChanged.connect(self.realTimeUpdate.updateWorkDir)
        #do not need label for work dir, as it is self-explanatory
        #self.WorkDir_entry.textChanged.connect(lambda: self.realTimeUpdate.updateWorkDir(self.WorkDir_entry.text()))
        #check ccClusterlog.txt when change the work dir:
        self.WorkDir_entry.textChanged.connect(self.check_ccCalLogStatus)
        #update the result list in ccCluster tab when work dir changed
        self.WorkDir_entry.textChanged.connect(self.realTimeUpdate.Tab_ccCluster.UpdateResultAndSyncTabs)
        self.WorkDir_entry.textChanged.connect(self.realTimeUpdate.Tab_ccCluster.auto_select_ccClusterLogPath)

        #select workdir button
        self.ChooseWorkDir = QtWidgets.QPushButton("Select work dir")
        self.ChooseWorkDir.clicked.connect(self.select_WorkDir)

        #create layout and add the content (button etc)
        GeneralLayout = QtWidgets.QVBoxLayout(self.general_widget)

        #workdir layout
        WorkDir_layout = QtWidgets.QHBoxLayout()
        WorkDir_layout.addWidget(self.WorkDir_entry, 4)
        WorkDir_layout.addWidget(self.ChooseWorkDir, 1)

        #setup title
        self.generalTitle = QtWidgets.QLabel("General settings")
        self.generalTitle.setStyleSheet("color: black; font-weight: bold; font-size: 16px;")

        #put in main layout
        GeneralLayout.addWidget(self.generalTitle, 1, alignment=QtCore.Qt.AlignCenter)
        GeneralLayout.addLayout(WorkDir_layout)


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

        #target folder for searching HKL files, default is pwd
        self.SearchDir = QtWidgets.QLineEdit()
        self.SearchDir.setPlaceholderText("Please select a HKL dir")
        self.FileName = QtWidgets.QLineEdit()
        self.FileName.setText("XSCALE.HKL")
        self.search_HKL_button = QtWidgets.QPushButton("Search HKL files")
        self.search_HKL_button.clicked.connect(self.search_HKL_files)

        #button to run ccCla
        self.run_ccCal = QtWidgets.QPushButton("Run ccCal")
        self.run_ccCal.setFixedSize(200, 50)
        self.run_ccCal.clicked.connect(self.submit_ccCal)

        #create main layout
        ccCalLayout = QtWidgets.QVBoxLayout(self.ccCal_widget)

        #put hkltext and button together
        self.hkladdfileswidget = QtWidgets.QWidget()
        self.hkladdfileslayout = QtWidgets.QHBoxLayout(self.hkladdfileswidget)   
        self.hkladdfileslayout.addWidget(self.insert_HKLPaths, 1)
        self.hkladdfileslayout.addWidget(self.SearchDir, 4)
        self.hkladdfileslayout.addWidget(self.FileName, 2)
        self.hkladdfileslayout.addWidget(self.search_HKL_button, 1)

        #setup title
        self.ccCalTitle = QtWidgets.QLabel("Check and generate ccClusterlog.txt file")
        self.ccCalTitle.setStyleSheet("color: black; font-weight: bold; font-size: 16px;")

        #other things layout
        ccCalLayout.addWidget(self.ccCalTitle, alignment=QtCore.Qt.AlignCenter)
        ccCalLayout.addWidget(self.ccCallog_status, alignment=QtCore.Qt.AlignCenter)

        #have hkllayout inside main layout
        ccCalLayout.addWidget(self.hkladdfileswidget)
        ccCalLayout.addWidget(self.HKLPaths_text)
        ccCalLayout.addWidget(self.run_ccCal, alignment=QtCore.Qt.AlignCenter)


    #get HKL file list for ccCal job
    def getHKLList (self):
        abs_file_list = []
        input_path = [line.strip() for line in self.HKLPaths_text.toPlainText().splitlines() if line.strip()]
        for path in input_path:
            hkl_files = glob.glob(path)
            if hkl_files:
                for hkl_file in hkl_files:
                    if os.path.isfile(hkl_file) and Path(hkl_file).suffix.lower() == ".hkl":
                        print(f"{colors.GREEN}Adding {hkl_file} to the HKL merge list{colors.ENDC}")
                        abs_file_list.append(os.path.abspath(hkl_file))
                    else:
                        print(f"{colors.RED}No HKL file: {hkl_file}, please check{colors.ENDC}")
            else:
                print(f"{colors.RED}No files or folder in {path}, please check{colors.ENDC}")
        
        return abs_file_list


    #submit ccCal jobs
    def submit_ccCal(self):
        HKL_list = self.getHKLList()
        if os.path.isdir(self.realTimeUpdate.shareWorkDir) and HKL_list:
            ccList(HKL_list, self.realTimeUpdate.shareWorkDir)
        elif not os.path.isdir(self.realTimeUpdate.shareWorkDir):
            print(f"{colors.RED}Working dir does not exist, please check: {self.realTimeUpdate.shareWorkDir}{colors.ENDC}")
        elif not HKL_list:
            print(f"{colors.RED}No HKL file list, please check HKL paths{colors.ENDC}")
        else:
            print(f"{colors.RED}Unknow problem with ccCal input HKLs or working dir, please check{colors.ENDC}")
        
        #check ccClusterLog.txt after generation
        self.check_ccCalLogStatus()
        self.realTimeUpdate.Tab_ccCluster.auto_select_ccClusterLogPath()        


    #select output folder
    def select_WorkDir(self):
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
            self.ccCallog_status.setStyleSheet("color: green; font-weight: bold; font-size: 14px")
        else:
            self.ccCallog_status.setText(f"ccClusterLog.txt not found in {self.realTimeUpdate.shareWorkDir}. Please generate one or select a different path")
            self.ccCallog_status.setStyleSheet("color: red; font-weight: bold; font-size: 14px")


    #search HKL files in the selected folder with the given suffix
    def search_HKL_files(self):
        matching_files = []
        FileName = self.FileName.text().strip()
        target_folder_text = self.SearchDir.text().strip()
        if not target_folder_text:
            print(f"{colors.RED}Please put a target folder for searching HKL files{colors.ENDC}")
            return
        
        if not FileName:
            print(f"{colors.RED}Please enter a file name or pattern to search{colors.ENDC}")
            return

        folder_paths = glob.glob(target_folder_text)
        if not folder_paths:
            print(f"{colors.RED}No folder found: {target_folder_text}{colors.ENDC}")
            return

        for folder in folder_paths:
            if os.path.isdir(folder):
                for file_path in Path(folder).glob(f"**/{FileName}"):
                    if file_path.is_file():
                        abs_path = str(file_path.absolute())
                        matching_files.append(abs_path)
                        print(f"{colors.BLUE}Found: {abs_path}{colors.ENDC}")
            else:
                print(f"{colors.RED}Folder does not exist: {folder}{colors.ENDC}")
        
        if matching_files:
            matching_files.sort()
            self.HKLPaths_text.clear()
            self.HKLPaths_text.setPlainText("\n".join(matching_files))
            print(f"{colors.GREEN}Added {len(matching_files)} file(s){colors.ENDC}")
        else:
            self.HKLPaths_text.clear()
            print(f"{colors.RED}No files found: {FileName}{colors.ENDC}")



#define Class for the tab for ccClustering
#tabs for result and summary are generated 
#through different modules
class tab_ccCluster(QWidget):
    def __init__(self, realTimeUpdates:MainWindow, **kwargs):
        super().__init__()

        #get real time update
        self.realTimeUpdate = realTimeUpdates

        #vertical layout
        self.ccCluster_layout = QtWidgets.QVBoxLayout(self)

        #set up buttons
        self.ccClusterSetup_area()
        #set up plot area
        self.plotDendroAndStatistic_area()

        #separation line:
        Septline = QtWidgets.QFrame()
        Septline.setFrameShape(QtWidgets.QFrame.HLine)
        Septline.setStyleSheet("background-color: #888888; max-height: 3px; border: none; margin: 5px 0;")

        #Add widget to the layout
        self.ccCluster_layout.addWidget(self.ccClusterSetup_widget, 1)
        self.ccCluster_layout.addWidget(Septline)
        self.ccCluster_layout.addWidget(self.plotDendroAndStatisticWidget, 9)

        #setup tabs if result exists
        self.UpdateResultAndSyncTabs()



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
        self.ShowThreshold.setPlaceholderText("Put threshold or click Auto Threshold")
        self.AutoThreshols = QtWidgets.QPushButton("Auto Threshold")
        self.AutoThreshols.clicked.connect(self.getAutoThreshold)

        #button to show the largest group number with current threshold
        self.ShowLargestGroup = QtWidgets.QLineEdit()
        self.ShowLargestGroup.setPlaceholderText("Show largest group with current threshold")
        self.CheckLargestGroup = QtWidgets.QPushButton("Check Largest Group")
        self.CheckLargestGroup.clicked.connect(self.getLargestGroup)

        #define which group to merge based on the pre-view Dendrogram
        self.mergeGroupLabel = QtWidgets.QLabel("Selected cluster:")
        self.mergeGroupLabel.setStyleSheet("color: black; font-weight: bold;;")

        self.mergeGroup = QtWidgets.QLineEdit()
        self.mergeGroup.setPlaceholderText("Put the cluster to merge with space or comma, e.g. 2 3, 5,8")

        #select reference HKL file for XSCALE merging
        self.reference_HKL = QtWidgets.QLineEdit()
        self.reference_HKL.setPlaceholderText("Put absolute path of reference HKL file for XSCALE merging, optional")
        self.select_reference_HKL = QtWidgets.QPushButton("Select reference HKL")
        self.select_reference_HKL.clicked.connect(self.select_reference_HKL_file)

        #anomalous flag
        self.anomBox = QtWidgets.QCheckBox("Anomalous data")
        self.anomBox.setChecked(False)

        #status bar to show information and button to run ccCluster job
        self.ccClusterStatusBarTitle = QtWidgets.QLabel("Monitor ccCluster job status and button to run the job")
        self.ccClusterStatusBarTitle.setStyleSheet("color: black; font-weight: bold; font-size: 12px;")
        self.ccClusterStatusBar = QtWidgets.QLineEdit()
        self.ccClusterStatusBar.setReadOnly(True)
        self.update_ccClusterStatusBar("ready to work")
        self.RunccCluster = QtWidgets.QPushButton("Run ccCluster")
        self.RunccCluster.setFixedSize(200, 50)
        self.RunccCluster.clicked.connect(self.submit_ccCluster)

        #create main layout
        ccClusterSetupLayout = QtWidgets.QVBoxLayout(self.ccClusterSetup_widget)

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
        AutoThresholsLayout.addWidget(self.mergeGroupLabel, 1)
        AutoThresholsLayout.addWidget(self.mergeGroup, 3)

        #put reference HKL and button together
        referenceHKLLayout = QtWidgets.QHBoxLayout()
        referenceHKLLayout.addWidget(self.reference_HKL, 4)
        referenceHKLLayout.addWidget(self.select_reference_HKL, 1)

        #put status bar and merge button together
        StatusBarLayout = QtWidgets.QVBoxLayout()
        StatusBarLayout.addWidget(self.ccClusterStatusBarTitle, alignment=QtCore.Qt.AlignCenter)
        StatusBarLayout.addWidget(self.ccClusterStatusBar)
        StatusBarLayout.addWidget(self.RunccCluster, alignment=QtCore.Qt.AlignCenter)

        #setup title
        self.ccClusterSetupTitle = QtWidgets.QLabel("Set and run ccCluster jobs")
        self.ccClusterSetupTitle.setStyleSheet("color: black; font-weight: bold; font-size: 16px;")

        #pack everything in the layout
        ccClusterSetupLayout.addWidget(self.ccClusterSetupTitle, alignment=QtCore.Qt.AlignCenter)
        ccClusterSetupLayout.addLayout(ccClusterLogLayout)
        ccClusterSetupLayout.addLayout(AutoThresholsLayout)
        ccClusterSetupLayout.addLayout(referenceHKLLayout)
        ccClusterSetupLayout.addLayout(StatusBarLayout)


    #select reference HKL file for XSCALE merging
    def select_reference_HKL_file(self):
        file_path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select reference HKL File", "", "HKL Files (*.HKL *.hkl);;All Files (*)")

        if file_path:
            self.reference_HKL.setText(os.path.abspath(file_path))


    #plot the Dendrogram and SCALE.LP statistics in different sub-tabs
    def plotDendroAndStatistic_area(self):
        #define widget
        self.plotDendroAndStatisticWidget = QtWidgets.QWidget()

        #Create layout
        self.plotDendroAndStatisticLayout = QtWidgets.QVBoxLayout(self.plotDendroAndStatisticWidget)

        #add title to the top of widget
        self.titleWidget = QtWidgets.QLabel("Show result Dendrogram and statistics for each merged result")
        #set up font
        self.titleWidget.setStyleSheet("color: black; font-weight: bold; font-size: 16px;")

        #Create tabs for Dendrogram and statistics
        self.DendroAndStatsPlot()

        #pack widgets into the layout
        self.plotDendroAndStatisticLayout.addWidget(self.titleWidget, 1, alignment=QtCore.Qt.AlignCenter)
        self.plotDendroAndStatisticLayout.addWidget(self.ResultDendroAndStatsTab, 19)


    #Add ccCluster log path is needed
    def select_ccClusterLogPath(self):
        file_path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select ccCluster log File", "", "ccCluster Log Files (*.txt);;All Files (*)")

        if file_path:
            self.ccClusterLogPath_text.setText(os.path.abspath(file_path))


    #Auto select ccCluster log file if exists in work Dir
    def auto_select_ccClusterLogPath(self):
        ccClusterLog = os.path.join(self.realTimeUpdate.shareWorkDir, "ccClusterLog.txt")
        print(f"ccClusterlog: {ccClusterLog}")

        #update the ccClusterLog.txt status
        if os.path.isfile(ccClusterLog):
            self.ccClusterLogPath_text.setText(ccClusterLog)
        else:
            self.ccClusterLogPath_text.setText(f"ccClusterLog.txt not found in {self.realTimeUpdate.shareWorkDir}. Please generate one or select a different path")


    #update ccCluster bar:
    def update_ccClusterStatusBar(self, status:str):
        print(f"{status}")
        self.ccClusterStatusBar.setText(f"{status}")
        #self.ccCallog_status.setStyleSheet("color: green; font-weight: bold")

    
    #auto define threshold
    def getAutoThreshold(self):
        CC, _, _, _, status_text = self.realTimeUpdate.setupCC(self.ccClusterLogPath_text.text())
        if CC is None:
            self.update_ccClusterStatusBar(status_text)
        else:
            Threshold = CC.thrEstimation()
            self.ShowThreshold.setText(str(Threshold))
            GroupNum, largestGroup, totalHKL = CC.checkMultiplicity(Threshold)
            self.getLargestGroup()


    #Show largest cluster
    def getLargestGroup(self):
        CC, _, _, _, status_text= self.realTimeUpdate.setupCC(self.ccClusterLogPath_text.text())
        if CC is None:
            self.update_ccClusterStatusBar(status_text)
        else:
            #remove white space and check if the threshold is a valid number
            ThresholdValue = self.ShowThreshold.text().strip()
            if ThresholdValue:
                try:
                    threshold_val = round(float(ThresholdValue), 2)
                    GroupNum, largestGroup, totalHKL = CC.checkMultiplicity(threshold_val)
                    self.update_ccClusterStatusBar(
                                                    f"Auto Threshold is {threshold_val}, "
                                                    f"the largest cluster number is {GroupNum} "
                                                    f"with {largestGroup}/{totalHKL} files"
                                                    )
                    self.ShowLargestGroup.setText(f"Largest Group: {GroupNum}; HKLs: {largestGroup}/{totalHKL}")
                except ValueError:
                    self.update_ccClusterStatusBar("Threshold must be a valid number (e.g., 0.5)")
            else:
                self.update_ccClusterStatusBar("Please input a threshold value")


    #submit ccCluster job
    def submit_ccCluster(self):
        CC, Tree, etiquets, _, status_text = self.realTimeUpdate.setupCC(self.ccClusterLogPath_text.text())
        threshold = self.ShowThreshold.text()

        #pass selected group as a str with space:
        SelectClusterText = self.mergeGroup.text().replace(',', ' ').strip()

        #check if empty, or try converting to integers
        if not SelectClusterText:
            SelectedCluster = None
        else:
            try:
                SelectedCluster = [int(x) for x in SelectClusterText.split()]
            except ValueError:
                SelectedCluster = None

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
                if not self.reference_HKL.text():
                    xscale_checker, xscale_path = CC.prepareXSCALE(anomlous, threshold, clusterList=SelectedCluster)
                elif os.path.isfile(self.reference_HKL.text()):
                    xscale_checker, xscale_path = CC.prepareXSCALE(anomlous, threshold, clusterList=SelectedCluster, refHKL=self.reference_HKL.text())
                else:
                    xscale_checker, xscale_path = CC.prepareXSCALE(anomlous, threshold, clusterList=SelectedCluster)
                if xscale_checker == True:
                    self.update_ccClusterStatusBar(f"Running XSCALE job in {xscale_path}")
                    CC.scaleAndMerge(anomlous, threshold, xscale_path)
                    #get json from XSCALE
                    CC.flatClusterPrinter(threshold, etiquets, anomlous, xscale_path)

                #prepare and run Pointless
                pointless_checker, pointless_path = CC.preparePointless(anomlous, threshold, clusterList=SelectedCluster)
                if pointless_checker == True:
                    self.update_ccClusterStatusBar(f"Running Pointless job in {pointless_path}")
                    CC.pointlessRun(anomlous, threshold, pointless_path)
                    #prepare and run aimless
                    self.update_ccClusterStatusBar(f"Running Aimless job in {pointless_path}, please be patient")
                    CC.aimlessRun(anomlous, threshold, pointless_path)

                #update result
                self.UpdateResultAndSyncTabs()
                self.ccClusterStatusBar.setText(f"ccCluster job finished, please check the result in {pointless_path}")

            #CC.passOInfoToGA(threshold, etiquets, anomlous)
            elif fileType=="mtz":
                #prepare and run Pointless
                pointless_checker, pointless_path = CC.preparePointless(anomlous, threshold, clusterList=SelectedCluster)
                if pointless_checker == True:
                    self.update_ccClusterStatusBar(f"Running Pointless job in {pointless_path}")
                    CC.pointlessRun(anomlous, threshold, pointless_path)
                    #prepare and run aimless
                    self.update_ccClusterStatusBar(f"Running Aimless job in {pointless_path}, please be patient")
                    CC.aimlessRun(anomlous, threshold, pointless_path)
                    CC.flatClusterPrinter(threshold, etiquets, anomlous, pointless_path)

                self.ccClusterStatusBar.setText(f"No statistcs as the input file is mtz, check results in: {pointless_path}")
            else:
                self.update_ccClusterStatusBar(f"Unknown input file format, please check ccCluster log file: {self.ccClusterLogPath_text.text()}")


    def UpdateResultAndSyncTabs(self):
        self.realTimeUpdate.CheckAndShowResult()
        self.SyncResultTabs()
        self.realTimeUpdate.Tab_plotStats.syncResultList()


    #prepare DendroGram tab from png for the result folder
    def DendrogramFromPNG(self, dendrogram_path:str):
        #set up the Dendrogram plot tab
        DendroprocessedWidget = QtWidgets.QWidget()
        Dendrolayout = QtWidgets.QVBoxLayout(DendroprocessedWidget)

            #Load and scale pixmap otherwise the image is TOO big
        pixmap = QtGui.QPixmap(dendrogram_path)
        
        #Scale to fit within max size
        ImageBox = QtWidgets.QLabel(self)
        ImageBox.setPixmap(pixmap)
        ImageBox.setScaledContents(True)
        ImageBox.setMinimumSize(300, 200)

        Dendrolayout.addWidget(ImageBox)

        return DendroprocessedWidget


    #prepare tabs for XSCALE statistics
    def CreateXSCALEStatTab(self, XSCALEFile:str):
        XSCALEstat = QtWidgets.QWidget()
        XSCALEstatlayout = QtWidgets.QVBoxLayout(XSCALEstat)

        #create a read-only text edit to show the XSCALE.LP statistics
        XSCALEText = QtWidgets.QTextEdit()
        XSCALEText.setReadOnly(True)

        #Get the content
        _, plotText = extractXSCALEStat(XSCALEFile)
        XSCALEText.setText(plotText)

        #put in to layout
        XSCALEstatlayout.addWidget(XSCALEText)

        #return the widget
        return XSCALEstat


    #function to add result tab in the self.PlottingTabWidget for realtime update when new result is generated or deleted
    def CreateResultTabs(self, result_name:str):
        result_folder = os.path.join(os.path.abspath(self.realTimeUpdate.shareWorkDir), result_name)
        dendrogram_path = os.path.join(f"{result_folder}/gallery", "Dendrogram.png")
        xscale_path = os.path.join(result_folder, "XSCALE.LP")

        #check if the dendrogram.png and XSCALE.LP exist in the result folder
        if not os.path.isfile(dendrogram_path):
            self.update_ccClusterStatusBar(f"No Dendrogram.png found in {result_folder}, please check")
            return

        if not os.path.isfile(xscale_path):
            self.update_ccClusterStatusBar(f"No XSCALE.LP found in {result_folder}, please check")
            return

        #set up the tab for Dendrogram and statistics
        resultTab = QtWidgets.QTabWidget()

        #set up the Dendrogram plot tab
        DendroprocessedWidget = self.DendrogramFromPNG(dendrogram_path)

        #add dendrogram plot for the result as a tab
        resultTab.addTab(DendroprocessedWidget, "Dendrogram")

        #set up the tab for statistics from XSCALE.LP
        resultTab.addTab(self.CreateXSCALEStatTab(xscale_path), "XSCALE Statistics")

        #plot The statistics from XSCALE.LP
        resultTab.addTab(SinglePlotTab(result_folder), "XSCALE Statistics Plot")

        #Return tabs for the result, will be added to the self.PlottingTabWidget
        return resultTab


    #remove non-existing result tabs from the self.PlottingTabWidget for realtime update when new result is generated or deleted
    def RemoveResultTabs(self, result_name:str):
        #check if self.PlottingTabWidget exists
        if not hasattr(self, 'PlottingTabWidget') or self.PlottingTabWidget is None:
            return

        #Find the tabs and remove it saftely
        for i in range(self.PlottingTabWidget.count()):
            if self.PlottingTabWidget.tabText(i) == result_name:
                widget = self.PlottingTabWidget.widget(i)
                #remove the tab from the PlottingTabWidget
                self.PlottingTabWidget.removeTab(i)
                if widget:
                    #delete the widget to free memory
                    widget.deleteLater()
                print(f"Removed tab: {result_name}")
                return


    #sync the result tabs with the self.MergeResult list, remove non-existing result tabs from the self.PlottingTabWidget
    def SyncResultTabs(self):
        print(f"Syncing result tabs with MergeResult list: {self.realTimeUpdate.MergeResult}")
        #check if self.PlottingTabWidget exists
        if not hasattr(self, 'PlottingTabWidget') or self.PlottingTabWidget is None:
            print(f"{colors.RED}PlottingTabWidget does not exist, cannot sync result tabs{colors.ENDC}")
            return

        #update the Tab list
        tablist = []
        for i in range(self.PlottingTabWidget.count()):
            tab_text = self.PlottingTabWidget.tabText(i)
            print(f"Existing tab: {tab_text}")
            # Skip the pre-plot tab
            if tab_text != "Pre-plot Dendrogram with threshold":
                tablist.append(tab_text)

        #get list of tabs to add and remove
        tabs_to_add = [f for f in self.realTimeUpdate.MergeResult if f not in tablist]
        tabs_to_remove = [f for f in tablist if f not in self.realTimeUpdate.MergeResult]

        #add ResultTabs
        for folder_name in tabs_to_add:
            resultPlotTab  = self.CreateResultTabs(folder_name)
            if resultPlotTab is not None:
                self.update_ccClusterStatusBar(f"Adding result tab for {folder_name}")
                self.PlottingTabWidget.addTab(resultPlotTab, folder_name)
                self.update_ccClusterStatusBar(f"Added result tab for {folder_name}")

        #remove ResultTabs
        for folder_name in tabs_to_remove:
            self.update_ccClusterStatusBar(f"Removing result tab for {folder_name}")
            self.RemoveResultTabs(folder_name)
            self.update_ccClusterStatusBar(f"Removed result tab for {folder_name}")

        #update log
        status_msg = f"Synced tabs: +{len(tabs_to_add)} added, -{len(tabs_to_remove)} removed"
        self.update_ccClusterStatusBar(status_msg)
        print(f"workdir: {self.realTimeUpdate.shareWorkDir}")


    #create tabs in the Result tab to show dendrogram and statistics for each merged result
    def DendroAndStatsPlot(self):            
        #set up the widget and layout
        self.ResultDendroAndStatsTab = QtWidgets.QWidget()
        self.ResultDendroAndStatsTabLayout = QtWidgets.QVBoxLayout(self.ResultDendroAndStatsTab)
  
        #plotting widget make it self so we can update tabs in real time when new result is generated/deleted
        self.PlottingTabWidget = QtWidgets.QTabWidget()

        #add the dendro pre plot tab to the plotting layout
        self.PlottingTabWidget.addTab(PrePlotDendrogram(self.ccClusterLogPath_text, self.ShowThreshold, self.realTimeUpdate.setupCC), "Pre-plot Dendrogram with threshold")
        
        #setup tabs for each merged result
        if not self.realTimeUpdate.MergeResult:
            self.update_ccClusterStatusBar(f"No merged result found in {self.realTimeUpdate.shareWorkDir}, please check")
        else:
            self.SyncResultTabs()

        #add tab to the main widget
        self.ResultDendroAndStatsTabLayout.addWidget(self.PlottingTabWidget)



class tab_plotStats(QtWidgets.QWidget):
    def __init__(self, realTimeUpdates:MainWindow, **kwargs):
        super().__init__()
        #get real time update from self of MainWindow
        self.realTimeUpdate = realTimeUpdates

        #selected file list for plotting, will be updated when new result is generated or deleted
        self.selectedList = []

        #vertical layout
        self.plotStats_layout = QtWidgets.QVBoxLayout(self)

        #set up buttons
        self.ResultSelection_area()
        #set up plot area
        self.CompareResults_area()

        #separation line:
        Septline = QtWidgets.QFrame()
        Septline.setFrameShape(QtWidgets.QFrame.HLine)
        Septline.setStyleSheet("background-color: #888888; max-height: 3px; border: none; margin: 5px 0;")

        #Add widget to the layout
        self.plotStats_layout.addWidget(self.ResultSelectionWidget, 1)
        self.plotStats_layout.addWidget(Septline)
        self.plotStats_layout.addWidget(self.CompareResultsWidget, 4)


    def ResultSelection_area(self):      
        self.ResultSelectionWidget = QtWidgets.QWidget()
        self.ResultSelectionLayout = QtWidgets.QVBoxLayout(self.ResultSelectionWidget)

        #set up selection box
        self.ResultSelectionBox = QtWidgets.QListWidget()
        self.ResultSelectionBox.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.ResultSelectionBox.setMinimumHeight(150)
        self.ResultSelectionBox.setMaximumHeight(300)

        #set up scroll bar for the list widget
        self.ResultSelectionBox.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.ResultSelectionBox.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)

        #connect the selection change signal to the plot function
        self.ResultSelectionBox.itemSelectionChanged.connect(self.ReadAndPlotSelectedResult)

        #setup title
        self.Resultselecttitle = QtWidgets.QLabel("Please select result in the box to compare them")
        self.Resultselecttitle.setStyleSheet("color: black; font-weight: bold; font-size: 16px;")

        #pack layout
        self.ResultSelectionLayout.addWidget(self.Resultselecttitle, 1, alignment=QtCore.Qt.AlignCenter)
        self.ResultSelectionLayout.addWidget(self.ResultSelectionBox, 19)



    def CompareResults_area(self):
        #create a widget to hold the plot area
        self.CompareResultsWidget = QtWidgets.QWidget()
        self.CompareResultsLayout = QtWidgets.QVBoxLayout(self.CompareResultsWidget)
        self.ComparetitleWidget = QtWidgets.QLabel("Compare statistics for selected merged result")
        self.ComparetitleWidget.setStyleSheet("color: black; font-weight: bold; font-size: 16px;")

        #create a widget to hold the plot area for XSCALE statistics comparison
        self.XSCALECompareWidget = QtWidgets.QTabWidget()
        self.XSCALEComparelayout = QtWidgets.QVBoxLayout(self.XSCALECompareWidget)

        XSCALECompareplaceholder = QtWidgets.QLabel("Select folders from the list above to compare XSCALE statistics")
        XSCALECompareplaceholder.setAlignment(QtCore.Qt.AlignCenter)
        self.XSCALECompareWidget.addTab(XSCALECompareplaceholder, "No Selection")

        self.CompareResultsLayout.addWidget(self.ComparetitleWidget, 1, alignment=QtCore.Qt.AlignCenter)
        self.CompareResultsLayout.addWidget(self.XSCALECompareWidget, 19)


    #sync the result list in the selection widget with the result list in the ccCluster tab
    #pass this function to the ccCluster tab to update the result list in real time when new result is generated or deleted
    #thus we need to have Tab_plotStats called before Tab_ccCluster in the MainWindow class
    def syncResultList(self):
        print(f"Syncing result list with MergeResult: {self.realTimeUpdate.MergeResult}")
        
        if not hasattr(self, 'ResultSelectionBox') or self.ResultSelectionBox is None:
            print(f"{colors.RED}ResultSelectionBox does not exist{colors.ENDC}")
            return
        
        #Save current selection
        current_selection = []
        for selecteditem in self.ResultSelectionBox.selectedItems():
            current_selection.append(selecteditem.text())
        print(f"Current selection preserved: {current_selection}")
        
        # Clear and repopulate
        self.ResultSelectionBox.clear()
        for resultName in self.realTimeUpdate.MergeResult:
            self.ResultSelectionBox.addItem(resultName)
        
        #Restore selection
        for i in range(self.ResultSelectionBox.count()):
            item = self.ResultSelectionBox.item(i)
            if item.text() in current_selection:
                item.setSelected(True)
                print(f"Restored selection: {item.text()}")
        
        # Update status
        count = self.ResultSelectionBox.count()
        selected_count = len(self.ResultSelectionBox.selectedItems())
        print(f"Result list updated: {count} folders, {selected_count} selected")


    #update the result list in the selection widget when new result is generated or deleted
    def UpdateResultList(self):
        self.selectedList = []
    
        for item in self.ResultSelectionBox.selectedItems():
            self.selectedList.append(item.text())

        print(f"Updated selectedList: {self.selectedList}")


    #Plot the selected result folders in the selection widget it will be one tab for refresh the result list and plot the selected result folders
    def PlotSelectedResult(self):
        #remove tabs:
        while self.XSCALECompareWidget.count() > 0:
            self.XSCALECompareWidget.removeTab(0)

        #add tabs for the selected result folders
        if self.selectedList:
            abs_ResultDirList = [f"{os.path.abspath(self.realTimeUpdate.shareWorkDir)}/{i}" for i in self.selectedList if os.path.isfile(f"{os.path.abspath(self.realTimeUpdate.shareWorkDir)}/{i}/XSCALE.LP")]
            self.XSCALECompareWidget.addTab(MultiPlotTab(abs_ResultDirList), f"Plotted {len(self.selectedList)} folder(s)")
            print(f"Plotted {len(self.selectedList)} folder(s): {self.selectedList}")
        else:
            placeholder = QtWidgets.QLabel("No folders selected. Please select at least one folder.")
            placeholder.setAlignment(QtCore.Qt.AlignCenter)
            self.XSCALECompareWidget.addTab(placeholder, "No Selection")
            print("No result folder selected")


    #read updated selectedList from the selection widget and plot the selected result folders
    def ReadAndPlotSelectedResult(self):
        self.UpdateResultList()
        self.PlotSelectedResult()



#put the tabs to gether in one GUI
class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, **kwargs):
        super().__init__()
        #make one or several linEdit or textEdit accessable to all tabs in realtime
        #set initial work dir as pwd and pass it to shared place
        self.WorkDir = os.getcwd()
        self.updateWorkDir(self.WorkDir)

        #put the result list in the self of MainWindow, so it can be accessed by all tabs
        self.MergeResult = []
        self.CheckAndShowResult()
        #self.Tab_ccCluster.SyncResultTabs()

        #pass the self of MainWindow to the tabs as argument
        #to enalbe realtime update on lineEdit
        #DO NOT change the order of the tabs, as the tab_ccCluster is used in tab_plotStats to get the result list
        self.Tab_plotStats = tab_plotStats(self)
        self.Tab_ccCluster = tab_ccCluster(self)
        self.Tab_ccCal = tab_ccCal(self)
        

        #add tabs
        self.tabWidget = QtWidgets.QTabWidget()
        self.setCentralWidget(self.tabWidget)
        self.tabWidget.addTab(self.Tab_ccCal, "ccCal tab")
        self.tabWidget.addTab(self.Tab_ccCluster, "ccCluster tab")
        self.tabWidget.addTab(self.Tab_plotStats, "Plot statistics tab")

        #Set up main window
        self.setObjectName("MainWindow")
        self.setWindowTitle("Cluster and merge")
        self.resize(1600, 1200)
    
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
        etiquets, etiquetsWithNum = CC.createLabels()
        text = f"CC setup successful: {correlationFile}"
        return CC, Tree, etiquets, etiquetsWithNum, text

        #Update results folder list, will be used by the result compare tab
    def CheckAndShowResult(self):
        abs_FolderPaths = Path(os.path.abspath(self.shareWorkDir))
        for folder_path in abs_FolderPaths.glob("cc_Cluster_*"):
            if folder_path.is_dir():
                if not (folder_path/"XSCALE.LP").is_file():
                    print(f"No XSCALE.LP found in {folder_path}, please check")
                    print(f"No XSCALE.LP found in {folder_path}, please check")
                    continue

                if not (folder_path/"gallery/Dendrogram.png").is_file():
                    print(f"No Dendrogram.png found in {folder_path}, please check")
                    print(f"No Dendrogram.png found in {folder_path}, please check")
                    continue

                if (folder_path/"XSCALE.LP").is_file() and (folder_path/"gallery/Dendrogram.png").is_file():
                    folder_name = folder_path.name
                    if folder_name not in self.MergeResult:
                        print(f"find result folder: {folder_path}")
                        self.MergeResult.append(folder_name)
            else:
                print(f"{colors.RED}Folder path does not exist: {folder_path}{colors.ENDC}")

        #remove the not existing result folder from the list
        Exist_results = []
        for result_folder in self.MergeResult:
            print(f"result folder name: {result_folder}")
            abs_result_folder = os.path.join(os.path.abspath(self.shareWorkDir), result_folder)
            if not os.path.isdir(abs_result_folder):
                print(f"Result folder {abs_result_folder} does not exist, remove it from the list")
            elif not os.path.isfile(os.path.join(abs_result_folder, "XSCALE.LP")):
                print(f"No XSCALE.LP found in {abs_result_folder}, remove it from the list")
            elif not os.path.isfile(os.path.join(f"{abs_result_folder}/gallery", "Dendrogram.png")):
                print(f"No Dendrogram.png found in {abs_result_folder}, remove it from the list")
            else:
                Exist_results.append(result_folder)

        self.MergeResult = Exist_results

        self.MergeResult.sort()



#prepare run the GUI
def main():
    app = QApplication(sys.argv)
    ex = MainWindow()
    ex.show()
    sys.exit(app.exec_())      


#run the GUI
if __name__== '__main__':
    main()
