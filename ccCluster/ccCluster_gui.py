#! /usr/bin/env python3
from __future__ import print_function, absolute_import

__author__ = "Gianluca Santoni"
__copyright__ = "Copyright 20150-2019"
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

#import CalcClass
from scipy.cluster import hierarchy
import collections
import operator
from time import sleep
import os
from .resultsTab import resultsTab
from .summary import resultsSummary
from .clustering import Clustering

#Insert parse  to change the file path from command line
import argparse

#Deal with wild card
import glob
from pathlib import Path'
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



#set up the parameters, maybe not necessary as you should be able to add them via the GUI
def process_args():
    input_args = argparse.ArgumentParser()
    input_args.add_argument("-i","--DISTfile", 
    default=None, 
    help="Distance file from ccCalc module"
    )
    
    #input_args.add_argument("-f", dest="structures", default= None ,  type=str, nargs='+', help='The list of refined structures to merge')
    #input_args.add_argument("-o","--outname", dest="outname", default='Dendrogram', help="output dendogram file name")
    
    input_args.add_argument("-t", "--threshold", 
    dest="threshold", 
    help="Distance threshold for clustering"
    )

    input_args.add_argument("-c", "--count",
    action="store_true", 
    dest="count", 
    default=False, 
    help="Counts datasets in the biggest cluster and exit"
    )

    input_args.add_argument("-e", "--estimation",
    action="store_true", 
    dest="est", 
    default=False, 
    help="Tries to guess an optimal threshold value"
    )

    input_args.add_argument("-wd", "--work_dir",
    type = str,
    help = "output directory, default is pwd",
    )
    
    #input_args.add_argument("-u", dest="cell", default= False , action="store_true" , help='Unit cell based clustering. requires list of input files')

    #input_args.print_help()
    args= input_args.parse_args()

    #set work dir to pwd if none
    if args.output_dir == None:
        args.output_dir = os.getcwd()


    #Suggest to run ccCalc if no correlation file is provided
    #Call to ccCalc if no distances found but files listed
    if args.DISTfile is None: 
        print('no inputs specified, please run ccCalc before')
    else:
        correlationFile=args.DISTfile
    
    return args



#here is the class that creates the tab for ccCal, 
#to generate ccClusterLog.txt in the target folder if it is not exist
class tab_ccCal(QWidget):
    def __init__(self, realTimeUpdates:MainWindow, **kwargs):
        #pass self to parent
        super().__init__()

        #get real time update
        self.realTimeUpdate = realTimeUpdates

        #set initial work dir as pwd
        self.WorkDir = os.getcwd()

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
        #update changes in the WorkDir entry in real time

        self.WorkDir_entry.textChanged.connect(self.updateWorkDir)
        #check ccClusterlog.txt when change the work dir:
        self.WorkDir_entry.textChanged.connect(self.check_ccCalLogStatus)

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
            print(f"{colors.RED}No HKL file list, please check HKL paths{color.ENDC}")
        else:
            print(f"{colors.RED}Unknow problem with ccCal input HKLs or working dir, please check{color.ENDC}")
        
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
    def check_ccCalLogStatus():
        ccClusterLog = os.path.join(self.workDir.text(), "ccClusterLog.txt")

        #update the ccClusterLog.txt status
        if os.path.isfile(ccClusterLog):
            self.ccCallog_status.setText(f"ccClusterLog.txt found: {ccClusterLog}")
            self.ccCallog_status.setStyleSheet("color: green; font-weight: bold")
        else:
            self.ccCallog_status.setText(f"ccClusterLog.txt not found in {self.workDir.text()}. Please generate one or select a diffenent path")
            self.ccCallog_status.setStyleSheet("color: red; font-weight: bold")


    #update workDir
    def updateWorkDir(self, text):
        self.realTimeUpdate.shareWorkDir = text



#define Class for the tab for ccClustering
#tabs for result and summary are generated 
#through different modules
class tab_ccCluster(QWidget):
    def __init__(self, realTimeUpdates:MainWindow, **kwargs):
        super()__init__()
        self.main_canvas = QHBoxLayout(self)

        #set up buttons
        self.ccCluster_buttons()
        #set up plot area
        self.ccCluster_plot_area()

        #other parameters
        self.counter =1


    #set up the buttons for ccCluster
    def ccCluster_buttons(self):
        #buttons create Dendrogram
        self.PlotButton = QtWidgets.QPushButton(self.verticalLayoutWidget)
        self.PlotButton.setObjectName(_fromUtf8("PlotButton"))
        self.PlotButton.clicked.connect(self.createDendrogram)

        #button to run ccCluster

        #button to show auto threshold result

        #button to show the largest group number with current threshold

    def setupUi(self, MainWindow):
        global Tree
        global threshold
        global etiquets
        self.counter = 1
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.setWindowTitle("Cluster and merge")
        MainWindow.resize(1215, 1000)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.processedValues = []
        self.threshold = threshold
        self.CurrentDir = os.getcwd()
        self.etiquets= etiquets


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

        ##########
        #Show previous results
        # to be removed from the setupUi class
        ###########
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


    def showSummary(self):
        self.sum = resultsSummary()
        self.sum.show()


    def createDendrogram(self):
        X = hierarchy.dendrogram(Tree, color_threshold=self.threshold)
        #self.textOutput.append('Plotted Dendrogram. Colored at a %s threshold for distance'%(threshold))
        self.TreeCanvas.draw()


    def getThreshold(self,event):
        self.coord = 0
        if event.ydata != None:
            self.coord = event.ydata
            #self.textOutput.append('Threshold is %s'%(self.coord))
            #self.thr = event.ydata
            self.threshold= float('%.2f'%(event.ydata))
            self.createDendrogram()
            self.statusbar.showMessage('New threshold value: %.2f'%(self.threshold))


    def onChanged(self, text):
        self.threshold = float(text)


    def processClusters(self):
        Log = open(self.CurrentDir+'/.cc_cluster.log', 'a')

        Clusters = hierarchy.fcluster(Tree, self.threshold, criterion='distance')
        counter=collections.Counter(Clusters)
        Best = max(counter.items(), key=operator.itemgetter(1))[0]
        #Chose if process all or just biggest cluster

        if self.checkBox.isChecked():
            ToProcess = [Best]    
            self.statusbar.showMessage('Processing the best cluster. It contains %s datasets'%(counter[Best]))
        else:
            ToProcess = set(Clusters)
            for key in ToProcess:
                if counter[key]==1:
                    ToProcess = [x for x in ToProcess if x != key]

        #flag anomalus process or not
        if self.anomBox.isChecked():
            self.anomFlag= 'ano'
        else:
            self.anomFlag= 'no_ano'

        # Delete previous processings and create working directories
        for x in ToProcess:
            if [self.threshold,x, self.anomFlag] not in  self.alreadyDone:
                os.mkdir(self.CurrentDir+'/cc_Cluster_%.2f_%s_%s'%(self.threshold,x, self.anomFlag))
                Xscale=open(self.CurrentDir+'/cc_Cluster_%.2f_%s_%s/XSCALE.INP'%(self.threshold,x, self.anomFlag), 'a')
                Pointless=open(self.CurrentDir+'/cc_Cluster_%.2f_%s_%s/launch_pointless.sh'%(self.threshold,x,self.anomFlag ), 'a')
                print('OUTPUT_FILE=scaled.hkl',file=Xscale)
                print('MERGE= TRUE', file=Xscale)
                print('pointless hklout clustered.mtz << eof', file=Pointless)
                if self.anomBox.isChecked():
                    print('FRIEDEL\'S_LAW= FALSE', file=Xscale)
                Xscale.close()
                Pointless.close()

        for cluster, filename in zip(Clusters,self.etiquets):
            if cluster in ToProcess:
                OUT = open(self.CurrentDir+'/cc_Cluster_%.2f_%s_%s/XSCALE.INP'%(self.threshold,cluster,self.anomFlag), 'a')
                Pointless=open(self.CurrentDir+'/cc_Cluster_%.2f_%s_%s/launch_pointless.sh'%(self.threshold,cluster,self.anomFlag), 'a')
                print ('INPUT_FILE= ../%s'%(filename), file=OUT)
                #print ('INCLUDE_RESOLUTION_RANGE=20, 2', file=OUT)
                print ('MINIMUM_I/SIGMA= 0', file=OUT)
                print ('XDSIN ../%s'%(filename), file= Pointless)
                OUT.close()
                Pointless.close()
        #optional run of XSCALE
        
        for x in ToProcess:
            #newProcesses=[]
            if [self.threshold,x, self.anomFlag] not in  self.alreadyDone:
                self.statusbar.showMessage('XSCALE is processing cluster %.2f %s %s/'%(self.threshold, x,self.anomFlag))
                plt.savefig(self.CurrentDir+'/cc_Cluster_%.2f_%s_%s/Dendrogram.png'%(self.threshold,x,self.anomFlag))
                process = QtCore.QProcess(self)
                process.setWorkingDirectory(self.CurrentDir+'/cc_Cluster_%.2f_%s_%s/'%(self.threshold, x,self.anomFlag))
                process.start('xscale_par')
                print('Cluster, %s , %s , %s, %s'%(self.threshold,x, self.anomFlag, counter[x]), file=Log)             
                Pointless=open(self.CurrentDir+'/cc_Cluster_%.2f_%s_%s/launch_pointless.sh'%(self.threshold,x,self.anomFlag), 'a')
                print('COPY \n  TOLERANCE 4 \n eof', file= Pointless)
                Pointless.close()
                #newProcesses.append([self.threshold,x, self.anomFlag])
                L=[self.threshold,x, self.anomFlag, counter[x]]
            process.waitForFinished()
            process.exitStatus()
            sleep(0.5)
            st = os.stat(self.CurrentDir+'/cc_Cluster_%.2f_%s_%s/launch_pointless.sh'%(self.threshold,x,self.anomFlag ))
            os.chmod(self.CurrentDir+'/cc_Cluster_%.2f_%s_%s/launch_pointless.sh'%(self.threshold,x,self.anomFlag ), st.st_mode | 0o111)
            self.tabWidget.addTab(resultsTab(float(L[0]), L[1], L[2], L[3]), ('%.2f %s %s')%(float(L[0]), L[1], L[2]))
            self.alreadyDone.append([L[0], L[1], L[2]])
        self.statusbar.showMessage('Ready!')
        
        # for L in newProcesses:
        #     self.tabWidget.addTab(resultsTab(float(L[0]), L[1], L[2]), ('%.2f %s %s')%(float(L[0]), L[1], L[2]))
        #     self.alreadyDone.append([L[0], L[1], L[2]])



class tab_compare():
    def __init__(self, realTimeUpdates:MainWindow, **kwargs):
        super().__init__()
        
        #setup buttons


#put the tabs to gether in one GUI
class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        pass

        #make one or several linEdit or textEdit accessable to all tabs in realtime
        self.shareWorkDir = ""

        ###pass the self of MainWindow to the tabs as argument
        ###to enalbe realtime update on lineEdit
        self.Tab_ccCal = tab_ccCal(self)
        self.Tab_ccCluster = tab.ccCluster(self)
        self.Tab_compare = tab_compare(self)

        #add tabs
        self.tabWidget.addTab(self.Tab_ccCal, "ccCal tab")
        self.tabWidget.addTab(self.Tab_ccCluster, "ccCluster tab")
        self.tabWidget.addTab(self.Tab_compare, "Compare merging result tab")


#Main part of the program
#with the different options, we can chose 
# to process through the shell,
#count the multiplicity of the highest cluster
#run the interface
def main():
    args = process_args()


    #get the processing scripts and initial settings
    WorkDir = os.getcwd()
    CC = Clustering(correlationFile, WorkDir)
    Tree = CC.avgTree()
    etiquets=CC.createLabels()
    threshold = CC.thrEstimation()
    if args.threshold:
        threshold = args.threshold
    else:
        threshold= CC.thrEstimation()

    if args.count:
        CC.checkMultiplicity(threshold)
    elif args.est:
        a = CC.thrEstimation()
        print(a)
    else:
        app = QApplication(sys.argv)
        ex = MainWindow()
        ex.show()
        sys.exit(app.exec_())      

if __name__== '__main__':
    main()
