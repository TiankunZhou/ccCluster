from __future__ import print_function

__author__ = "Rita Giordano, Gianluca Santoni, Tiankun Zhou"
__copyright__ = "Copyright 2015-2026"
__credits__ = ["Rita Giordano, Gianluca Santoni, Tiankun Zhou, Alexander Popov"]
__license__ = ""
__version__ = "2.0"
__maintainer__ = "Tiankun Zhou"
__email__ = "tiankun.zhou@esrf.fr"
__status__ = "Beta"


from scipy.cluster import hierarchy
#import scipy
import matplotlib as mpl
from matplotlib.lines import Line2D
from matplotlib.figure import Figure
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
import numpy as np
import subprocess
import collections
import operator
import json
import random
import textwrap
import re
import matplotlib.patches as mpatches
from pathlib import Path


#Set color for printing
class colors:
    BOLD = '\033[1m'
    GREEN = '\033[32m'
    BLUE = '\033[34m'
    RED = '\033[31m'
    ENDC = '\033[m'


# find which cluster an index belongs to
def find_cluster_for_index(cluster_to_indices:dict, index:int):
    """Find which cluster an index belongs to"""
    for cluster, indices in cluster_to_indices.items():
        if index in indices:
            return cluster
    return None

#Read XSCALE.LP and extract information to plot the statistics
def extractXSCALEStat(XSCALEFile):
    plotList = []
    filteredList = []
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
        #Statistics table should have 14 columns, if not, print a warning and return
        if plotList:
            #check the length of the first line, if not 14, print a warning and return
            if len(plotList[0]) != 14:
                print(f"{colors.RED}Warning: The statistics table may not have 14 columns, it has {len(plotList[0])} columns\nUse this column number to align the table{colors.ENDC}")
                columnNum = len(plotList[0])
            else:
                columnNum = 14
            for line in plotList:
                if len(line) != columnNum:
                    print(f"{colors.RED}Warning: The line {line} does not have {columnNum} elements\nit has {len(line)} elements and will be skipped, merged dataset has issues at this resolution{colors.ENDC}")
                else:
                    filteredList.append(line)
            col_widths = [max(len(row[col_idx]) for row in filteredList) + 3 for col_idx in range(len(filteredList[0]))]
            for line in filteredList:
                aligned_line = "".join(f"{token:>{col_widths[i]}}" for i, token in enumerate(line))
                plotText += f"{aligned_line}\n"
    
    return filteredList, plotText


#function to check the indices if correct with the data number in the label list, to avoid any mismatch between the two
#only checks the FlatC list with the labelList ([[dataNum, filename], ...]) to see if the dataNum is correct with the index in the FlatC list
#make sure you pass the correct arguments to the function, otherwise it will not work as expected
def checkIndices(FlatCList:list, labelList:list):
    if len(FlatCList) != len(labelList):
        print(f"Warning: The length of FlatCList ({len(FlatCList)}) does not match the length of labelList ({len(labelList)}). Please check the input arguments.")
        return False

    for i, (cluster, label) in enumerate(zip(FlatCList, labelList)):
        dataNum = int(label[0])
        if i != dataNum:
            print(f"Warning: Mismatch at index {i}: FlatCList has cluster {cluster}, but labelList has data number {dataNum}.")
            return False
    return True


class Clustering():
    """
    parse cc_calc output and perform HCA
    at each call, it generates the distance matrix
    You get the dendrogram through Clustering.tree()
    """
    def __init__(self, ccCalcOutput, run_dir:str):
        self.ccFile= ccCalcOutput
        #self.CurrentDir = os.getcwd()
        self.RunDir = os.path.abspath(run_dir) # maybe better to use absolut path for better output control
        self.ccTable, self.Dimension = self.parseCCFile()
        self.createLabels()


    def parseCCFile(self):
        """
        Gets data from ccCalc ouput file and populates a numpy array with the distances
        """
        with open(self.ccFile, 'r') as f:
            dataArr = None
            data=[]
            Index = []
            #ignore everything before "Correlation coefficients"
            for line in f:
                if line.strip() == 'Correlation coefficients':
                    break

            for line in f:
                dataline= line.rstrip().split()
                data.append(dataline)
                Index.append(int(dataline[0])+1)
                Index.append(int(dataline[1])+1)
        Dimension=max(Index)
        dataArr = np.array(data,dtype=(float))

        return dataArr, Dimension


    def createLabels(self):
        """
        Gets the labels from the ccCalc output with the input file names
        """        
        self.labelList= []
        with open(self.ccFile) as f:   
            for line in f:
                if line.strip() == 'Labels':
                    break
            for line in f:
                if line.strip() == 'Correlation coefficients':
                    break
                goodLine = line.split()
                self.labelList.append(["%s"%(goodLine[1].strip('\n')), "%s"%(goodLine[2].strip('\n'))])

        #return an extra lable list with data number to double check the cluster umber is correct
        return self.labelList
        

    def inputType(self):
        """
        return input file type. Either mtz or HLK
        """        
        element = self.labelList[0][1]
        extension = element.split('.')[-1]
        print(extension)
        return extension


    def tree(self):
        """
        Returns the HCA dendrogrm, using the complete linkage method
        """        
        data = self.ccTable
        Matrix=np.zeros((self.Dimension,self.Dimension))

        reducedArray=[]
        for line in data:
                #print line
            if line is not None and len(line) != 0:
                 Matrix[int(line[0]),int(line[1])]= line[2]
                 Matrix[int(line[1]),int(line[0])]= line[2]
        for x in range(0,self.Dimension):
            for y in range(x+1,self.Dimension):
                reducedArray.append(Matrix[x,y])

        Distances = np.array(reducedArray, dtype=(float))
        self.Tree =hierarchy.linkage(Distances, 'complete')
        return self.Tree


    def avgTree(self):
        """
        Returns the HCA dendrogrm, using the average linkage method
        """        
        data = self.ccTable
        Matrix=np.zeros((self.Dimension,self.Dimension))

        reducedArray=[]
        for line in data:
                #print line
            if line is not None and len(line) != 0:
                 Matrix[int(line[0]),int(line[1])]= line[2]
                 Matrix[int(line[1]),int(line[0])]= line[2]
        for x in range(0,self.Dimension):
            for y in range(x+1,self.Dimension):
                reducedArray.append(Matrix[x,y])

        Distances = np.array(reducedArray, dtype=(float))

        self.Tree =hierarchy.linkage(Distances, 'average')

        return self.Tree


    def flatClusterPrinter(self, thr, labelsList, run_dir:str):
        """
        Prints the flat cluster at a chosen threshold to a .json file
        """        
        FlatC=hierarchy.fcluster(self.Tree, thr, criterion='distance')
        clusterToJson={}
        clusterToJson['HKL']=[]
        abs_run_dir = os.path.abspath(run_dir)
        print(f"run_dir: {abs_run_dir}")
        #check run dir
        if os.path.isdir(abs_run_dir):
            with open(f"{abs_run_dir}/flatCluster.json", 'w') as clusterFile:
                for cluster, hkl in zip(FlatC, labelsList):
                    clusterToJson['HKL'].append({
                        'input_file':hkl,
                        'cluster':str(cluster)
                        })
                print(f"prepare to convent cluster information to flatCluster.json: \n{clusterToJson}")
                json.dump(clusterToJson, clusterFile, indent=4)
        else:
            print(f"Run dir for flatClusterPrinter does not exist, please check: {abs_run_dir}")


    #function to pipe HCA into codgas for subsequent GA analysis.
    def passOInfoToGA(self, thr, labelsList, anomFlag):
            FlatC=hierarchy.fcluster(self.Tree, thr, criterion='distance')
            counter=collections.Counter(FlatC)
            Best = max(counter.items(), key=operator.itemgetter(1))[0]
            with open(self.RunDir+'/GA/codgas.INP', 'a') as codgasINP:
                for cluster, dataset in zip(FlatC, labelsList):
                    codgasINP.write(f"{dataset} {cluster}\n")


    def thrEstimation(self):
        """
        Estimates the threshold for optimal clustering, based on the multiplicity of the biggest cluster
        """             
        x = 0.00
        dx = 0.05
        countsList = []
        x_list = []
        while x < 1:
            
            FlatC = hierarchy.fcluster(self.Tree, x, criterion='distance')
            counter=collections.Counter(FlatC)
            Best = max(counter.items(), key=operator.itemgetter(1))[0]
            countsList.append(counter[Best])        
            x+= dx
            x_list.append(x)
        dy = np.diff(countsList)

        for a, b in zip (x_list, dy):
            if b == max(dy):
                return round(a, 2)


    def checkMultiplicity(self, thr):
        """
        Prints the multiplicity of the biggest cluster at a given threshold
        """                
        FlatC = hierarchy.fcluster(self.Tree, thr, criterion='distance')
        counter=collections.Counter(FlatC)
        Best = max(counter.items(), key=operator.itemgetter(1))[0]
        print('You are clustering with a threshold of %s'%(thr))
        print('The biggest cluster contains %s datasets from a total of %s'%(counter[Best], len(self.labelList)))
        return Best, counter[Best], len(self.labelList)


    def completenessEstimation(self):
        x = 0.00
        dx = 0.05
        while x > 1:
            FlatC = hierarchy.fcluster(self.Tree, x, criterion='distance')
            counter=collections.Counter(FlatC)
            Best = max(counter.items(), key=operator.itemgetter(1))[0]


    #the list self.ToProcess is needed by the scaling routines
    #fix all this new mess!
    #Tk: do we still need this function?
    def whatToProcess(self):
        FlatC = hierarchy.fcluster(self.Tree, thr, criterion='distance')        
        counter=collections.Counter(FlatC)
        Best = max(counter.items(), key=operator.itemgetter(1))[0]
        Process = True
        #change checkboxes to standard variables
        if Process:
            self.ToProcess = [Best]    
        else:
            self.ToProcess = set(Clusters)
            for key in self.ToProcess:
                if counter[key]==1:
                    self.ToProcess = [x for x in self.ToProcess if x != key]
        return self.ToProcess


    #save XSCALE.LP statistics plot in the gallery folder for data porte
    def SaveXscalePlot(self, ProcessDir:str, res, value, title:str, fileName:str):
        statsPlot = Figure()
        Ax = statsPlot.add_subplot(111)
        statsPlot.set_size_inches(16, 9)
        plotDataX= []
        plotDataY= []
        plotList, _ = extractXSCALEStat(f"{ProcessDir}/XSCALE.LP")

        #check if plotList is empty, if so, print a warning and return
        if not plotList:
            print(f"Warning: No data found in {ProcessDir}/XSCALE.LP")
            return 

        #Setup resolution and prepare the figure
        LowestRes = float(plotList[0][res]) + 0.25
        HighestRes = float(plotList[-1][res]) - 0.25 if float(plotList[-1][res]) > 0.5 else 0

        for line in plotList:
            plotDataX.append(float(line[res])) 
            plotDataY.append(float(line[value].strip('*').strip('%')))
        Ax.plot(plotDataX, plotDataY, 'r-^')

        #Set x and y labels, as well as the title
        Ax.set_xlabel("Resolution (Å)", fontsize=12, fontweight='bold', color='black')
        y_label = title.split(" vs ")[0]
        Ax.set_ylabel(y_label, fontsize=12, fontweight='bold', rotation=90, labelpad=20, color='black')
        Ax.set_title(title)

        Ax.set_xlim(LowestRes, HighestRes)
        Ax.set_title(title)

        #Check if gallery folder exists, if not create it, we need to save the dendrogram in the gallery folder for data porte
        PlotFile  = Path(f"{ProcessDir}/gallery/{fileName}.png")
        PlotFile.parent.mkdir(parents=True, exist_ok=True)

        #Save the plot
        statsPlot.savefig(PlotFile, bbox_inches="tight", dpi=300)
        plt.close(statsPlot)


    #Run XSCALE to merge the biggest cluster
    #input files
    #!!!! Will need to define the processes to run externally
    #renaming function! Edit the calls in ccCluster accordingly
    def prepareXSCALE(self, anomFlag, thr, clusterList:list=None, **kwargs):
        #check if thr is or can be converted to float, if not, print a warning and return
        try:
            threshold = round(float(thr), 2)
        except ValueError:
            print(f"{colors.RED}Error: Threshold value {thr} is not a valid float. Please check the input.{colors.ENDC}")
            return False, None

        FlatC = hierarchy.fcluster(self.Tree, threshold, criterion='distance') #takes threshold here to clustter the files
        print(f"FlatC is: {FlatC}")
        counter=collections.Counter(FlatC)
        print(f"counter is: {counter}")

        #set up selected cluster for merging
        if clusterList == None:
            self.SelectedCluster = [max(counter.items(), key=operator.itemgetter(1))[0]] #returns one item (group number) as list
        else:
            self.SelectedCluster = clusterList

        print(f"Checking selected cluster for XSCALE: {self.SelectedCluster}")
        #check to make sure at least one of the selected cluster is in the FlatC list, otherwise it will not be processed
        Checkitem = all(num not in FlatC for num in self.SelectedCluster)
        if Checkitem:
            print(f"{colors.RED}Warning: None of the selected clusters exist in the FlatC list, XSCALE will not start, please check{colors.ENDC}")
            return False, None

        for num in self.SelectedCluster:
            if num in FlatC:
                print(f"cluster number {num} exist, will be merged")
            else:
                print(f"cluster number {num} does not exist, will be ignored, please check")

        #Setup running directory
        clusterStr = "n".join(str(x) for x in self.SelectedCluster)
        processing_dir_XSCALE = f"{self.RunDir}/cc_Cluster_{threshold}_{clusterStr}_{anomFlag}"
        XSCALE_file = f"{processing_dir_XSCALE}/XSCALE.INP"
        
        #check working dir
        if os.path.isdir(processing_dir_XSCALE):
            print(f"Processing folder exists, checking content: {processing_dir_XSCALE}")
        else:
            os.mkdir(processing_dir_XSCALE)

        #check whether XSCALE.INP exists and skik the job if exists
        if os.path.isfile(XSCALE_file):
            print(f"XSCALE.INP exist: {XSCALE_file}\nWill pass the XSCALE process. Please reomve the file/folder if you want to re-run the job")
            return False, None
        else:
            #Create XSCALE.INP setups
            with open(XSCALE_file, 'a') as Xscale:
                Xscale.write(f"OUTPUT_FILE= xscale_scaled.hkl\n")
                Xscale.write(f"MERGE= TRUE\n")
                if anomFlag=='ano':
                    Xscale.write(f"FRIEDEL\'S_LAW= FALSE\n")
                elif anomFlag=='no_ano':
                    Xscale.write(f"FRIEDEL\'S_LAW= TRUE\n")
                if kwargs.get("refHKL"):
                    reference_HKL = os.path.abspath(kwargs["refHKL"])
                    print(f"Reference HKL file exists, added to XSCALE.INP: {reference_HKL}")
                    Xscale.write(f"REFERENCE_DATA_SET= {reference_HKL}\n")
                else:
                    print(f"No optional reference HKL for XSCALE, continue")

                #put HKL files for merging in the XSCALE.INP It should be OK to line in the self.Toprocess loop
                # as there is only one item in [Best]
                for cluster, filename in zip(FlatC, self.labelList):
                    if cluster in self.SelectedCluster:
                        Xscale.write(f"INPUT_FILE= {filename[1]}\n")
                        if kwargs.get("resolutionRange") and len(kwargs["resolutionRange"]) == 2:
                            Xscale.write(f"INCLUDE_RESOLUTION_RANGE= {kwargs['resolutionRange'][0]} {kwargs['resolutionRange'][1]}\n")
                        elif kwargs.get("resolutionRange") and len(kwargs["resolutionRange"]) != 2:
                            print(f"{colors.RED}Warning: Invalid resolution range provided: {kwargs['resolutionRange']}. It should be a list of two numbers.{colors.ENDC}")
                        else:
                            print(f"Resolution range not provided, skipping INCLUDE_RESOLUTION_RANGE in XSCALE.INP")
                        #Xscale.write(f"INCLUDE_RESOLUTION_RANGE= 20, 1.8\n")
                        #Xscale.write(f"MINIMUM_I/SIGMA= 0\n")

            return True, processing_dir_XSCALE

    #read the f2mtz and cad commands from xdsconv generatet F2MTZ.INP
    def readXDSCONVLP(self, XDSCONVLP:str, outputName:str):
        f2mtz_cmd = ""
        cad_cmd = ""
        cad_input = ""
        if not os.path.isfile(XDSCONVLP):
            print(f"{colors.RED}Error: XDSCONV log file {XDSCONVLP} does not exist.{colors.ENDC}")
            return None, None
        else:
            with open(XDSCONVLP, "r") as f:
                for line in f:
                    if line.strip().startswith("f2mtz"):
                        print(f"Found f2mtz command: {line.strip()}")
                        f2mtz_cmd = line.replace("<", " ").strip().split()
                        f2mtz_cmd.pop()  # remove the last element which is the input file
                        break
                for line in f:
                    if line.strip().startswith("cad"):
                        print(f"Found cad command: {line.strip()}")
                        cad_cmd = line.strip().replace("<", " ").split()
                        cad_cmd.pop()  # remove the last element which is the input file
                        cad_cmd[-1] = outputName  # add the output file name as the last argument
                        break
                for line in f:
                    print(f"Found cad command: {line.strip()}")
                    cad_input += f"{line.strip()}\n"
        return f2mtz_cmd, cad_cmd, cad_input

    
    #convert xscale output to mtz using xdsconv
    def xscaleToMtz(self, XscaleHKLPath:str, run_dir:str, anomFlag:str, **kwargs):
        abs_run_dir = os.path.abspath(run_dir)
        if not os.path.isdir(abs_run_dir):
            print(f"{colors.RED}Error: Run directory {abs_run_dir} does not exist.{colors.ENDC}")
            return

        xdsconv_file = os.path.join(abs_run_dir, "XDSCONV.INP")
        if os.path.isfile(xdsconv_file):
            print(f"XDSCONV.INP already exists: {xdsconv_file}\nWill pass the XDSCONV process. Please reomve the file/folder if you want to re-run the job")
            return
        else:
            if not os.path.isfile(XscaleHKLPath):
                print(f"{colors.RED}Error: Xscale HKL file {XscaleHKLPath} does not exist.{colors.ENDC}")
                return
            
            with open(xdsconv_file, 'w') as Xconv:
                Xconv.write(f"INPUT_FILE= {XscaleHKLPath}\n")
                Xconv.write(f"OUTPUT_FILE= tmp.hkl CCP4\n")
                if anomFlag=='ano':
                    Xconv.write(f"FRIEDEL\'S_LAW= FALSE\n")
                elif anomFlag=='no_ano':
                    Xconv.write(f"FRIEDEL\'S_LAW= TRUE\n")

            print(f"Running XDSCONV with input {XscaleHKLPath}")
            subprocess.run('xdsconv',cwd=abs_run_dir)

            #make sure the necessary files exist before proceeding
            if not os.path.isfile(os.path.join(abs_run_dir, "XDSCONV.LP")):
                print(f"{colors.RED}Error: XDSCONV.LP not found after running xdsconv in {abs_run_dir}.{colors.ENDC}")
                return
            if not os.path.isfile(os.path.join(abs_run_dir, "F2MTZ.INP")):
                print(f"{colors.RED}Error: F2MTZ.INP not found after running xdsconv in {abs_run_dir}.{colors.ENDC}")
                return
            
            f2mtz_cmd, cad_cmd, cad_input = self.readXDSCONVLP(os.path.join(abs_run_dir, "XDSCONV.LP"), "xscale_scaled.mtz")
            print(f"Extracted commands from XDSCONV.LP:\nF2MTZ: {f2mtz_cmd}\nCAD: {cad_cmd}\nCAD Input: {cad_input}")

            #run the f2mtz command
            with open(os.path.join(abs_run_dir, "F2MTZ.INP"), 'r') as f:
                f2mtz_input = f.read()
            subprocess.run(f2mtz_cmd, input=f2mtz_input, text=True, cwd=abs_run_dir)

            #run the cad command
            subprocess.run(cad_cmd, input=cad_input, text=True, cwd=abs_run_dir)
        

    def preparePointless(self, anomFlag, thr, clusterList:list=None, **kwargs):
        FlatC = hierarchy.fcluster(self.Tree, thr, criterion='distance')
        counter=collections.Counter(FlatC)

        #check if thr is or can be converted to float, if not, print a warning and return
        try:
            threshold = round(float(thr), 2)
        except ValueError:
            print(f"{colors.RED}Error: Threshold value {thr} is not a valid float. Please check the input.{colors.ENDC}")
            return False, None

        #set up selected cluster for merging
        if clusterList == None:
            self.SelectedCluster = [max(counter.items(), key=operator.itemgetter(1))[0]] #returns one item (group number) as list
        else:
            self.SelectedCluster = clusterList

        print(f"Checking selected cluster for XSCALE: {self.SelectedCluster}")
        #check to make sure at least one of the selected cluster is in the FlatC list, otherwise it will not be processed
        Checkitem = all(num not in FlatC for num in self.SelectedCluster)
        if Checkitem:
            print(f"{colors.RED}Warning: None of the selected clusters exist in the FlatC list, Aimless and Pointless will not run, please check{colors.ENDC}")
            return False, None
        
        for num in self.SelectedCluster:
            if num in FlatC:
                print(f"cluster number {num} exist, will be merged")
            else:
                print(f"cluster number {num} does not exist, will be ignored, please check")

        #Check whether folder/file exists, and prepare it if not
        #Setup running folder
        clusterStr = "n".join(str(x) for x in self.SelectedCluster)
        processing_dir_Pointless = f"{self.RunDir}/cc_Cluster_{threshold}_{clusterStr}_{anomFlag}"
        Pointless_file = f"{processing_dir_Pointless}/launch_pointless.sh"
        if os.path.isdir(processing_dir_Pointless):
            print(f"Processing folder exists, checking content: {processing_dir_Pointless}")
        else:
            os.mkdir(processing_dir_Pointless)

        #Check whether launch_pointless exists
        if os.path.isfile(Pointless_file):
            print(f"launch_pointless.sh exist: {Pointless_file}\nWill pass the Pointless process. Please reomve the file/folder if you want to re-run the job")
            return False, None
        else:
            with open(Pointless_file, 'a') as Pointless:
                Pointless.write(f"pointless hklout pointless_clustered.mtz << EOF\n")
                Pointless.write(f"XMLOUT pointlessLog.xml\n")

                #put HKL files for merging in the XSCALE.INP It should be OK to line in the self.Toprocess loop
                for cluster, filename in zip(FlatC,self.labelList):
                    #if cluster in self.ToProcess:
                    if cluster in self.SelectedCluster:
                        Pointless.write(f"HKLIN {filename[1]}\n")
                        #Pointless.write(f"EOF\n")

            return True, processing_dir_Pointless


    #Run XSCALE in the pre-determined folders, not self.Rundir ().
    def scaleAndMerge(self, anomFlag, thr, run_dir:str):
        print(f"Selected cluster number: {self.SelectedCluster}")
        abs_run_dir = os.path.abspath(run_dir)
        xscale_file = f"{abs_run_dir}/XSCALE.INP"

        #threshold should be a float, check it just in case
        try:
            threshold = round(float(thr), 2)
        except ValueError:
            print(f"{colors.RED}Error: Threshold value {thr} is not a valid float. Please check the input.{colors.ENDC}")
            return

        if os.path.isdir(abs_run_dir):
            if os.path.isfile(xscale_file):
                #self.createDendrogram(thr)
                X = hierarchy.dendrogram(self.Tree, color_threshold=threshold)
                plt.gcf().set_size_inches(16, 9)

                #try to match the color to cluster number using count
                #As they contains same amount od datasets
                #get the FlatC with the same threshold
                FlatC = hierarchy.fcluster(self.Tree, threshold, criterion='distance')
                cluster_counts = collections.Counter(FlatC)
                color_counts = collections.Counter(X['leaves_color_list'])

                #remove unwanted color
                for remove_color in ['C0', 'k', 'grey']:
                    if remove_color in color_counts:
                        del color_counts[remove_color]

                #remove cluster with only 1 dataset
                filtered_clusters = {cluster: count for cluster, count in cluster_counts.items() if count > 1}

                #Sort by count
                sorted_clusters = sorted(filtered_clusters.items(), key=lambda x: x[1], reverse=True)
                sorted_colors = sorted(color_counts.items(), key=lambda x: x[1], reverse=True)
                print(f"sorted cluster: {sorted_clusters}")
                print(f"sorted color: {sorted_colors}")

                #Match clusters to colors by count
                cluster_to_color = {}
                for (cluster, cluster_count), (color, color_count) in zip(sorted_clusters, sorted_colors):
                    if cluster_count == color_count:
                        cluster_to_color[cluster] = color
                    else:
                        print(f"Warning: Count mismatch! Cluster {cluster} has {cluster_count}, Color {color} has {color_count}")
                        cluster_to_color[cluster] = color

                #Create legend using matched clusters
                legend_handles = []
                
                #add x and y labels, as well as the threshold line
                plt.xlabel("Dateset numbers", fontsize=12, fontweight='bold', color='black')
                plt.ylabel("Correlation coefficients", fontsize=12, fontweight='bold', rotation=90, labelpad=20, color='black')
                plt.axhline(y=threshold, color='grey', linestyle='--', alpha=0.5, linewidth=2)
                threshold_handle = Line2D([0], [0], color='grey', linestyle='--', alpha=1, linewidth=2, label=f'$\\mathbf{{Thresh:\\ {threshold}}}$')
                legend_handles.append(threshold_handle)
    
                #match cluster to color and add to legend
                for cluster in sorted(cluster_to_color.keys()):
                    color = cluster_to_color[cluster]
                    legend_handles.append(mpatches.Patch(color=color, label=f"$\\mathbf{{Cluster\\ {cluster}}}$ :\n{cluster_counts[cluster]} datasets"))

                #Add legend
                legend = plt.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(1.02, 1),
                            borderaxespad=0, fontsize="small", handleheight=3, handlelength=3,
                            title_fontsize="medium")

                #Check if gallery folder exists, if not create it, we need to save the dendrogram in the gallery folder for data porte
                Dendrogramfile  = Path(f"{abs_run_dir}/gallery/Dendrogram.png")
                Dendrogramfile.parent.mkdir(parents=True, exist_ok=True)

                #save dendrogram with legend
                plt.savefig(Dendrogramfile, bbox_inches="tight", dpi=300, bbox_extra_artists=[legend], pad_inches=0.3)
                plt.close()

                #run XSCALE in the pre-determined folders, not self.Rundir ().
                print(f"Running XSCALE in {abs_run_dir}, please be patient\n")
                subprocess.run('xscale_par',cwd=abs_run_dir)

                #save the XSCALE.LP statistics plot in the gallery folder for data porte
                #plot from XSCLE.LP, position 0 is resolution, 4 is COMPLETENESS
                #5 is R factor observed, 8 is <I/sigma>, 10 is CC1/2, 12 is SigAno
                self.SaveXscalePlot(abs_run_dir, 0, 10, "CC_vs_Res", "CC_vs_Res")
                self.SaveXscalePlot(abs_run_dir, 0, 4, "comp_vs_Res", "comp_vs_Res")
                self.SaveXscalePlot(abs_run_dir, 0, 5, "Robs_vs_Res", "Robs_vs_Res")
                self.SaveXscalePlot(abs_run_dir, 0, 8, "<I/\u03C3I>_vs_Res", "I_SigmaI_vs_Res")
                self.SaveXscalePlot(abs_run_dir, 0, 12, "Sig_Ano_vs_Res", "Sig_Ano_vs_Res")

                #run XDSCONV to convert the XSCALE output to mtz format
                XscaleHKLPath = os.path.join(abs_run_dir, "xscale_scaled.hkl")

                #check if the XSCALE output file exists before running XDSCONV
                if not os.path.isfile(XscaleHKLPath):
                    print(f"{colors.RED}XSCALE output file is missing: {XscaleHKLPath}, please check the XSCALE run{colors.ENDC}")
                    return
                #convert the XSCALE output to mtz format using XDSCONV
                self.xscaleToMtz(XscaleHKLPath, abs_run_dir, anomFlag)
            else:
                print(f"ccCluster XSCALE.INP file does not exist, please check: {abs_run_dir}")
        else:
            print(f"ccCluster XSCALE run dir is missing, please check: {abs_run_dir}")


    #run Pointless in each folder from the processing List
    def pointlessRun(self, anomFlag, thr, run_dir:str):
        print(f"Selected cluster number: {self.SelectedCluster}")
        abs_run_dir=os.path.abspath(run_dir)
        pointless_file = f"{abs_run_dir}/launch_pointless.sh"
        Pointless_log = f"{abs_run_dir}/pointless.log"
        if os.path.isdir(abs_run_dir):
            if os.path.isfile(pointless_file):
                with open(pointless_file, 'a') as Pointless:
                    Pointless.write(f"COPY\nbg\nTOLERANCE 4\nEOF\n")
                print(f"Running Pointless, please be patient\n")
                with open(Pointless_log, "a") as Pointlesslog:
                    Pointless_command = ["bash", f"{pointless_file}"]
                    subprocess.run(Pointless_command, stdout=Pointlesslog, check=True, cwd=abs_run_dir)
            else:
                print(f"Looks like no lunch_pointless.sh in: {abs_run_dir}\nPlease check")
        else:
            print(f"ccCluster Pointless run dir is missing, please check: {abs_run_dir}")


    #extract unit cell from poineless.log for aimless, an ugly way but works
    def ExtracUnitCell(self, aimless_log:str):
        lines = open(aimless_log).readlines()
        lines_strip = [a.strip() for a in lines] #remove the space at the beginning of the line for searching keywords
        lines_tidy = [' '.join(re.sub(r"(\[)|(\])", " ", x.replace(" - ", " ")).split()) for x in lines_strip]

        counter = 0

        for line in lines_tidy:
            if "Space group =" in line:
                sg_line= line.split("\'")
                sg = sg_line[1]
            elif "Cell Dimensions " in line:
                counter += 3

            if counter > 1:
                counter -= 1

            elif counter == 1:
                counter -= 1
                uc_line = line.split("\'")
                uc = uc_line[0]
        return sg, uc


    #run aimless on the output from pointless
    #will run in folders with clustered.mtz file available.
    #TBD: fix directories paths into the aimless.inp file
    #also set all the proper input values into the function call
    #path to aimless executable to be verified.
    def aimlessRun(self, anomFlag, thr, run_dir:str, **kwargs):
        #Get variables, maybe can just put them in the file?
        infile = "pointless_clustered.mtz"
        setname = "aimless_clustered"
        if kwargs.get("resolutionRange") and len(kwargs["resolutionRange"]) == 2:
            resLow = str(kwargs["resolutionRange"][0])
            resHigh = str(kwargs["resolutionRange"][1])
        elif kwargs.get("resolutionRange") and len(kwargs["resolutionRange"]) != 2:
            print(f"{colors.RED}Warning: Invalid resolution range provided: {kwargs['resolutionRange']}. It should be a list of two numbers.\n"
                  f"Using default resolution range for aimless: 1.0 - 60.0{colors.ENDC}")
            resHigh = "1.0"
            resLow = "60"
        else:
            print(f"Resolution range not provided, use default for aimless: 1.0 - 60.0")
            resHigh = "1.0"
            resLow = "60"
        if anomFlag=='ano':
            anomflag = "ON"
        elif anomFlag=='no_ano':
            anomflag = "OFF"
        abs_run_dir = os.path.abspath(run_dir)
        aimless_file = f"{abs_run_dir}/aimless.inp"
        aimless_log = f"{abs_run_dir}/aimless.log"
        #get uc param
        pointless_log = f"{abs_run_dir}/pointless.log"
        sg, cell = self.ExtracUnitCell(pointless_log)
        SpaceGroup = f"\"{sg}\"" #need to have the ""
        if os.path.isdir(abs_run_dir):
            if os.path.isfile(aimless_file):
                print(f"aimless.inp already exists in rundir, stop re-run the job: {abs_run_dir}")
            else:
                with open(aimless_file, "a") as f1:
                    f1.write(textwrap.dedent(f"""\
                                            #!/bin/bash

                                            aimless HKLIN {infile} << EOF
                                            HKLOUT {setname}_aimless.mtz
                                            RESOLUTION LOW {resLow} HIGH {resHigh}
                                            OUTPUT MERGED
                                            anomalous {anomflag}
                                            EOF

                                            #truncate: generate Fs
                                            truncate hklin {setname}_aimless.mtz hklout {setname}_tr.mtz <<EOF-trunc
                                            truncate yes
                                            EOF-trunc


                                            #unique: generate unique reflection set for rfree
                                            unique HKLOUT {setname}_unq.mtz << EOF
                                            CELL {cell}
                                            SYMMETRY {SpaceGroup}
                                            LABOUT F=FUNI SIGF=SIGFUNI
                                            RESOLUTION {resHigh}
                                            EOF

                                            #freerflag: generate free reflections
                                            freerflag HKLIN {setname}_unq.mtz HKLOUT {setname}_FreeR_unq.mtz <<EOF
                                            FREERFRAC 0.05
                                            END
                                            EOF

                                            #cad: combine free reflections with data
                                            cad HKLIN1 {setname}_FreeR_unq.mtz HKLIN2 {setname}_tr.mtz HKLOUT {setname}_cad.mtz<<EOF
                                            LABI FILE 1 E1=FreeR_flag
                                            LABI FILE 2 ALLIN
                                            END
                                            EOF

                                            freerflag HKLIN {setname}_cad.mtz HKLOUT {setname}_scaled.mtz <<EOF
                                            COMPLETE FREE=FreeR_flag
                                            END
                                            EOF
                                        """))
                print(f"Running aimless, please be patient\n")
                with open(aimless_log, "a") as aimlesslog:
                    aimless_command = ["bash", f"{aimless_file}"]
                    subprocess.run(aimless_command, stdout=aimlesslog, check=True, cwd=abs_run_dir)
        else:
            print(f"Aimless run dir (pointless outputdir) does not exists, please check: {abs_run_dir}")


    #A function to investigate the influence of reference file in merging results
    def shuffleXscale(self, anomFlag, thr):
        FlatC = hierarchy.fcluster(self.Tree, thr, criterion='distance')
        run_dir_shuffle = f"{self.RunDir}/xscale_shuffle"
        Log = open(run_dir_shuffle+'/.cc_cluster.log', 'a') #not sure anything is written into it leave it like this for now
        counter=collections.Counter(FlatC)
        self.Best = [max(counter.items(), key=operator.itemgetter(1))[0]]
        print(self.Best)
        Process = True
        xscaleInputFiles=[]

        #Prepare list of filenames to shuffle over
        for cluster, filename in zip(FlatC, self.labelList):
            if cluster == self.Best[0]:
                xscaleInputFiles.append(filename[1])
        print(xscaleInputFiles)

        #run XSCALE with random ordered files 20 times
        #check whether run_dir already exists
        for x in range(0,20):
            if os.path.isdir(run_dir_shuffle):
                print(f"XSCALE shuffel folder exists, checking content: {run_dir_shuffle}")
            else:
                os.makedirs(run_dir_shuffle)
            #check whether run dir is empty
            if os.listdir(run_dir_shuffle):
                print(f"XSCALE shuffle folder exists and is not empty, shuffle stopped, Please check: {run_dir_shuffle}")
            else:
                os.makedirs(run_dir_shuffle+'/thr_%.2f_run_%s'%(float(thr),x))
                with open(run_dir_shuffle+'/thr_%.2f_run_%s/XSCALE.INP'%(float(thr),x), 'a') as Xscale:
                    Xscale.write(f"OUTPUT_FILE=scaled.hkl\n")
                    Xscale.write(f"MERGE= TRUE\n")
                    if anomFlag == "ano":
                        Xscale.write(f"FRIEDEL\'S_LAW=FALSE\n")
                    elif anomFlag == "no_ano":
                         Xscale.write(f"FRIEDEL\'S_LAW=TRUE\n")

                    random.shuffle(xscaleInputFiles)
                    for hkl in xscaleInputFiles:
                        Xscale.write(f"INPUT_FILE= {hkl}\n")
                        #Xscale.write(f"MINIMUM_I/SIGMA= 0\n")
                subprocess.run('xscale_par',cwd=self.CurrentDir+'/thr_%.2f_run_%s'%(float(thr),x))

#seems we do not need to run the main function, as we will call the clustering function from the ccCluster.py
"""def main():
    from optparse import OptionParser
    parser = OptionParser(usage="%prog --XSCALEfile=<LP filename> --outname=<output dendogram>")

    parser.add_option("-o","--outname", dest="outname", default='Dendrogram', help="output dendogram file name")
    parser.add_option("-t", "--threshold", dest="threshold", default='0.4', help="Distance threshold for clustering")
    parser.add_option("-c", "--count",action="store_true", dest="count", default=False, help="Counts datasets in the biggest cluster and exit")
    (options, args) = parser.parse_args()

    thr = float(options.threshold)
    CC = Clustering('Cluster_log.txt')
    link = CC.tree()
    if options.count:
        CC.checkMultiplicity(thr)
        print(CC.thrEstimation())
    else:
        CC.checkMultiplicity(thr) 
        CC.merge('ano', thr)

if __name__== '__main__':
    main()"""
