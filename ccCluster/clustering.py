from __future__ import print_function
__author__ = "Gianluca Santoni"
__copyright__ = "Copyright 20150-2019"
__credits__ = ["Gianluca Santoni, Alexander Popov"]
__license__ = ""
__version__ = "1.0"
__maintainer__ = "Gianluca Santoni"
__email__ = "gianluca.santoni@esrf.fr"
__status__ = "Beta"




from scipy.cluster import hierarchy
import scipy
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
import numpy as np
import subprocess
import collections
import operator
import stat
import json
import random
import textwrap
import re

# from .report  import WorkflowStepReport

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
        self.previousProcess()


    def previousProcess(self):
        """
        Lists all the clusters which have already been processed from a log file.
        Updates the global variable alreadyDone
        """
        self.alreadyDone= []
        if os.path.isfile(os.getcwd()+'/.cc_cluster.log'):
            with open(os.getcwd()+'/.cc_cluster.log') as log:
                for line in log:
                    L = line.split(',')
                    self.alreadyDone.append([L[1], L[2].strip(), L[3].strip()])


    def parseCCFile(self):
        """
        Gets data from ccCalc ouput file and populates a numpy array with the distances
        """
        with open(self.ccFile, 'r') as f:
            dataArr = None
            data=[]
            Index = []
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
                self.labelList.append("%s"%(goodLine[2].strip('\n')))
        return self.labelList
        

    def inputType(self):
        """
        return input file type. Either mtz or HLK
        """        
        element = self.labelList[0]
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


    def flatClusterPrinter(self, thr, labelsList, anomFlag, run_dir:str):
        """
        Prints the flat cluster at a chosen threshold to a .json file
        """        
        FlatC=hierarchy.fcluster(self.Tree, thr, criterion='distance')
        counter=collections.Counter(FlatC)
        clusterToJson={}
        clusterToJson['HKL']=[]
        abs_run_dir = os.path.abspath(run_dir)
        print(f"run_dir: {abs_run_dir}")
        #check run dir
        Best = max(counter.items(), key=operator.itemgetter(1))[0]
        if os.path.isdir(abs_run_dir):
            with open(f"{abs_run_dir}/flatCluster.json", 'w') as clusterFile:
                for cluster, hkl in zip(FlatC, labelsList):
                    clusterToJson['HKL'].append({
                        'input_file':hkl,
                        'cluster':str(cluster)
                        })
                print(f"prepare to convent cluster information toflatCluster.json: \n{clusterToJson}")
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
                return a


    def checkMultiplicity(self, thr):
        """
        Prints the multiplicity of the biggest cluster at a given threshold
        """                
        FlatC = hierarchy.fcluster(self.Tree, thr, criterion='distance')
        counter=collections.Counter(FlatC)
        Best = max(counter.items(), key=operator.itemgetter(1))[0]
        print('You are clustering with a threshold of %s'%(thr))
        print('The biggest cluster contains %s datasets from a total of %s'%(counter[Best], len(self.labelList)))


    def completenessEstimation(self):
        x = 0.00
        dx = 0.05
        while x > 1:
            FlatC = hierarchy.fcluster(self.Tree, x, criterion='distance')
            counter=collections.Counter(FlatC)
            Best = max(counter.items(), key=operator.itemgetter(1))[0]


    #the list self.ToProcess is needed by the scaling routines
    #fix all this new mess!
    #Tk: do we still need this function? may be use self.Best instead, as it will be a single element anyway
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


    #Run XSCALE to merge the biggest cluster
    #input files
    #!!!! Will need to define the processes to run externally
    #renaming function! Edit the calls in ccCluster accordingly
    def prepareXSCALE(self, anomFlag, thr):
        FlatC = hierarchy.fcluster(self.Tree, thr, criterion='distance') #takes threshold here to clustter the files
        print(f"FlactC is: {FlatC}")
        counter=collections.Counter(FlatC)
        print(f"counter is: {counter}")
        self.Best = max(counter.items(), key=operator.itemgetter(1))[0] #returns onle one item (group number)
        #Process = True # Maybe not needed?
        #change checkboxes to standard variables
        #Do we need this as Process will be True anyway it sets up self.Toprocess but it is same as best
        #Keep this one for now but remove te rest, in case we need to revive it
        """
        if Process:
            self.ToProcess = [Best]    
        else:
            self.ToProcess = set(Clusters)
            for key in self.ToProcess:
                if counter[key]==1:
                    self.ToProcess = [x for x in self.ToProcess if x != key]
        #for x in self.ToProcess:
        """

        #Setup running directory
        processing_dir_XSCALE = self.RunDir+'/cc_Cluster_%.2f_%s_%s'%(float(thr),self.Best, anomFlag)
        XSCALE_file = f"{processing_dir_XSCALE}/XSCALE.INP"
        if [thr, self.Best, anomFlag] not in self.alreadyDone:
            #check working dir
            if os.path.isdir(processing_dir_XSCALE):
                print(f"Processing folder exists, checking content: {processing_dir_XSCALE}")
            else:
                os.mkdir(processing_dir_XSCALE)

            #check whether XSCALE.INP exists and skik the job if exists
            if os.path.isfile(XSCALE_file):
                print(f"XSCALE.INP exist: {XSCALE_file}\nWill pass the XSCALE process. Please reomve the file/folder is you want to re-run the job")
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

                #put HKL files for merging in the XSCALE.INP It should be OK to line in the self.Toprocess loop
                # as there is only one item in [Best]
                    for cluster, filename in zip(FlatC, self.labelList):
                        if cluster == self.Best:
                            Xscale.write(f"INPUT_FILE= {filename}\n")
                            #Xscale.write(f"INCLUDE_RESOLUTION_RANGE= 20, 1.8\n")
                            #Xscale.write(f"MINIMUM_I/SIGMA= 0\n")

                return True, processing_dir_XSCALE


    def preparePointless(self, anomFlag, thr):
        FlatC = hierarchy.fcluster(self.Tree, thr, criterion='distance')
        counter=collections.Counter(FlatC)
        self.Best = max(counter.items(), key=operator.itemgetter(1))[0]

        #Check whether folder/file exists, and prepare it if not
        #Setup running folder
        processing_dir_Pointless = self.RunDir+'/cc_Cluster_%.2f_%s_%s'%(float(thr), self.Best, anomFlag)
        Pointless_file = f"{processing_dir_Pointless}/launch_pointless.sh"
        if [thr, self.Best, anomFlag] not in self.alreadyDone: #need to check what is this
            if os.path.isdir(processing_dir_Pointless):
                print(f"Processing folder exists, checking content: {processing_dir_Pointless}")
            else:
                os.mkdir(processing_dir_Pointless)

            #Check whether launch_pointless exists
            if os.path.isfile(Pointless_file):
                print(f"launch_pointless.sh exist: {Pointless_file}\nWill pass the XSCALE process. Please reomve the file/folder is you want to re-run the job")
                return False, None
            else:
                with open(Pointless_file, 'a') as Pointless:
                    Pointless.write(f"pointless hklout pointless_clustered.mtz << EOF\n")
                    Pointless.write(f"XMLOUT pointlessLog.xml\n")

                    #put HKL files for merging in the XSCALE.INP It should be OK to line in the self.Toprocess loop
                    # as there is only one item in [Best]
                    for cluster, filename in zip(FlatC,self.labelList):
                        #if cluster in self.ToProcess:
                        if cluster == self.Best:
                            Pointless.write(f"HKLIN {filename}\n")
                            #Pointless.write(f"EOF\n")

                return True, processing_dir_Pointless


    #Run XSCALE in the pre-determined folders, not self.Rundir ().
    def scaleAndMerge(self, anomFlag, thr, run_dir:str):
        print(f"Best cluster number: {self.Best}")
        newProcesses=[] #what does it do? - maybe used in other functions?
        abs_run_dir = os.path.abspath(run_dir)
        xscale_file = f"{abs_run_dir}/XSCALE.INP"
        if os.path.isdir(abs_run_dir):
            if os.path.isfile(xscale_file):
                if [thr, self.Best, anomFlag] not in self.alreadyDone:
                    #self.createDendrogram(thr)
                    X = hierarchy.dendrogram(self.Tree, color_threshold=float(thr))
                    plt.savefig(abs_run_dir+'/Dendrogram.png')
                    subprocess.run('xscale_par',cwd=abs_run_dir)
                    newProcesses.append([thr, self.Best, anomFlag])
            else:
                print(f"ccCluster XSCALE.INP file does not exist, please check: {abs_run_dir}")
        else:
            print(f"ccCluster XSCALE run dir is missing, please check: {abs_run_dir}")


    #run Pointless in each folder from the processing List
    def pointlessRun(self, anomFlag, thr, run_dir:str):
        print(f"Best cluster number: {self.Best}")
        newProcesses=[]
        abs_run_dir=os.path.abspath(run_dir)
        pointless_file = f"{abs_run_dir}/launch_pointless.sh"
        Pointless_log = f"{abs_run_dir}/pointless.log"
        if os.path.isdir(abs_run_dir):
            if [thr, self.Best, anomFlag] not in self.alreadyDone:
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
    def aimlessRun(self, anomFlag, thr, run_dir:str):
        #Get variables, maybe can just put them in the file?
        infile = "pointless_clustered.mtz"
        setname = "aimless_clustered"
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
                if [thr, self.Best, anomFlag] not in  self.alreadyDone:
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
        Log = open(run_dir+'/.cc_cluster.log', 'a') #not sure anything is written into it leave it like this for now
        counter=collections.Counter(FlatC)
        self.Best = max(counter.items(), key=operator.itemgetter(1))[0]
        print(self.Best)
        Process = True
        xscaleInputFiles=[]

        #change checkboxes to standard variables
        """
        if Process:
            self.ToProcess = [Best]    
        else:
            self.ToProcess = set(Clusters)
            for key in self.ToProcess:
                if counter[key]==1:
                    self.ToProcess = [x for x in self.ToProcess if x != key]
        """

        #Prepare list of filenames to shuffle over
        for cluster, filename in zip(FlatC, self.labelList):
            if cluster == self.Best:
                xscaleInputFiles.append(filename)
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


def main():
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
    main()
