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




import matplotlib.pyplot as plt
import sys
import os
import time
import subprocess
from os.path import join , isfile
#import ccCluster classes to run the soft
from .ccCalc import ccList
from .clustering import Clustering
import argparse


#pass arguements
def process_args():
    input_args = argparse.ArgumentParser()

    input_args.add_argument("-i","--input_folder",
    nargs="+",
    help="give the input folder that contains the processed files as list, make sure you have output.composite_output=False for data processing"
    )
    

    input_args.add_argument("-g", "--grenade", 
    action="store_true",  
    help="Whether use grenade folder, default is False"
    )


    input_args.add_argument("-f", "--file_list", 
    type=str, 
    nargs="+", 
    help="The list of refined structures (files) to merge"
    )

    input_args.add_argument("-o", "--output_dir",
    type = str,
    help = "output directory, default is pwd.",
    )

    input_args.add_argument("-t", "--threshold",
    type = str,
    help = "Distance threshold for clustering",
    )

    input_args.add_argument("-ano", "--anomalous_scattering",
    action="store_true",
    help = "Whether it is anomalous data, default is False",
    )

    #save args
    args = input_args.parse_args()

    #set working folder
    if args.output_dir == None:
        args.output_dir = os.getcwd()

    #return args    
    return args


def main():
    args = process_args()

    #set up abs file list for non-grenade files
    abs_file_list = [os.path.abspath(f) for f in args.file_list if os.path.isfile(f)]
    print(f"abs_path_list: {abs_file_list}")

    #set up anomalous
    if args.anomalous_scattering == False:
        anomlous = "no_ano"
    else:
        anomlous = "ano"

    #create working dir
    workdir = os.path.abspath(f"{args.output_dir}/HCA")
    print(f"ccCluster output directory would be : {workdir}")
    if os.path.isdir(workdir):
        print(f"Output folder exists: {workdir}")
        print(f"moving into output folder")
        os.chdir(workdir)
    else:
        print(f"Creating output folder: {workdir}")
        os.makedirs(workdir)
        print(f"moving into output folder")
        os.chdir(workdir)

    #set up the job
    if args.grenade == True:        
        grenades_runs = [join(workdir,x) for x in os.listdir(workdir) if 'grenades' in x]
        print(f"grenade_run: {grenades_runs}")
        success = []
        failed = []
        shouldContinue = True
        while shouldContinue==True:
            time.sleep(5)
            for path in grenades_runs :
                if isfile(path+'/.SUCCESS') :
                    success.append(path)
                elif isfile(path+'/.FAILED'):
                    failed.append(path)
            if len(success)+len(failed)==len(grenades_runs):
                shouldContinue = False
        #print(f"successlist: {success}")
        abs_success = [os.path.abspath(f) for f in success if os.path.isfile(f) ]
        print(f"abs_path_list: {abs_success}")
        ccList(success)

    elif args.grenade == False:
        ccList(abs_file_list)
    
    correlationFile='ccClusterLog.txt'
    CC = Clustering(correlationFile)
    Tree = CC.avgTree()
    etiquets=CC.createLabels()
    threshold = CC.thrEstimation()
    fileType = CC.inputType()
    if fileType=="HKL":
        CC.prepareXSCALE(anomlous,threshold)
        CC.scaleAndMerge(anomlous,threshold)
        CC.flatClusterPrinter(threshold, etiquets, anomlous)
    elif fileType=="mtz":
        CC.preparePointless(anomlous,threshold)
        CC.pointlessRun(anomlous,threshold)
        CC.flatClusterPrinter(threshold, etiquets, anomlous)
    


if __name__ =='__main__' :
    main()

    
