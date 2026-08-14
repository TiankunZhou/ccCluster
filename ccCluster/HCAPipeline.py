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
    
    input_args.add_argument("-ref", "--reference_HKL",
    type = str,
    help="Give an optional reference HKL file for XSCALE merging, recommend to give absolute path"
    )

    input_args.add_argument("-clu", "--clusters",
    nargs='*',
    type=int,
    help="pass Selected cluster number as a list: e.g. -clu 1 5 3 2"
    )

    #save args
    args = input_args.parse_args()

    #set working folder
    if args.output_dir == None:
        args.output_dir = os.getcwd()

    #setup selected cluster:
    if args.clusters is None or len(args.clusters) == 0:
        args.clusters == None

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

    #create working dir, make it abspath
    workdir = os.path.abspath(f"{args.output_dir}/HCA")
    print(f"ccCluster output directory would be : {workdir}")
    if os.path.isdir(workdir):
        print(f"Output folder exists: {workdir}")
        #print(f"moving into output folder")
        #os.chdir(workdir)
    else:
        print(f"Creating output folder: {workdir}")
        os.makedirs(workdir)
        # print(f"moving into output folder")
        # os.chdir(workdir)

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
        ccList(success, workdir)

    elif args.grenade == False:
        ccList(abs_file_list, workdir)
    
    #check correlation file and use abs path
    if os.path.isfile(f"{workdir}/ccClusterLog.txt"):
        correlationFile=f"{workdir}/ccClusterLog.txt" # it should be already an abs path
        print(f"use correlation file: {correlationFile}")
    else:
        print(f"Correlation file does not exist, please check")

    #set up the job
    CC = Clustering(correlationFile, workdir)
    Tree = CC.avgTree()
    etiquets = CC.createLabels()

    #check threshold
    if args.threshold == None:
        threshold = CC.thrEstimation()
    else:
        threshold = args.threshold

    #check file type
    CC.checkMultiplicity(threshold)
    fileType = CC.inputType()
    if fileType=="HKL":
        #prepare and run XSCALE job
        if args.reference_HKL == None:
            xscale_checker, xscale_path = CC.prepareXSCALE(anomlous, threshold, clusterList=args.clusters)
        elif os.path.isfile(args.reference_HKL):
            xscale_checker, xscale_path = CC.prepareXSCALE(anomlous, threshold, clusterList=args.clusters, refHKL=args.reference_HKL)
        else:
            xscale_checker, xscale_path = CC.prepareXSCALE(anomlous,threshold, clusterList=args.clusters)
        if xscale_checker == True:
            CC.scaleAndMerge(anomlous, threshold, xscale_path)
            #get jason from XSCALE
            CC.flatClusterPrinter(threshold, etiquets, xscale_path)

        #prepare and run Pointless
        pointless_checker, pointless_path = CC.preparePointless(anomlous, threshold, clusterList=args.clusters)
        if pointless_checker == True:
            CC.pointlessRun(anomlous, threshold, pointless_path)
            #prepare and run aimless
            CC.aimlessRun(anomlous, threshold, pointless_path)

        #CC.passOInfoToGA(threshold, etiquets, anomlous)
    elif fileType=="mtz":
        #prepare and run Pointless
        pointless_checker, pointless_path = CC.preparePointless(anomlous,threshold, clusterList=args.clusters)
        if pointless_checker == True:
            CC.pointlessRun(anomlous, threshold, pointless_path)
            #prepare and run aimless
            CC.aimlessRun(anomlous, threshold, pointless_path)
            CC.flatClusterPrinter(threshold, etiquets, pointless_path)
    else:
        print(f"Unknown input file format, please check the distance file: {correlationFile}")
    


if __name__ =='__main__' :
    main()

    
