#! /usr/bin/env python3
from __future__ import print_function, absolute_import

__author__ = "Gianluca Santoni"
__copyright__ = "Copyright 2015-2019"
__credits__ = ["Gianluca Santoni, Alexander Popov"]
__license__ = ""
__version__ = "1.0"
__maintainer__ = "Tiankun Zhou"
__email__ = "tiankun.zhou@esrf.fr"
__status__ = "Beta"


import matplotlib.pyplot as plt
import sys
import os
import subprocess
#import ccCluster classes
from .clustering import Clustering
# Insert parse  to change the file path from command line
import argparse
import textwrap
#Startup message

print(r"""ccCluster - HCA for protein crystallography
G. Santoni and A. Popov, 2015-2019
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

def process_args():
    input_args = argparse.ArgumentParser()

    input_args.add_argument("-i", "--dist_file",
    type=str,
    help="Distance file from ccCalc module"
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

    input_args.add_argument("-p", "--process",
    action="store_true",
    default=False,
    help="Launch program in shell mode. Need to specify the threshold value"
    )

    input_args.add_argument("-c", "--count",
    action="store_true",
    default=False,
    help="Counts datasets in the biggest cluster and exit"
    )

    input_args.add_argument("-e", "--estimation",
    action="store_true",
    default=False,
    help="Tries to guess an optimal threshold value"
    )

    input_args.add_argument("-ref", "--reference_HKL",
    type = str,
    help="Give an optional reference HKL file for XSCALE merging, recommend to give absolute path"
    )

    #save args
    args = input_args.parse_args()

    #set working folder
    if args.output_dir == None:
        args.output_dir = os.getcwd()

    #check whether processing type is there
    if not args.process and not args.count and not args.estimation:
        print(textwrap.dedent(f"""\
                                    Did not select how to process the data
                                    use -p for processing
                                    -c to calculate the biggest cluster
                                    and -e to estimate threshold"""
                            ))
        exit()

    #return args
    return args


#Main part of the program
#with the different options, we can chose
# to process through the shell,
#count the multiplicity of the highest cluster
def main():
    args = process_args()

    #set up abs path for ccClusteringlog.txt
    if os.path.isfile(args.dist_file):
        correlationFile = os.path.abspath(args.dist_file)
        print(f"abs_path_list: {correlationFile}")
    else:
        print(f"distance file does not exist, please check: args.dist_file")

    #set up anomalous
    if args.anomalous_scattering == False:
        anomlous = "no_ano"
    else:
        anomlous = "ano"

    #create working dir
    workdir = os.path.abspath(f"{args.output_dir}")
    print(f"ccCluster output directory would be : {workdir}")
    if os.path.isdir(workdir):
        print(f"Output folder exists: {workdir}")
        #print(f"moving into output folder")
        #os.chdir(workdir)
    else:
        print(f"Creating output folder: {workdir}")
        os.makedirs(workdir)
        #print(f"moving into output folder")
        #os.chdir(workdir)

    #set up the job
    CC = Clustering(correlationFile, workdir)
    Tree = CC.avgTree()
    etiquets=CC.createLabels()

    #check threshold
    if args.threshold == None:
        threshold = CC.thrEstimation()
    else:
        threshold = args.threshold

    #check file type
    if args.process:
        CC.checkMultiplicity(threshold)
        fileType = CC.inputType()
        if fileType=="HKL":
            #prepare and run XSCALE job
            if args.reference_HKL == None:
                xscale_checker, xscale_path = CC.prepareXSCALE(anomlous,threshold)
            elif os.path.isfile(args.reference_HKL):
                xscale_checker, xscale_path = CC.prepareXSCALE(anomlous,threshold, refHKL=args.reference_HKL)
            else:
                xscale_checker, xscale_path = CC.prepareXSCALE(anomlous,threshold)
            if xscale_checker == True:
                CC.scaleAndMerge(anomlous, threshold, xscale_path)
                #get jason from XSCALE
                CC.flatClusterPrinter(threshold, etiquets, anomlous, xscale_path)

            #prepare and run Pointless
            pointless_checker, pointless_path = CC.preparePointless(anomlous, threshold)
            if pointless_checker == True:
                CC.pointlessRun(anomlous, threshold, pointless_path)
                #prepare and run aimless
                CC.aimlessRun(anomlous, threshold, pointless_path)

        #CC.passOInfoToGA(threshold, etiquets, anomlous)
        elif fileType=="mtz":
            #prepare and run Pointless
            pointless_checker, pointless_path = CC.preparePointless(anomlous,threshold)
            if pointless_checker == True:
                CC.pointlessRun(anomlous, threshold, pointless_path)
                #prepare and run aimless
                CC.aimlessRun(anomlous, threshold, pointless_path)
                CC.flatClusterPrinter(threshold, etiquets, anomlous, pointless_path)
        else:
            print(f"Unknown input file format, please check the distance file: {correlationFile}")

    #check datasets in the biggest cluster and exit
    elif args.count:
        CC.checkMultiplicity(threshold)

    # check estimated threshold
    elif args.estimation:
        a = CC.thrEstimation()
        print(f"Estimated threshold is {a}, you can input another value use -t if needed")


"""
    #Old:
    args = process_args()

    #Suggest to run ccCalc if no correlation file is provided
    if args.DISTfile is None:
        print('no inputs specified, please run ccCalc before')
    else:
        correlationFile=os.path.abspath(args.DISTfile)

    CC = Clustering(correlationFile)
    Tree = CC.avgTree()
    etiquets=CC.createLabels()
    threshold = CC.thrEstimation()
    fileType = CC.inputType()

    if args.threshold:
        threshold = args.threshold
    else:
        threshold= CC.thrEstimation()
    if args.shell:
        CC.checkMultiplicity(threshold)
        if fileType=="HKL":
            CC.prepareXSCALE('ano',threshold)
            CC.scaleAndMerge('ano',threshold)
            CC.flatClusterPrinter(threshold, etiquets, 'ano')
            CC.passOInfoToGA(threshold, etiquets, 'ano')
        elif fileType=="mtz":
            CC.preparePointless('ano',threshold)
            CC.pointlessRun('ano',threshold)
            CC.flatClusterPrinter(threshold, etiquets, 'ano')
        else:            print("Unknown input file format.")
    elif args.count:
        CC.checkMultiplicity(threshold)
    elif args.est:
        a = CC.thrEstimation()
        print(a)
"""

if __name__== '__main__':
    main()
