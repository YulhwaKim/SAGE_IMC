import argparse
import os
import csv
import math
import pandas as pd
#import matplotlib.pyplot as plt
#import seaborn as sns
import numpy as np

def read_arguments():
    parser = argparse.ArgumentParser(
                description = 'read pickle and plot the result')
    
    parser.add_argument('--net', default='Network_VGG8', type=str)
    parser.add_argument('--netfile', default='Network_VGG8.csv', type=str)
    parser.add_argument('--basefolder', default='Network_VGG8_Simul', type=str)
    parser.add_argument('--max-numHierarchy', default=3, type=int)
    parser.add_argument('--numRowCIMArray', default=128, type=int)
    parser.add_argument('--numColCIMArray', default=128, type=int)
    parser.add_argument('-w', '--wbits', default=4, type=int)
    parser.add_argument('-a', '--abits', default=4, type=int)
    parser.add_argument('--cellBit', default=2, type=int)
    parser.add_argument('--numCellPerSynpase', default=2, type=int)
    parser.add_argument('--compact-mapping', default=1, type=int)
    args = parser.parse_args()

    args.numCellPerSynapse = math.ceil( args.wbits / args.cellBit )
    args.netfile = args.net + '.csv' 
    args.basefolder = args.net + '_Simul_prev' 

    print("----------args----------")
    print(args)

    return args

def testConventionalArch(args):

    # other design params for IC and BU (hBus, SRAM buffer)
    designParams_NeuroSim = [[4, 4, 3, 128, 3, 128, 0, 128, 128], 
                             [2, 2, 3, 128, 3, 128, 0, 128, 128], 
                             [1000, 1000, 3, 128, 3, 128, 1, 128, 0]]
    designParams_ISAAC =    [[4, 2, 3, 128, 3, 128, 0, 128, 128], 
                             [3, 4, 3, 128, 3, 128, 0, 128, 128], 
                             [1000, 1000, 3, 128, 3, 128, 1, 128, 0]]
    designParams_PIMCA =    [[3, 6, 3, 128, 3, 128, 0, 128, 128], 
                             [1000, 1000, 3, 128, 3, 128, 1, 128, 0]]
    designParams_PUMA =     [[2, 1, 3, 128, 3, 128, 0, 128, 128], 
                             [8, 1, 3, 128, 3, 128, 0, 128, 128], 
                             [1000, 1000, 3, 128, 3, 128, 1, 128, 0]]
    designParams_list = [designParams_NeuroSim, designParams_ISAAC, designParams_PIMCA,
                         designParams_PUMA, designParams_CMP]

    # make basefolder (folder for gathering simulation data)
    if os.path.exists(args.basefolder):
        print(f"remove folder {args.basefolder}")
        os.system(f"rm -rf {args.basefolder}")
    os.makedirs(args.basefolder)

    # make folders to categorize data
    designParam_folder = os.path.join(args.basefolder, "designParam")
    designArch_folder = os.path.join(args.basefolder, "designArch")
    performanceChip_folder = os.path.join(args.basefolder, "performanceChip")
    performanceHObj_folder = os.path.join(args.basefolder, "performanceHObj")

    os.makedirs(designParam_folder)
    os.makedirs(designArch_folder)
    os.makedirs(performanceChip_folder)
    os.makedirs(performanceHObj_folder)

    design_counter = 0

    # simulate different architectures
    # get arch design params (variable: numSubObjects)
    num_design = len(designParams_list)

    for designIdx in range(0, num_design):
        # get designParams
        designParams = designParams_list[designIdx]
        # save designParams as designParams.csv file
        filename_designParam = os.path.join(designParam_folder, "designParam_{}.csv".format(design_counter))
        with open(filename_designParam, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerows(designParams)

        # generate architecture with design Params
        filename_designArch = os.path.join(designArch_folder, "designArch_{}.csv".format(design_counter))
        print(f"./arch_generator_for_net {filename_designParam} {args.netfile} {args.wbits} {args.abits} {args.compact_mapping} {filename_designArch}")
        os.system(f"./arch_generator_for_net {filename_designParam} {args.netfile} {args.wbits} {args.abits} {args.compact_mapping} {filename_designArch}")

        # do the simulation with generatred architecture
        print(f"./main_iter {filename_designArch} {args.netfile} {args.wbits} {args.abits} {args.compact_mapping} {design_counter} {args.basefolder}")
        os.system(f"./main_iter {filename_designArch} {args.netfile} {args.wbits} {args.abits} {args.compact_mapping} {design_counter} {args.basefolder}")
        
        design_counter += 1

if __name__ == '__main__':
    args = read_arguments()
    testConventionalArch(args)
