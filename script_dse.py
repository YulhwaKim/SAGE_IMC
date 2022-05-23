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
    args.basefolder = args.net + '_Simul' 

    print("----------args----------")
    print(args)

    return args

def scanNumHierarchy(args):

    # other design params for IC and BU (hBus, SRAM buffer)
    designParams_ic_bu = [3, 128, 3, 128, 0, 128, 128]
    designParams_ic_bu_top = [3, 128, 3, 128, 1, 128, 0]

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
    for num_hierarchy in range(1, args.max_numHierarchy+1):
        # get arch design params (variable: numSubObjects)
        numHObj_list = get_numSubObject_list(num_hierarchy)
        num_design = len(numHObj_list)

        for designIdx in range(0, num_design):
            # get design 
            numHObj = numHObj_list[designIdx]

            # generate designParams
            designParams = []
            for hIdx in range(0, num_hierarchy-1):
                designParams_hObj = numHObj[hIdx] + designParams_ic_bu
                designParams.append(designParams_hObj)
            designParams.append(numHObj[num_hierarchy-1] + designParams_ic_bu_top)
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

def get_numSubObject_list(num_hierarchy):
    # get arch design params (variable: numSubObjects)
    list_numSubObjectRow = [ i for i in range(1,10) ]
    list_numSubObjectCol = [ i for i in range(1,9) ]

    numHObj_list = []

    for h in range(0, num_hierarchy-1):
        h_numHObj_list = numHObj_list
        numHObj_list = []
        for numSubObjectRow in list_numSubObjectRow:
            for numSubObjectCol in list_numSubObjectRow:
                if not ((numSubObjectRow == 1) & (numSubObjectCol ==1)):
                    # get numSubObject design
                    numSubObject = [numSubObjectRow, numSubObjectCol]
                    # update numHObj_list
                    if ( len(h_numHObj_list) == 0):
                        numHObj_list.append([numSubObject])
                    else:
                        for i in range(0, len(h_numHObj_list)):
                            tmp_numHObj = h_numHObj_list[i].copy()
                            tmp_numHObj.append(numSubObject)
                            numHObj_list.append(tmp_numHObj)

    # set the numSubObject of the top object as (1,1) as it would be tuned by the mapper
    if ( len(numHObj_list) == 0):
        numHObj_list.append([[10000, 10000]])
    else:
        for i in range(0, len(numHObj_list)):
            numHObj_list[i].append([10000,10000])

    #print(numHObj_list)
    return numHObj_list

if __name__ == '__main__':
    args = read_arguments()
    scanNumHierarchy(args)
