/*******************************************************************************
* Copyright (c) 2015-2017
* School of Electrical, Computer and Energy Engineering, Arizona State University
* PI: Prof. Shimeng Yu
* All rights reserved.
* 
* This source code is part of NeuroSim - a device-circuit-algorithm framework to benchmark 
* neuro-inspired architectures with synaptic devices(e.g., SRAM and emerging non-volatile memory). 
* Copyright of the model is maintained by the developers, and the model is distributed under 
* the terms of the Creative Commons Attribution-NonCommercial 4.0 International Public License 
* http://creativecommons.org/licenses/by-nc/4.0/legalcode.
* The source code is free and you can redistribute and/or modify it
* by providing that the following conditions are met:
* 
*  1) Redistributions of source code must retain the above copyright notice,
*     this list of conditions and the following disclaimer.
* 
*  2) Redistributions in binary form must reproduce the above copyright notice,
*     this list of conditions and the following disclaimer in the documentation
*     and/or other materials provided with the distribution.
* 
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
* 
* Developer list: 
*   Pai-Yu Chen	    Email: pchen72 at asu dot edu 
*                    
*   Xiaochen Peng   Email: xpeng15 at asu dot edu
********************************************************************************/

#include <cstdio>
#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <chrono>
#include <algorithm>
#include "constant.h"
#include "formula.h"
#include "HierarchyRoot.h"
#include "HierarchyObject.h"
#include "Param.h"
#include "Definition.h"
#include "util.h"
#include "NetworkScheduler.h"

using namespace std;


int main(int argc, char * argv[]) {

	auto start = chrono::high_resolution_clock::now();
	
	gen.seed(0);

    /* get architecture information */
	vector<vector<double>> designArch;
    vector<vector<double>> networkStructure;
	designArch = readCSV(argv[1]);
    networkStructure = readCSV(argv[2]);

    int numHierarchy = designArch.size();
    int numLayer = networkStructure.size();

	// define weight/input/memory precision from wrapper
	param->synapseBit = atoi(argv[3]);              // precision of synapse weight
	param->numBitInput = atoi(argv[4]);             // precision of input neural activation
	if (param->cellBit > param->synapseBit) {
		cout << "ERROR!: Memory precision is even higher than synapse precision, please modify 'cellBit' in Param.cpp!" << endl;
		param->cellBit = param->synapseBit;
	}

    int scheduler_type = atoi(argv[5]);
    int archIdx = atoi(argv[6]);
    string basefolder = argv[7];

    double numComputation = 0;
    for (int i=0; i<networkStructure.size(); i++) {
        //numComputation += 2*( networkStructure[i][0] * networkStructure[i][1] * networkStructure[i][2] 
        //                        * networkStructure[i][3] * networkStructure[i][4] * networkStructure[i][5] )
        //                  * param->synapseBit * param->numBitInput;
        numComputation += 2*( networkStructure[i][0] * networkStructure[i][1] * networkStructure[i][2] 
                                * networkStructure[i][3] * networkStructure[i][4] * networkStructure[i][5] );
    }
	
	/*** initialize operationMode as default ***/
	param->conventionalParallel = 0;
	param->conventionalSequential = 0;
	param->BNNparallelMode = 0;                // parallel BNN
	param->BNNsequentialMode = 0;              // sequential BNN
	param->XNORsequentialMode = 0;           // Use several multi-bit RRAM as one synapse
	param->XNORparallelMode = 0;         // Use several multi-bit RRAM as one synapse
	switch(param->operationmode) {
		case 6:	    param->XNORparallelMode = 1;               break;     
		case 5:	    param->XNORsequentialMode = 1;             break;     
		case 4:	    param->BNNparallelMode = 1;                break;     
		case 3:	    param->BNNsequentialMode = 1;              break;    
		case 2:	    param->conventionalParallel = 1;           break;     
		case 1:	    param->conventionalSequential = 1;         break;    
		case -1:	break;
		default:	exit(-1);
	}
	
	param->numColPerSynapse = ceil((double)param->synapseBit/(double)param->cellBit); 

    /* Architecture Design Initialization (Initialization include area calculation) */
    HierarchyRoot *hRoot;
    HierarchyObject *prevObject;
    HierarchyObject *lastObject;
    vector<HierarchyObject*> hObjectVector;
    printf("---------- Start Object Initialization ----------\n");
    for (int h=0; h < (numHierarchy + 1); h++) {
        printf("Initialize %d-level object\n",h);
        if ( h==0 ) {
            hRoot = new HierarchyRoot(inputParameter, tech, cell);
            hRoot->Initialize();
        } else {
            lastObject = new HierarchyObject(inputParameter, tech, cell, h, hRoot, prevObject, designArch[h-1]);
            lastObject->Initialize(param->clkFreq);

            hObjectVector.push_back(lastObject);
            prevObject = lastObject;
        }
    }
    printf("---------- Finish Object Initialization ----------\n\n");

    /* Architecture CLK Period Calculation */
    double clkPeriod;
    hRoot->CalculateLatency(true, (double)hRoot->numCol, (double)hRoot->numRow, 1, 1, 1, &clkPeriod);
    if ( param->clkFreq > 1/clkPeriod ) {
        param->clkFreq = 1/clkPeriod;
    }

    /* Get Architecture leakage information */
    double chipLeakage = lastObject->leakage; 

    /* Get Architecture area information */
    vector<double> chipAreaVector = lastObject->areaVector;
    
    vector<vector<double>> chipAreaVector2;
    for (int h=0; h < numHierarchy; h++) {
        auto tmpObject = hObjectVector.at(h);
        chipAreaVector2.push_back(tmpObject->areaVector2);
    }

    /* Network Scheduling */
    NetworkScheduler *networkScheduler = new NetworkScheduler();
    networkScheduler->Initialize(networkStructure, hObjectVector[numHierarchy-1]);
    vector<vector<double>> networkInfoRead = networkScheduler->Scheduling(scheduler_type);

    /* Architecture latency Measurement */
    vector<vector<double>> networkLatencyVector, networkEnergyVector;
    vector<vector<double>> layerLatencyVector2, layerEnergyVector2, networkLatencyVector2, networkEnergyVector2;
    vector<double> chipLatencyVector, chipEnergyVector;

    networkScheduler->CalculatePerformance(&networkLatencyVector, &networkEnergyVector, 
                                           &layerLatencyVector2, &layerEnergyVector2,
                                           &networkLatencyVector2, &networkEnergyVector2, 
                                           hRoot, hObjectVector, networkInfoRead);

    double clkPeriod_ns = clkPeriod * 1e9;

    /* Mergy latency / energy info */
    for (int layerIdx=0; layerIdx < networkLatencyVector.size(); layerIdx++) {

        vector<double> latencyVector = networkLatencyVector[layerIdx];
        vector<double> readDynamicEnergyVector = networkEnergyVector[layerIdx];

        // merge latency / energy info
        if ( chipLatencyVector.size() == 0 ) {
            chipLatencyVector.assign(latencyVector.begin(), latencyVector.end());
            chipEnergyVector.assign(readDynamicEnergyVector.begin(), readDynamicEnergyVector.end());
        } else {
            for (int latencyIdx=0; latencyIdx < latencyVector.size(); latencyIdx++) {
                chipLatencyVector[latencyIdx] += latencyVector[latencyIdx];
            }
            for (int energyIdx=0; energyIdx < readDynamicEnergyVector.size(); energyIdx++) {
                chipEnergyVector[energyIdx] += readDynamicEnergyVector[energyIdx];
            }
        }

    }
    double chipLeakageEnergy = chipLeakage * chipLatencyVector[0] * clkPeriod;   

    ///* Print Performance (latency/power) Information (breakdown type1) */
    //printLatencyVector(&networkLatencyVector2, chipLatencyVector[0], clkPeriod_ns, numHierarchy, numLayer, true, false);
    //printEnergyVector(&networkEnergyVector2, chipEnergyVector[0], numHierarchy, numLayer, true, false);

    //printAreaVector2(&chipAreaVector2, chipAreaVector[0], numHierarchy, &hObjectVector);
    //
    ///* Print Summary */
    //printf("\n-------------------- Summary --------------------\n");

    //cout << "Chip clock period is: " << clkPeriod_ns << "ns" << endl;
    //printf("%-20s %15.4e ns\n", "Chip clock period", clkPeriod_ns);
    //printf("%-20s %15.4e uW\n", "Chip leakagePower", chipLeakage*1e6);
    //printf("%-20s %15.4e pJ\n", "Chip leakageEnergy", chipLeakageEnergy*1e12);

    //printf("---------- chipArea Breakdown ----------\n");
    //printf("%-20s %15.4e um^2\n", "Chip Area", chipAreaVector[0]*1e12);
    //printf("%-20s %15.4e um^2 (%10.2f%%)\n", "SubArray", chipAreaVector[1]*1e12, chipAreaVector[1] / chipAreaVector[0] * 100);
    //printf("%-20s %15.4e um^2 (%10.2f%%)\n", "ADC", chipAreaVector[2]*1e12, chipAreaVector[2] / chipAreaVector[0] * 100);
    //printf("%-20s %15.4e um^2 (%10.2f%%)\n", "Accumulation", chipAreaVector[3]*1e12, chipAreaVector[3] / chipAreaVector[0] * 100);
    //printf("%-20s %15.4e um^2 (%10.2f%%)\n", "Buffer", chipAreaVector[4]*1e12, chipAreaVector[4] / chipAreaVector[0] * 100);
    //printf("%-20s %15.4e um^2 (%10.2f%%)\n", "IC", chipAreaVector[5]*1e12, chipAreaVector[5] / chipAreaVector[0] * 100);
    //printf("%-20s %15.4e um^2 (%10.2f%%)\n", "Other", chipAreaVector[6]*1e12, chipAreaVector[6] / chipAreaVector[0] * 100);


    //printf("---------- readLatency Breakdown (per image) ----------\n");
    //printf("%-20s %15.4e ns                [%5.1f cycle]\n", "Chip Latency", 
    //                                                        chipLatencyVector[0]*clkPeriod_ns, 
    //                                                        chipLatencyVector[0]);
    //printf("%-20s %15.4e ns (%10.2f%%)  [%5.1f cycle]\n", "SubArray", 
    //                                                    chipLatencyVector[1]*clkPeriod_ns, 
    //                                                    chipLatencyVector[1] / chipLatencyVector[0] * 100, 
    //                                                    chipLatencyVector[1]);
    //printf("%-20s %15.4e ns (%10.2f%%)  [%5.1f cycle]\n", "Accumulation", 
    //                                                    chipLatencyVector[2]*clkPeriod_ns, 
    //                                                    chipLatencyVector[2] / chipLatencyVector[0] * 100, 
    //                                                    chipLatencyVector[2]);
    //printf("%-20s %15.4e ns (%10.2f%%)  [%5.1f cycle]\n", "Buffer", 
    //                                                    chipLatencyVector[3]*clkPeriod_ns, 
    //                                                    chipLatencyVector[3] / chipLatencyVector[0] * 100, 
    //                                                    chipLatencyVector[3]);
    //printf("%-20s %15.4e ns (%10.2f%%)  [%5.1f cycle]\n", "IC", 
    //                                                    chipLatencyVector[4]*clkPeriod_ns, 
    //                                                    chipLatencyVector[4] / chipLatencyVector[0] * 100, 
    //                                                    chipLatencyVector[4]);
    //printf("%-20s %15.4e ns (%10.2f%%)  [%5.1f cycle]\n", "Other", 
    //                                                    chipLatencyVector[5]*clkPeriod_ns, 
    //                                                    chipLatencyVector[5] / chipLatencyVector[0] * 100, 
    //                                                    chipLatencyVector[5]);
    //printf("(Sum of the latency components exceed the total latency)\n");

    //printf("---------- readDynamicEnergy Breakdown (per image) ----------\n");
    //printf("%-20s %15.4e pJ\n", "Chip DynamicEnergy", chipEnergyVector[0]*1e12);
    //printf("%-20s %15.4e pJ (%10.2f%%)\n", "SubArray", 
    //                                            chipEnergyVector[1]*1e12, chipEnergyVector[1] / chipEnergyVector[0] * 100);
    //printf("%-20s %15.4e pJ (%10.2f%%)\n", "ADC", 
    //                                            chipEnergyVector[2]*1e12, chipEnergyVector[2] / chipEnergyVector[0] * 100);
    //printf("%-20s %15.4e pJ (%10.2f%%)\n", "Accumulation",
    //                                            chipEnergyVector[3]*1e12, chipEnergyVector[3] / chipEnergyVector[0] * 100);
    //printf("%-20s %15.4e pJ (%10.2f%%)\n", "Buffer",
    //                                            chipEnergyVector[4]*1e12, chipEnergyVector[4] / chipEnergyVector[0] * 100);
    //printf("%-20s %15.4e pJ (%10.2f%%)\n", "IC",
    //                                            chipEnergyVector[5]*1e12, chipEnergyVector[5] / chipEnergyVector[0] * 100);
    //printf("%-20s %15.4e pJ (%10.2f%%)\n", "Other",
    //                                            chipEnergyVector[6]*1e12, chipEnergyVector[6] / chipEnergyVector[0] * 100);
    //printf("(Accumulation Circuits - subarray level: adders, shiftAdds; PE/Tile/Global level: accumulation units)\n");
    //printf("(Other Peripheries     - pooling and activation units)\n");


    printf("---------- [ArchIdx: %8d] Performance ----------\n", archIdx);
    double topsw = numComputation/(chipEnergyVector[0]+chipLeakageEnergy)/1e12;
    double tops = numComputation/(chipLatencyVector[0]*clkPeriod)/1e12;
    printf("%-20s %15.4e mm^2\n", "Chip Area", chipAreaVector[0]*1e12/1e6);
    printf("%-20s %15.4f TOPS/W\n", "Energy Efficiency", topsw);
    printf("%-20s %15.4f TOPS\n", "Throughput", tops);


    printf("[START] Saving Simulation Results to CSV file \n");

    /* NOTE: indicators for data logging */
    string indicator_header;
    string indicator;

    //// we assume single layer simulation
    //if ( networkStructure.size() > 1 ) {
    //    cerr << "Indicator does not work for multi-layer simulation, but we got multiple layers!" << endl;
    //    exit(-1);
    //}
    indicator_header = "IC,OC,Wbit,Abit,numHierarchy,busType";
    indicator = to_string((int)networkStructure[0][2]) + "," + to_string((int)networkStructure[0][5]) + ","
                + to_string(param->synapseBit) + "," + to_string(param->numBitInput) + "," + to_string(numHierarchy+1);

    // get busType
    int busType = 0; // 0 - bus, 1 - sys1, 2 - sys2
    for ( int h=0; h < numHierarchy; h++ ) {
        if ( hObjectVector[h]->interConnect->inType == 1 ) {
            busType = h + 1;
        }
    }
    indicator += "," + to_string(busType);


    string filename = basefolder + "/performanceHObj/performanceHObj_" + to_string(archIdx) + ".csv";
    savePerformanceVector2(0, filename, indicator_header, indicator, &networkLatencyVector2, clkPeriod_ns, numHierarchy, 1);
    savePerformanceVector2(1, filename, indicator_header, indicator, &networkEnergyVector2, 1e12, numHierarchy, 1);
    vector<vector<double>> newAreaVector2 = updateAreaVector2(&chipAreaVector2, chipAreaVector[0], numHierarchy, &hObjectVector);
    savePerformanceVector2(2, filename, indicator_header, indicator, &newAreaVector2, 1e12, numHierarchy, 0);

    filename = basefolder + "/performanceChip/performanceChip_" + to_string(archIdx) + ".csv";
    savePerformanceVector(0, filename, indicator_header, indicator, &chipLatencyVector, clkPeriod_ns);
    savePerformanceVector(1, filename, indicator_header, indicator, &chipEnergyVector, 1e12);
    savePerformanceVector(2, filename, indicator_header, indicator, &chipAreaVector, 1e12);

    filename = basefolder + "/energy.csv";
    savePerformanceMetric(archIdx, filename, numHierarchy, chipEnergyVector[0]*1e12);

    filename = basefolder + "/energy_with_leakage.csv";
    savePerformanceMetric(archIdx, filename, numHierarchy, (chipEnergyVector[0]+chipLeakageEnergy)*1e12);

    filename = basefolder + "/latency.csv";
    savePerformanceMetric(archIdx, filename, numHierarchy, (chipLatencyVector[0]*clkPeriod)*1e12);
	
    filename = basefolder + "/topsw.csv";
    savePerformanceMetric(archIdx, filename, numHierarchy, topsw);

    filename = basefolder + "/tops.csv";
    savePerformanceMetric(archIdx, filename, numHierarchy, tops);

    filename = basefolder + "/area.csv";
    savePerformanceMetric(archIdx, filename, numHierarchy, chipAreaVector[0]*1e12);


    printf("[FINISH] Saving Simulation Results to CSV file \n");

}
