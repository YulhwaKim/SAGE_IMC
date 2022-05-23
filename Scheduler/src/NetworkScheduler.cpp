#include <cmath>
#include <iostream>
#include <fstream>
#include <map>
#include "LayerScheduler.h"
#include "NetworkScheduler.h"
#include "Param.h"

extern Param *param;

NetworkScheduler::NetworkScheduler() {
    int hlevelMappingUnit = 2; // TILE
}

void NetworkScheduler::Initialize(const vector<vector<double>> _networkStructure,
                                 const HierarchyObject *_hTop) {
    // update network info
    networkStructure = _networkStructure;
    numLayer = networkStructure.size();
    // update hardware info
    hTop = _hTop;
    // initialize numUsedSubObject_top
    numUsedSubObject_top = 0;
}

vector<vector<double>> NetworkScheduler::Scheduling(int scheduler_type) {

    // clear scheduling result
    networkInfoRead.clear();
    // clear numUsedSubObject_top
    numUsedSubObject_top = 0;

    // define layer scheduler
    LayerScheduler *layerScheduler = new LayerScheduler();
    vector<vector<int>> idxOffsetVector;
    bool offsetObject = false;

    if ( scheduler_type == 0 ) {
        printf("non-compact mapping!\n");
        for (int layerIdx=0; layerIdx < numLayer; layerIdx++) {

            // initialzie layer scheduler
            vector<double> layerStructure = networkStructure[layerIdx];
            layerScheduler->Initialize(layerIdx, layerStructure, hlevelMappingUnit, idxOffsetVector, hTop);
        
            // layer scheduling 
            vector<vector<double>> layerInfoRead;
            layerInfoRead = layerScheduler->HObjectScheduling_00(hTop, 0, 0, 
                                                            layerScheduler->kH, layerScheduler->kW, layerScheduler->inC,
                                                            layerScheduler->outC * layerScheduler->numCellPerSynapse,
                                                            offsetObject, true);
            // get used #top subObject
            numUsedSubObject_top = MAX(numUsedSubObject_top, layerScheduler->numUsedSubObject_top);

            // get layer info
            networkInfoRead.insert(networkInfoRead.end(), layerInfoRead.begin(), layerInfoRead.end()); // update
            layerInfoRead.clear(); // clear updated result

            // get offset info
            idxOffsetVector.clear();
            if ( layerScheduler->nextIdxOffsetVector.size() > 0 ) {
                offsetObject = true;
                idxOffsetVector.assign(layerScheduler->nextIdxOffsetVector.begin(), layerScheduler->nextIdxOffsetVector.end());
                layerScheduler->nextIdxOffsetVector.clear();
            } else {
                offsetObject = false;
            }
        
        }
    } else if ( scheduler_type == 1 ) {
        printf("compact mapping!\n");
        for (int layerIdx=0; layerIdx < numLayer; layerIdx++) {

            //printf("layer %d\n", layerIdx);

            // initialzie layer scheduler
            vector<double> layerStructure = networkStructure[layerIdx];
            layerScheduler->Initialize(layerIdx, layerStructure, hlevelMappingUnit, idxOffsetVector, hTop);
        
            // layer scheduling 
            vector<vector<double>> layerInfoRead;
            layerInfoRead = layerScheduler->HObjectScheduling_01(hTop, 0, 0, -1, -1,
                                                            layerScheduler->kH, layerScheduler->kW, layerScheduler->inC,
                                                            layerScheduler->outC * layerScheduler->numCellPerSynapse);

            //// print scheduling data
            //if ( layerIdx < 3 ) {
            //    int numRow = layerInfoRead.size();
            //    printf("length layerInfoRead: %d\n", numRow);
            //    vector<double> tmpVector;
            //    int numCol; 
            //    double tmpData;
            //    for ( int vecIdx=0; vecIdx < numRow; vecIdx++ ) {

            //        tmpVector = layerInfoRead.at(vecIdx);
            //        int numCol = tmpVector.size(); // data of vector on Col dim

            //        // get & print components
            //        for ( int dataIdx=0; dataIdx < numCol - 1; dataIdx++ ) {
            //            tmpData = tmpVector[dataIdx];
            //            cout << tmpData << ", ";
            //        }
            //        tmpData = tmpVector.back();
            //        cout << tmpData << endl;
            //    }
            //}

            // get used #top subObject
            numUsedSubObject_top = MAX(numUsedSubObject_top, layerScheduler->numUsedSubObject_top);

            // get layer info
            networkInfoRead.insert(networkInfoRead.end(), layerInfoRead.begin(), layerInfoRead.end()); // update
            layerInfoRead.clear(); // clear updated result

            // get offset info
            idxOffsetVector.clear();
            if ( layerScheduler->nextIdxOffsetVector.size() > 0 ) {
                vector<vector<int>> tmp_idxOffsetVector;
                vector<int> tmp;
                tmp_idxOffsetVector.assign(layerScheduler->nextIdxOffsetVector.begin(), 
                                           layerScheduler->nextIdxOffsetVector.end());
                layerScheduler->nextIdxOffsetVector.clear();
                // reorder the idxOffstVector
                while ( tmp_idxOffsetVector.size() > 0 ) {
                    tmp = tmp_idxOffsetVector.back();
                    tmp_idxOffsetVector.pop_back();
                    idxOffsetVector.push_back(tmp);
                }
            }

            //// print next offset data
            //int numRow = idxOffsetVector.size();
            //vector<int> tmpVector;
            //int tmpData, numCol; 
            //for ( int vecIdx=0; vecIdx < numRow; vecIdx++ ) {

            //    tmpVector = idxOffsetVector.at(vecIdx);
            //    int numCol = tmpVector.size(); // data of vector on Col dim

            //    // get & print components
            //    for ( int dataIdx=0; dataIdx < numCol - 1; dataIdx++ ) {
            //        tmpData = tmpVector[dataIdx];
            //        cout << tmpData << ", ";
            //    }
            //    tmpData = tmpVector.back();
            //    cout << tmpData << endl;
            //    
            //}
        
        }
    }


    // write scheduling results to the csv file
    ofstream out("scheduling_result.csv");
    for (auto& row : networkInfoRead) {
        for (auto col : row) {
            out << col << ',';
        }
        out << '\n';
    }

    return networkInfoRead;

}

void NetworkScheduler::CalculatePerformance(
                        vector<vector<double>> *networkLatencyVector,
                        vector<vector<double>> *networkEnergyVector,
                        vector<vector<double>> *layerLatencyVector2,
                        vector<vector<double>> *layerEnergyVector2,
                        vector<vector<double>> *networkLatencyVector2,
                        vector<vector<double>> *networkEnergyVector2,
                        HierarchyRoot *hRoot, vector<HierarchyObject*> hObjectVector,
                        const vector<vector<double>> networkInfoRead) {

    // space for keeping network latency results
    networkLatencyVector->clear(); 
    networkEnergyVector->clear(); 
    layerLatencyVector2->clear();
    layerEnergyVector2->clear();
    networkLatencyVector2->clear();
    networkEnergyVector2->clear();

    // for keeping latency & energy information
    vector<int> hlevelHeap; // layer-wise info
    vector<vector<double>> latencyHeap, energyHeap; // layer-wise info
    vector<double> latencyVector, subLatencyVector, tmpLatencyVector; // hlevel-wise info
    vector<double> energyVector, subEnergyVector, tmpEnergyVector; // hlevel-wise info

    vector<double> latencyVector2, tmpLatencyVector2;
    vector<double> energyVector2, tmpEnergyVector2;
    map<int, vector<double>> layerLatencyDict2, layerEnergyDict2; // hlevel, latency/energy breakdown
    map<int, vector<double>> totalLatencyDict2, totalEnergyDict2; // hlevel, latency/energy breakdown

    // get required info
    int numInfoRead = networkInfoRead.size(); // number of calculation
    int layerIdx = (int)networkInfoRead[0][0]; // layer index of the first infoRead
    int hlevel = (int)networkInfoRead[0][1]; // hlevel of the first infoRead
    int layerHlevel = 0;

    // get latency information for each operation
    for ( int infoIdx=0; infoIdx < numInfoRead; infoIdx++ ) {

        //printf("infoIdx: %d\n", infoIdx);

        vector<double> infoRead = networkInfoRead[infoIdx];

        //printf("1\n");

        // check if moved to the other hierarchy
        if ( (int)infoRead[1] != hlevel ) {

            // clear sub latency/energy
            subLatencyVector.clear();
            subEnergyVector.clear();

            // update layer-wise & sub latency/energy info
            if ( (int)infoRead[1] == 0 ) { // go back to root
                // update latency/energy heap
                if ( (hlevelHeap.size() == 0) || (hlevelHeap.back() > hlevel) ) { // update new hlevel object (size==0 -> highest)
                    // push_back heap
                    hlevelHeap.push_back(hlevel);
                    latencyHeap.push_back(latencyVector);
                    energyHeap.push_back(energyVector);
                } else if ( hlevelHeap.back() == hlevel ) { // update the same hlevel object
                    // pop_back heap to update heap
                    tmpLatencyVector = latencyHeap.back();
                    tmpEnergyVector = energyHeap.back();
                    latencyHeap.pop_back();
                    energyHeap.pop_back();
                    // update vector
                    for ( int latencyIdx=0; latencyIdx < latencyVector.size(); latencyIdx++ ) { // update max latency
                        tmpLatencyVector[latencyIdx] = MAX( tmpLatencyVector[latencyIdx], latencyVector[latencyIdx] );
                    }
                    for ( int energyIdx=0; energyIdx < energyVector.size(); energyIdx++ ) { // addup energy
                        tmpEnergyVector[energyIdx] = tmpEnergyVector[energyIdx] + energyVector[energyIdx];
                    }
                    // update heap
                    latencyHeap.push_back(tmpLatencyVector);
                    energyHeap.push_back(tmpEnergyVector);
                } else { 
                    cerr << "[CalculateNetworkPerformance] Unexpected Pattern for Heap" << endl;
                    exit(-1);
                }
            } else { // go to the next hierarhcy
                if ( (hlevelHeap.size() == 0) || (hlevelHeap.back() > hlevel) ) { // local scanning (initial || hlevel=0)
                    // update sub
                    subLatencyVector.assign(latencyVector.begin(), latencyVector.end());
                    subEnergyVector.assign(energyVector.begin(), energyVector.end());
                } else if ( hlevelHeap.back() == hlevel ) { // merge hlevel info
                    // pop_back heap to generate subLatency/Energy
                    tmpLatencyVector = latencyHeap.back();
                    tmpEnergyVector = energyHeap.back();
                    hlevelHeap.pop_back();
                    latencyHeap.pop_back();
                    energyHeap.pop_back();
                    // update sub
                    for ( int latencyIdx=0; latencyIdx < latencyVector.size(); latencyIdx++ ) { // update max latency
                        double tmpLatency = MAX( tmpLatencyVector[latencyIdx], latencyVector[latencyIdx] );
                        subLatencyVector.push_back( tmpLatency );
                    }
                    for ( int energyIdx=0; energyIdx < energyVector.size(); energyIdx++ ) { // addup energy
                        double tmpEnergy = tmpEnergyVector[energyIdx] + energyVector[energyIdx];
                        subEnergyVector.push_back( tmpEnergy );
                    }
                } else {
                    cerr << "[CalculateNetworkPerformance] Unexpected Pattern for Heap" << endl;
                    exit(-1);
                }
            }

            // clear hlevel-wise vectors
            latencyVector.clear();
            energyVector.clear();

            // update breakdown type2
            if ( hlevel > 0 ) {
                // update latency
                if ( layerLatencyDict2.find(hlevel) == layerLatencyDict2.end() ) { // no hlevel info yet
                    layerLatencyDict2[hlevel] = latencyVector2;
                } else {
                    layerLatencyDict2[hlevel] = MAX(layerLatencyDict2[hlevel], latencyVector2);
                }
                // update energy
                if ( layerEnergyDict2.find(hlevel) == layerEnergyDict2.end() ) { // no hlevel info yet
                    layerEnergyDict2[hlevel] = energyVector2;
                } else {
                    for ( int energyIdx=0; energyIdx < energyVector2.size(); energyIdx++ ) { // addup energy
                        layerEnergyDict2[hlevel][energyIdx] = layerEnergyDict2[hlevel][energyIdx] + energyVector2[energyIdx];
                    }
                }
            }
            // clear temporary vectors
            latencyVector2.clear();
            energyVector2.clear();
            

            // move to the next hierarchy
            hlevel = infoRead[1];
        }

        //printf("2\n");

        // check if moved to next layer operation
        if ( (int)infoRead[0] != layerIdx ) {
            //printf("2-1\n");
            //printf("size: %d %d %d %d\n", networkLatencyVector->size(), latencyHeap.size(), networkEnergyVector->size(), energyHeap.size());
            // update latency/Energy info
            networkLatencyVector->push_back(latencyHeap.back());
            networkEnergyVector->push_back(energyHeap.back());

            //printf("2-2\n");
            // clear/initialize layer info
            latencyHeap.clear();
            energyHeap.clear();
            hlevelHeap.clear();

            //printf("2-3\n");
            // update breakdown type2
            // update latency2
            for ( auto it = layerLatencyDict2.begin(); it != layerLatencyDict2.end(); it++ ) {
                int key = it->first;
                vector<double> val = it->second;
                // update total latency
                if ( totalLatencyDict2.find(key) == totalLatencyDict2.end() ) { // no hlevel info yet
                    totalLatencyDict2[key] = val;
                } else {
                    for ( int latencyIdx=0; latencyIdx < val.size(); latencyIdx++ ) { // add latency
                        totalLatencyDict2[key][latencyIdx] += val[latencyIdx];
                    }
                }
                // update layer latency
                val.insert(val.begin(), (double)key);
                val.insert(val.begin(), (double)layerIdx);
                layerLatencyVector2->push_back(val);
            }
            layerLatencyDict2.clear();
            //printf("2-4\n");
            // update energy2
            for ( auto it = layerEnergyDict2.begin(); it != layerEnergyDict2.end(); it++ ) {
                int key = it->first;
                vector<double> val = it->second;
                // update total energy
                if ( totalEnergyDict2.find(key) == totalEnergyDict2.end() ) { // no hlevel info yet
                    totalEnergyDict2[key] = val;
                } else {
                    for ( int energyIdx=0; energyIdx < val.size(); energyIdx++ ) { // add energy
                        totalEnergyDict2[key][energyIdx] += val[energyIdx];
                    }
                }
                // update layer energy
                val.insert(val.begin(), (double)key);
                val.insert(val.begin(), (double)layerIdx);
                layerEnergyVector2->push_back(val);
            }
            layerEnergyDict2.clear();

            //printf("2-5\n");
            // move to next layer
            layerIdx = infoRead[0];
        }

        //printf("3\n");

        // get Performance info
        if ( hlevel == 0 ) {
            // get infoRead
            vector<double> infoReadHRoot {&infoRead[4], &infoRead[9]};
            // Calculate Latency
            hRoot->CalculateLatency(infoReadHRoot);
            tmpLatencyVector = hRoot->latencyVector;
            // Calculate Energy
            hRoot->CalculatePower(infoReadHRoot);
            tmpEnergyVector = hRoot->readDynamicEnergyVector;
        } else {
            HierarchyObject *hObject = hObjectVector[hlevel-1];
            if ( infoRead.size() == param->lengthInfoReadForIC ) { // get IC performance
                while ( infoRead.size() == param->lengthInfoReadForIC) {
                    vector<double> infoReadIC {&infoRead[4], &infoRead[10]};
                    // Calculate IC Latency
                    hObject->CalculateICLatency(infoReadIC); 
                    // Calculate IC Energy
                    hObject->CalculateICPower(infoReadIC);
                    // get next infoRead
                    infoIdx += 1; 
                    infoRead = networkInfoRead[infoIdx];
                }
            }
            // get info Read
            vector<double> infoReadHObject {&infoRead[9], &infoRead[25]};
            // Calculate Latency
            hObject->CalculateLatency(infoReadHObject, subLatencyVector);
            tmpLatencyVector = hObject->latencyVector;
            tmpLatencyVector2 = hObject->latencyVector2;
            // Calculate Energy
            hObject->CalculatePower(infoReadHObject, subEnergyVector);
            tmpEnergyVector = hObject->readDynamicEnergyVector;
            tmpEnergyVector2 = hObject->readDynamicEnergyVector2;
        }


        //printf("4\n");

        // update Performance for breakdown type1
        if ( latencyVector.size() == 0 ) { // hRoot/hObject lie in the same parent hObject is not calculated yet
            latencyVector.assign(tmpLatencyVector.begin(), tmpLatencyVector.end());
            energyVector.assign(tmpEnergyVector.begin(), tmpEnergyVector.end());
        } else {
            for ( int latencyIdx=0; latencyIdx < latencyVector.size(); latencyIdx++ ) { // update max latency
                latencyVector[latencyIdx] = MAX(latencyVector[latencyIdx], tmpLatencyVector[latencyIdx]);
            }
            for ( int energyIdx=0; energyIdx < energyVector.size(); energyIdx++ ) { // addup energy
                energyVector[energyIdx] = energyVector[energyIdx] + tmpEnergyVector[energyIdx];
            }
        }
        tmpLatencyVector.clear();
        tmpEnergyVector.clear();

        //printf("5\n");

        // update Performance for breakdown type2
        if ( latencyVector2.size() == 0 ) { // no hRoot/hObject lie in the same parent hObject is calculated yet
            latencyVector2.assign(tmpLatencyVector2.begin(), tmpLatencyVector2.end());
            energyVector2.assign(tmpEnergyVector2.begin(), tmpEnergyVector2.end());
        } else {
            for ( int latencyIdx=0; latencyIdx < latencyVector2.size(); latencyIdx++ ) { // update max latency
                latencyVector2[latencyIdx] = MAX(latencyVector2[latencyIdx], tmpLatencyVector2[latencyIdx]);
            }
            for ( int energyIdx=0; energyIdx < energyVector2.size(); energyIdx++ ) { // addup energy
                energyVector2[energyIdx] = energyVector2[energyIdx] + tmpEnergyVector2[energyIdx];
            }
        }
        tmpLatencyVector2.clear();
        tmpEnergyVector2.clear();

        //printf("6\n");
    }

    // update latency/energy heap
    if ( (hlevelHeap.size() == 0) || (hlevelHeap.back() > hlevel) ) { // update new hlevel object (size==0 -> highest)
        // push_back heap
        hlevelHeap.push_back(hlevel);
        latencyHeap.push_back(latencyVector);
        energyHeap.push_back(energyVector);
    } else if ( hlevelHeap.back() == hlevel ) { // update the same hlevel object
        // pop_back heap to update heap
        tmpLatencyVector = latencyHeap.back();
        tmpEnergyVector = energyHeap.back();
        latencyHeap.pop_back();
        energyHeap.pop_back();
        // update vector
        for ( int latencyIdx=0; latencyIdx < latencyVector.size(); latencyIdx++ ) { // update max latency
            tmpLatencyVector[latencyIdx] = MAX( tmpLatencyVector[latencyIdx], latencyVector[latencyIdx] );
        }
        for ( int energyIdx=0; energyIdx < energyVector.size(); energyIdx++ ) { // addup energy
            tmpEnergyVector[energyIdx] = tmpEnergyVector[energyIdx] + energyVector[energyIdx];
        }
        // update heap
        latencyHeap.push_back(tmpLatencyVector);
        energyHeap.push_back(tmpEnergyVector);
    } else { 
        cerr << "[CalculateNetworkPerformance] Unexpected Pattern for Heap" << endl;
        exit(-1);
    }
    // clear temporary vectors
    subLatencyVector.clear();
    subEnergyVector.clear();
    latencyVector.clear();
    energyVector.clear(); 

    // update breakdwon type2
    if ( hlevel > 0 ) {
        if ( layerLatencyDict2.find(hlevel) == layerLatencyDict2.end() ) { // no hlevel info yet
            layerLatencyDict2[hlevel] = latencyVector2;
            layerEnergyDict2[hlevel] = energyVector2;
        } else {
            layerLatencyDict2[hlevel] = MAX(layerLatencyDict2[hlevel], latencyVector2);
            for ( int energyIdx=0; energyIdx < energyVector2.size(); energyIdx++ ) { // addup energy
                layerEnergyDict2[hlevel][energyIdx] = layerEnergyDict2[hlevel][energyIdx] + energyVector2[energyIdx];
            }
        }
        latencyVector2.clear();
        energyVector2.clear();
    }

    /* finish by updating */
    // update latency/Energy info
    networkLatencyVector->push_back(latencyHeap.back());
    networkEnergyVector->push_back(energyHeap.back());

    // clear/initialize layer info
    latencyHeap.clear();
    energyHeap.clear();
    hlevelHeap.clear();

    // update breakdown type2
    // update latency2
    for ( auto it = layerLatencyDict2.begin(); it != layerLatencyDict2.end(); it++ ) {
        int key = it->first;
        vector<double> val = it->second;
        // update total latency
        if ( totalLatencyDict2.find(key) == totalLatencyDict2.end() ) { // no hlevel info yet
            totalLatencyDict2[key] = val;
        } else {
            for ( int latencyIdx=0; latencyIdx < val.size(); latencyIdx++ ) { // add latency
                totalLatencyDict2[key][latencyIdx] += val[latencyIdx];
            }
        }
        // update layer latency
        val.insert(val.begin(), (double)key);
        val.insert(val.begin(), (double)layerIdx);
        layerLatencyVector2->push_back(val);
    }
    layerLatencyDict2.clear();
    // update energy2
    for ( auto it = layerEnergyDict2.begin(); it != layerEnergyDict2.end(); it++ ) {
        int key = it->first;
        vector<double> val = it->second;
        // update total energy
        if ( totalEnergyDict2.find(key) == totalEnergyDict2.end() ) { // no hlevel info yet
            totalEnergyDict2[key] = val;
        } else {
            for ( int energyIdx=0; energyIdx < val.size(); energyIdx++ ) { // add energy
                totalEnergyDict2[key][energyIdx] += val[energyIdx];
            }
        }
        // update layer energy
        val.insert(val.begin(), (double)key);
        val.insert(val.begin(), (double)layerIdx);
        layerEnergyVector2->push_back(val);
    }
    layerEnergyDict2.clear();

    // update total latency2
    for ( auto it = totalLatencyDict2.begin(); it != totalLatencyDict2.end(); it++ ) {
        int key = it->first;
        vector<double> val = it->second;
        val.insert(val.begin(), (double)key);
        networkLatencyVector2->push_back(val);
    }
    totalLatencyDict2.clear();
    // update total energy2
    for ( auto it = totalEnergyDict2.begin(); it != totalEnergyDict2.end(); it++ ) {
        int key = it->first;
        vector<double> val = it->second;
        val.insert(val.begin(), (double)key);
        networkEnergyVector2->push_back(val);
    }
    totalEnergyDict2.clear();
}

