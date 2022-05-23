#include "util.h"
#include <fstream>
#include <string>
#include <iostream>
#include <unistd.h>

vector<vector<double>> readCSV(const string &inputfile) {
    ifstream infile(inputfile.c_str());
    string inputline;
    string inputval;

    int numInRow=0, numInCol=0;
    if (!infile.good()) {
        cerr << "Error: Input file [" << inputfile << "] cannot be opened!!" << endl;
        exit(1);
    } else {
        while (getline(infile, inputline, '\n')) {
            numInRow++;
        }
        infile.clear();
        infile.seekg(0, ios::beg);
        if (getline(infile, inputline, '\n')) {
            istringstream iss (inputline);
            while (getline(iss, inputval, ',')) {
                numInCol++;
            }
        }
    }
    infile.clear();
    infile.seekg(0, ios::beg);

    vector<vector<double>> archInfo;
    for ( int row=0; row < numInRow; row++ ) {
        vector<double> archInfoRow;
        getline(infile, inputline, '\n');
        istringstream iss;
        iss.str(inputline);
        for ( int col=0; col < numInCol; col++ ) {
            while( getline(iss, inputval, ',') ) {
                istringstream fs;
                fs.str(inputval);   
                double f=0;
                fs >> f;
                archInfoRow.push_back(f);
            }
        }
        archInfo.push_back(archInfoRow);
    }
    infile.close();

    return archInfo;
}

vector<vector<int>> readCSVint(const string &inputfile) {
    ifstream infile(inputfile.c_str());
    string inputline;
    string inputval;

    int numInRow=0, numInCol=0;
    if (!infile.good()) {
        cerr << "Error: Input file [" << inputfile << "] cannot be opened!!" << endl;
        exit(1);
    } else {
        while (getline(infile, inputline, '\n')) {
            numInRow++;
        }
        infile.clear();
        infile.seekg(0, ios::beg);
        if (getline(infile, inputline, '\n')) {
            istringstream iss (inputline);
            while (getline(iss, inputval, ',')) {
                numInCol++;
            }
        }
    }
    infile.clear();
    infile.seekg(0, ios::beg);

    vector<vector<int>> readInfo;
    for ( int row=0; row < numInRow; row++ ) {
        vector<int> readInfoRow;
        getline(infile, inputline, '\n');
        istringstream iss;
        iss.str(inputline);
        for ( int col=0; col < numInCol; col++ ) {
            while( getline(iss, inputval, ',') ) {
                istringstream fs;
                fs.str(inputval);   
                int f=0;
                fs >> f;
                readInfoRow.push_back(f);
            }
        }
        readInfo.push_back(readInfoRow);
    }
    infile.close();

    return readInfo;
}

void saveIntVector2(const string &filename,
        const vector<vector<int>> *intVector2) {

    // open file
    fstream fout;
    fout.open(filename, ios::out | ios::trunc /*remove contents if file exist*/);

    // get Vector2 size (#vector = #row of output file)
    int numRow = intVector2->size();

    // write intVector2 data
    vector<int> tmpVector;
    int tmpData, numCol; 
    for ( int vecIdx=0; vecIdx < numRow; vecIdx++ ) {

        tmpVector = intVector2->at(vecIdx);
        int numCol = tmpVector.size(); // data of vector on Col dim

        // get & write components
        for ( int dataIdx=0; dataIdx < numCol - 1; dataIdx++ ) {
            tmpData = tmpVector[dataIdx];
            fout << to_string(tmpData) << ","; 
        }
        tmpData = tmpVector.back();
        fout << to_string(tmpData) << "\n"; 
    }

    // close the file
    fout.close();
   
}

void printIntVector2(const vector<vector<int>> *intVector2) {
    
    // get Vector2 size (#vector = #row of print)
    int numRow = intVector2->size();

    // print intVector2 data
    vector<int> tmpVector;
    int tmpData, numCol; 
    for ( int vecIdx=0; vecIdx < numRow; vecIdx++ ) {

        tmpVector = intVector2->at(vecIdx);
        int numCol = tmpVector.size(); // data of vector on Col dim

        // get & print components
        for ( int dataIdx=0; dataIdx < numCol - 1; dataIdx++ ) {
            tmpData = tmpVector[dataIdx];
            cout << tmpData << ", ";
        }
        tmpData = tmpVector.back();
        cout << tmpData << endl;
        
    }

}

// save performance metric for iterative simulation
void savePerformanceMetric(int archIdx, const string &filename, 
                        int numHierarchy, double performance) {
    // get names of components
    string components = "archIdx,numHierarchy,metric";
    int numComponents = 2; // num components of the file

    // check if the file already exist
    bool fileExist = ( access(filename.c_str(), F_OK) != -1 );

    // open csv file
    fstream fout;
    fout.open(filename, ios::out | ios::app);

    // write header if neccesary
    if ( !fileExist ) {
        fout << components << "\n";
    }

    // write indicator (archIdx, numHierarchy)
    fout << archIdx << "," << numHierarchy << ",";

    // write performance
    char charPerformance[15];
    sprintf(charPerformance, "%.4e", performance);
    fout << charPerformance << "\n";

    // close the file
    fout.close();
}


// performanceType: 0 - latency, 1 - dynamicEnergy, 2 - area
void savePerformanceVector(int performanceType, const string &filename,
        const string &indicator_header, const string &indicator,
        const vector<double> *performanceVector, double scalingFactor) { // chip performance
    // get names of components
    string components = "Total,SubArray,ADC,Accum,Buffer,InterConnect,Other";
    int numComponents = 7; // num components from performanceVector

    // check if the file already exist
    bool fileExist = ( access(filename.c_str(), F_OK) != -1 );

    // open csv file
    fstream fout;
    fout.open(filename, ios::out | ios::app);

    // write header if neccesary
    if ( !fileExist ) {
        fout << indicator_header << ",";
        fout << "performanceType" << ",";
        fout << components << "\n";
    }

    string localIndicator;
    // wirte performanceType
    switch (performanceType) {
        case 0:
            localIndicator = indicator + ",latency";
            break;
        case 1:
            localIndicator = indicator + ",dynamicE";
            break;
        case 2:
            localIndicator = indicator + ",area";
            break;
        default:
            localIndicator = indicator + ",Unknown";
            break;
    }
   
    // write local indicator
    fout << localIndicator << ",";

    // get & write components
    double tmpPerformance; 
    for ( int compIdx=0; compIdx < numComponents - 1; compIdx++ ) {
        if ( performanceType == 0 ) {
            if ( compIdx < 2 ) {
                tmpPerformance = performanceVector->at(compIdx);
            } else if ( compIdx == 2 ) {
                tmpPerformance = 0;
            } else {
                tmpPerformance = performanceVector->at(compIdx-1);
            }
        } else {
            tmpPerformance = performanceVector->at(compIdx);
        }
        char charPerformance[15];
        sprintf(charPerformance, "%.4e", tmpPerformance * scalingFactor);
        //fout << str(charPerformance) << ",";
        fout << charPerformance << ",";
    }
    tmpPerformance = performanceVector->back();
    char charPerformance[15];
    sprintf(charPerformance, "%.4e", tmpPerformance * scalingFactor);
    fout << charPerformance << "\n";

    // close the file
    fout.close();
}

// performanceType: 0 - latency, 1 - dynamicEnergy, 2 - area
void savePerformanceVector2(int performanceType, const string &filename,
        const string &indicator_header, const string &indicator,
        const vector<vector<double>> *performanceVector,
        double scalingFactor, int numHierarchy, int bias) { // support only hObjBreakdown
    // get names of components
    string components = "hObj,Total,SubObject,Accum,Buffer,InterConnect,Other";
    int numComponents = 6; // num components from performanceVector (exclude hObj)

    // check if the file already exist
    bool fileExist = ( access(filename.c_str(), F_OK) != -1 );

    // open csv file
    fstream fout;
    fout.open(filename, ios::out | ios::app);
    
    // write header if neccesary
    if ( !fileExist ) {
        fout << indicator_header << ",";
        fout << "performanceType" << ",";
        fout << components << "\n";
    }

    string localIndicator;
    // wirte performanceType
    switch (performanceType) {
        case 0:
            localIndicator = indicator + ",latency";
            break;
        case 1:
            localIndicator = indicator + ",dynamicE";
            break;
        case 2:
            localIndicator = indicator + ",area";
            break;
        default:
            localIndicator = indicator + ",Unknown";
            break;
    }
   
    // write performance data
    vector<double> tmpPerformanceVector;
    double tmpPerformance; 
    for ( int hIdx=0; hIdx < numHierarchy; hIdx++ ) {
        tmpPerformanceVector = performanceVector->at(hIdx);
        // write local indicator
        fout << localIndicator << ",";
        // get & write components
        fout << to_string(hIdx+1) << ","; // hObj
        for ( int compIdx=0; compIdx < numComponents - 1; compIdx++ ) {
            tmpPerformance = tmpPerformanceVector[compIdx + bias];
            char charPerformance[15];
            sprintf(charPerformance, "%.4e", tmpPerformance * scalingFactor);
            fout << charPerformance << ",";
        }
        tmpPerformance = tmpPerformanceVector.back();
        char charPerformance[15];
        sprintf(charPerformance, "%.4e", tmpPerformance * scalingFactor);
        fout << charPerformance << "\n";
    }

    // close the file
    fout.close();
}

void printLatencyVector(const vector<vector<double>> *latencyVector, 
                double totalLatency, double clkPeriod_ns,
                int numHierarchy, int numLayer,
                bool hObjectBreakdown, bool layerBreakdown) {
    // get names of components
    vector<string> components;
    if ( hObjectBreakdown ) {
        components.insert(components.begin(), {"Total", "SubObject", "Accum", "Buffer", "IC", "Other"});
    } else{
        components.insert(components.begin(), {"Total", "SubArray", "Accum", "Buffer", "IC", "Other"});
    }

    // organize latencyVector information
    vector<double> tmpLatencyVector;
    if ( hObjectBreakdown ) {
        if ( layerBreakdown ) {
            for ( int layerIdx=0; layerIdx < numLayer; layerIdx ++ ) {
                printf("\n------------------- [Layer %2d] --------------------\n", layerIdx);
                for ( int hIdx=0; hIdx < numHierarchy; hIdx++ ) {
                    tmpLatencyVector = latencyVector->at(layerIdx * numHierarchy + hIdx);
                    if ( layerIdx != (int)tmpLatencyVector[0] ) {
                        cerr << "[printLatencyVector] layerIdx miss match: " << layerIdx <<" " << tmpLatencyVector[0] << endl;
                        exit(-1);
                    }
                    printf("--------- [hObj %2d] readLatency Breakdown (per image) ----------\n", (int)tmpLatencyVector[1]);
                    printf("%-20s %15.4e ns {%10.2f%%} [%5.1f cycle]\n", components[0].c_str(),
                                tmpLatencyVector[2]*clkPeriod_ns, tmpLatencyVector[2] / totalLatency * 100, tmpLatencyVector[2]);
                    for ( int compIdx=1; compIdx < components.size(); compIdx++ ) {
                        printf("%-20s %15.4e ns (%10.2f%%) [%5.1f cycle]\n", components[compIdx].c_str(),
                            tmpLatencyVector[compIdx+2]*clkPeriod_ns, tmpLatencyVector[compIdx+2] / tmpLatencyVector[2] * 100, 
                            tmpLatencyVector[compIdx+2]);
                    }
                }
            }
        } else {
            for ( int hIdx=0; hIdx < numHierarchy; hIdx++ ) {
                tmpLatencyVector = latencyVector->at(hIdx);
                printf("--------- [hObj %2d] readLatency Breakdown (per image) ----------\n", (int)tmpLatencyVector[0]);
                printf("%-20s %15.4e ns {%10.2f%%} [%5.1f cycle]\n", components[0].c_str(),
                            tmpLatencyVector[1]*clkPeriod_ns, tmpLatencyVector[1] / totalLatency * 100, tmpLatencyVector[1]);
                for ( int compIdx=1; compIdx < components.size(); compIdx++ ) {
                    printf("%-20s %15.4e ns (%10.2f%%) [%5.1f cycle]\n", components[compIdx].c_str(),
                            tmpLatencyVector[compIdx+1]*clkPeriod_ns, tmpLatencyVector[compIdx+1] / tmpLatencyVector[1] * 100, 
                                        tmpLatencyVector[compIdx+1]);
                }
            }
        }
    } else {
        for ( int layerIdx=0; layerIdx < numLayer; layerIdx++ ) {
            tmpLatencyVector = latencyVector->at(layerIdx);
            printf("--------- [Layer %2d] readLatency Breakdown (per image) ----------\n", layerIdx);
            printf("%-20s %15.4e ns {%10.2f%%} [%5.1f cycle]\n", components[0].c_str(),
                                tmpLatencyVector[0]*clkPeriod_ns, tmpLatencyVector[0] / totalLatency * 100, tmpLatencyVector[0]);
            for ( int compIdx=1; compIdx < components.size(); compIdx++ ) {
                printf("%-20s %15.4e ns (%10.2f%%) [%5.1f cycle]\n", components[compIdx].c_str(),
                                    tmpLatencyVector[compIdx]*clkPeriod_ns, tmpLatencyVector[compIdx] / tmpLatencyVector[0] * 100, 
                                    tmpLatencyVector[compIdx]);
            }
        }
    }
}


void printEnergyVector(const vector<vector<double>> *energyVector, 
                double totalEnergy, 
                int numHierarchy, int numLayer,
                bool hObjectBreakdown, bool layerBreakdown) {
    // get names of components
    vector<string> components;
    if ( hObjectBreakdown ) {
        components.insert(components.begin(), {"Total", "SubObject", "Accum", "Buffer", "IC", "Other"});
    } else{
        components.insert(components.begin(), {"Total", "SubArray", "ADC", "Accum", "Buffer", "IC", "Other"});
    }

    // organize energyVector information
    vector<double> tmpEnergyVector;
    if ( hObjectBreakdown ) {
        if ( layerBreakdown ) {
            for ( int layerIdx=0; layerIdx < numLayer; layerIdx ++ ) {
                printf("\n------------------- [Layer %2d] --------------------\n", layerIdx);
                for ( int hIdx=0; hIdx < numHierarchy; hIdx++ ) {
                    tmpEnergyVector = energyVector->at(layerIdx * numHierarchy + hIdx);
                    if ( layerIdx != (int)tmpEnergyVector[0] ) {
                        cerr << "[printEnergyVector] layerIdx miss match: " << layerIdx << " " << tmpEnergyVector[0] << endl;
                        exit(-1);
                    }
                    printf("--------- [hObj %2d] readDynamicEnergy Breakdown (per image) ----------\n", (int)tmpEnergyVector[1]);
                    printf("%-20s %15.4e pJ {%10.2f%%}\n", components[0].c_str(),
                                        tmpEnergyVector[2]*1e12, tmpEnergyVector[2] / totalEnergy * 100);
                    for ( int compIdx=1; compIdx < components.size(); compIdx++ ) {
                        printf("%-20s %15.4e pJ (%10.2f%%)\n", components[compIdx].c_str(),
                                        tmpEnergyVector[compIdx+2]*1e12, tmpEnergyVector[compIdx+2] / tmpEnergyVector[2] * 100);
                    }
                }
            }
        } else {
            for ( int hIdx=0; hIdx < numHierarchy; hIdx++ ) {
                tmpEnergyVector = energyVector->at(hIdx);
                printf("--------- [hObj %2d] readDynamicEnergy Breakdown (per image) ----------\n", (int)tmpEnergyVector[0]);
                printf("%-20s %15.4e pJ {%10.2f%%}\n", components[0].c_str(),
                                    tmpEnergyVector[1]*1e12, tmpEnergyVector[1] / totalEnergy * 100);
                for ( int compIdx=1; compIdx < components.size(); compIdx++ ) {
                    printf("%-20s %15.4e pJ (%10.2f%%)\n", components[compIdx].c_str(),
                                        tmpEnergyVector[compIdx+1]*1e12, tmpEnergyVector[compIdx+1] / tmpEnergyVector[1] * 100);
                }
            }
        }
    } else {
        for ( int layerIdx=0; layerIdx < numLayer; layerIdx++ ) {
            tmpEnergyVector = energyVector->at(layerIdx);
            printf("--------- [Layer %2d] readDynamicEnergy Breakdown (per image) ----------\n", layerIdx);
            printf("%-20s %15.4e pJ {%10.2f%%}\n", components[0].c_str(),
                                tmpEnergyVector[0]*1e12, tmpEnergyVector[0] / totalEnergy * 100);
            for ( int compIdx=1; compIdx < components.size(); compIdx++ ) {
                printf("%-20s %15.4e pJ (%10.2f%%)\n", components[compIdx].c_str(),
                                    tmpEnergyVector[compIdx]*1e12, tmpEnergyVector[compIdx] / tmpEnergyVector[0] * 100);
            }
        }
    }
}

void printAreaVector2(const vector<vector<double>> *areaVector, 
                double totalArea, int numHierarchy,
                const vector<HierarchyObject*> *hObjectVector) {
    // get names of components
    vector<string> components;
    components.insert(components.begin(), {"Total", "SubObject", "Accum", "Buffer", "IC", "Other"});
    
    // get number of object
    vector<int> numSubObjectVector, numObjectVector;
    for ( int hIdx=0; hIdx < numHierarchy; hIdx++ ) {
       numSubObjectVector.push_back(hObjectVector->at(hIdx)->numSubObject); 
    }
    for ( int hIdx=0; hIdx < numHierarchy; hIdx++ ) {
        int objIdx = hIdx+1;
        int numObject = 1;
        while ( objIdx < numHierarchy ) {
            numObject *= numSubObjectVector[objIdx];
            objIdx++;
        }
       numObjectVector.push_back(numObject); 
    }

    // organize areaVector information
    vector<double> tmpAreaVector;
    for ( int hIdx=0; hIdx < numHierarchy; hIdx++ ) {
        tmpAreaVector = areaVector->at(hIdx);
        int numObject = numObjectVector[hIdx];
        printf("--------- [hObj %2d] chipArea Breakdown ----------\n", hIdx+1);
        printf("%-20s %15.4e um^2 {%10.2f%%}\n", components[0].c_str(),
                            tmpAreaVector[0]*1e12*numObject, tmpAreaVector[0]*numObject / totalArea * 100);
        for ( int compIdx=1; compIdx < components.size(); compIdx++ ) {
            printf("%-20s %15.4e um^2 (%10.2f%%)\n", components[compIdx].c_str(),
                                tmpAreaVector[compIdx]*numObject*1e12, tmpAreaVector[compIdx] / tmpAreaVector[0] * 100);
        }
    }
}

vector<vector<double>> updateAreaVector2(const vector<vector<double>> *areaVector, 
                double totalArea, int numHierarchy,
                const vector<HierarchyObject*> *hObjectVector) {
    // get names of components
    vector<string> components;
    components.insert(components.begin(), {"Total", "SubObject", "Accum", "Buffer", "IC", "Other"});
    
    // get number of object
    vector<int> numSubObjectVector, numObjectVector;
    for ( int hIdx=0; hIdx < numHierarchy; hIdx++ ) {
       numSubObjectVector.push_back(hObjectVector->at(hIdx)->numSubObject); 
    }
    for ( int hIdx=0; hIdx < numHierarchy; hIdx++ ) {
        int objIdx = hIdx+1;
        int numObject = 1;
        while ( objIdx < numHierarchy ) {
            numObject *= numSubObjectVector[objIdx];
            objIdx++;
        }
       numObjectVector.push_back(numObject); 
    }

    // organize areaVector information
    vector<vector<double>> newAreaVector;
    for ( int hIdx=0; hIdx < numHierarchy; hIdx++ ) {
        vector<double> tmpAreaVector = areaVector->at(hIdx);
        int numObject = numObjectVector[hIdx];
        for ( int compIdx=0; compIdx < components.size(); compIdx++ ) {
            tmpAreaVector[compIdx] = tmpAreaVector[compIdx] * numObject;
        }
        newAreaVector.push_back(tmpAreaVector);
    }
    return newAreaVector;
}
