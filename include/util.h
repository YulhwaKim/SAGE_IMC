#ifndef NEUROHUB_UTIL_H_
#define NEUROHUB_UTIL_H_

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
#include "HierarchyObject.h"


using namespace std;

vector<vector<double>> readCSV(const string &inputfile);
vector<vector<int>> readCSVint(const string &inputfile);
void saveIntVector2(const string &filename,
        const vector<vector<int>> *intVector2);
void printIntVector2(const vector<vector<int>> *intVector2);

// save performance metric for iterative simulation
void savePerformanceMetric(int archIdx, const string &filename, 
                        int numHierarchy, double performance);

// performanceType: 0 - latency, 1 - dynamicEnergy, 2 - area
void savePerformanceVector(int performanceType, const string &filename,
        const string &indicator_header, const string &indicator,
        const vector<double> *performanceVector, double scalingFactor); // chip performance

// performanceType: 0 - latency, 1 - dynamicEnergy, 2 - area
void savePerformanceVector2(int performanceType, const string &filename,
        const string &indicator_header, const string &indicator,
        const vector<vector<double>> *performanceVector, 
        double scalingFactor, int numHierarchy, int bias); // support only hObjBreakdown!


void printLatencyVector(const vector<vector<double>> *latencyVector, 
                double totalLatency, double clkPeriod_ns,
                int numHierarchy, int numLayer,
                bool hObjectBreakdown, bool layerBreakdown);

void printEnergyVector(const vector<vector<double>> *energyVector, 
                double totalEnergy, 
                int numHierarchy, int numLayer,
                bool hObjectBreakdown, bool layerBreakdown);

void printAreaVector2(const vector<vector<double>> *areaVector, 
                double totalArea, int numHierarchy, 
                const vector<HierarchyObject*> *hObjectVector);

vector<vector<double>> updateAreaVector2(const vector<vector<double>> *areaVector, 
                double totalArea, int numHierarchy, 
                const vector<HierarchyObject*> *hObjectVector);

#endif /* NEUROHUB_UTIL_H_ */
