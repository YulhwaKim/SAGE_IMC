/*********************************************************************************************
* Definition of the Network Scheduler
*********************************************************************************************/

#ifndef NETWORKCHEDULER_H_
#define NETWORKCHEDULER_H_

#include <vector>
#include "HierarchyObject.h"

using namespace std;

class NetworkScheduler {
//private:
public:
    NetworkScheduler();
    virtual ~NetworkScheduler() {}

    /* Functions */
    void Initialize(const vector<vector<double>> _networkStructure,
                    const HierarchyObject *_hTop);

    vector<vector<double>> Scheduling(int scheduler_type); // calculate #top level objects

    void CalculatePerformance(vector<vector<double>> *networkLatencyVector,
                              vector<vector<double>> *networkEnergyVector,
                              vector<vector<double>> *layerLatencyVector2,
                              vector<vector<double>> *layerEnergyVector2,
                              vector<vector<double>> *networkLatencyVector2,
                              vector<vector<double>> *networkEnergyVector2,
                              HierarchyRoot *hRoot, vector<HierarchyObject*> hObjectVector,
                              const vector<vector<double>> networkInfoRead);

    /* Network Properties */
    vector<vector<double>> networkStructure;
    int numLayer;
    const HierarchyObject* hTop;
    int hlevelMappingUnit;
    int numUsedSubObject_top;

    /* Scheduling Result*/
    vector<vector<double>> networkInfoRead;

};

#endif /* NETWORKSCHEDULER_H_ */
