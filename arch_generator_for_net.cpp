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

vector<int> get_network_max(vector<vector<int>> networkStructure) {
                    //int* max_fanIn, int* max_fanOut, int* max_fmap) {
    // initialize values
    int max_fanIn = 0;
    int max_fanOut = 0;
    int max_fmap = 0;
    int max_mpWindow = 0;

    // scan network
    for (int i = 0; i < networkStructure.size(); i++ ) {
        /* get layerSturcture */
        vector<int> layerStructure = networkStructure[i];

        /* (0-inW, 1-inH, 2-inC, 3-kW, 4-kH, 5-outC, 6-maxPool, 7-padding, 8-stride) */
        int inW = layerStructure[0];
        int inH = layerStructure[1];
        int inC = layerStructure[2];
        int kW = layerStructure[3];
        int kH = layerStructure[4];
        int outC = layerStructure[5];
        int mpWindow = (layerStructure[6]==1)? 4 : 0;
        int fanIn = kW * kH * inC;   // fanIn/Out - number of in/out processed in array
        int fanOut = outC * param->numColPerSynapse;
        int outW = (int)ceil( ( inW + 2.0*layerStructure[7] - kW ) / layerStructure[8] ) + 1;
        int outH = (int)ceil( ( inH + 2.0*layerStructure[7] - kH ) / layerStructure[8] ) + 1;
        int numConv = outW * outH;

        int in_fmap = inW * inH * inC;
        int out_fmap = outW * outH * outC;
        int fmap = in_fmap + out_fmap;
        
        max_fanIn = (max_fanIn < fanIn)? fanIn : max_fanIn;
        max_fanOut = (max_fanOut < fanOut)? fanOut : max_fanOut;
        max_fmap = (max_fmap < fmap)? fmap : max_fmap;
        max_mpWindow = (max_mpWindow < mpWindow)? mpWindow : max_mpWindow;
    }

    vector<int> results = {max_fanIn, max_fanOut, max_fmap};
    return results;
}

int main(int argc, char * argv[]) {
    
    /* set const info */
    const int maxNumDU = 128;

    /* get input information */
    vector<vector<int>> networkStructure;
    vector<vector<int>> archParams;

    archParams = readCSVint(argv[1]);
    networkStructure = readCSVint(argv[2]);

    int numHierarchy = archParams.size();

    // define weight/input/memory precision from wrapper
    param->synapseBit = atoi(argv[3]);      // precision of synapse weight
    param->numBitInput = atoi(argv[4]);     // precision of input neural activation
    int numBitInput = param->numBitInput;
    if (param->cellBit > param->synapseBit) {
        cout << "ERROR!: Memory precision is even higher than synpase precision, please modify 'cellBit' in Param.cpp!" << endl;
        param->cellBit = param->synapseBit;
    }
    int numColPerSynapse = (int)ceil((double)param->synapseBit/(double)param->cellBit);
    param->numColPerSynapse = numColPerSynapse;

    int scheduler_type = atoi(argv[5]);

    // define filename for saving architecture design
    string filename = argv[6];

    /* get properties of CIM array output */
    HierarchyRoot *hRoot;
    hRoot = new HierarchyRoot(inputParameter, tech, cell);
    hRoot->Initialize();
    int bitArrayOut = (int)hRoot->numOutBit;


    /* generate architecture design with given info */
    vector<vector<int>> designArch;
    int bitSubObjectOut = bitArrayOut;
    int bitObjectOut;
    int numFanInSubObject = param->numRowCIMArray;
    int numFanOutSubObject = param->numColCIMArray;
    int numFanInObject, numFanOutObject;
    int numSubObjectInExt, numSubObjectOutExt;

    /* update non-Top hierarchy design */
    for (int h=0; h < numHierarchy-1; h++) {
        vector<int> designHObj;

        /* get hObject param & imporant design properties */
        vector<int> hObjParams = archParams[h];
        int buType = hObjParams[6];
        bool hasInputBuffer = ((buType!=0) & (hObjParams.back()==0))? false : true;
        bool colSystolic = (hObjParams[2]==1/*LinearArray*/)? true : false;
        bool rowSystolic = (hObjParams[4]==1/*LinearArray*/)? true : false;
        // check the type of subObject address, 0: row/col for input/output extension, 1: linear address
        int idxType = ( hasInputBuffer || colSystolic || rowSystolic )? 0 : 1;

        /* get & update number of subObject */
        int numSubObjectRow = hObjParams[0];
        int numSubObjectCol = hObjParams[1];
        int numSubObject = numSubObjectRow * numSubObjectCol;
        designHObj.push_back(numSubObjectRow);
        designHObj.push_back(numSubObjectCol);

        /* get the number of arrays used for in/out extension */
        if ( idxType == 0 ) { // row/col for input/output extension
            numSubObjectInExt = numSubObjectRow;
            numSubObjectOutExt = numSubObjectCol;
        } else {
            cerr << "[ERROR] Non-top object should hve idxType 0!!" << endl;
            exit(-1);
        }
        numFanInObject = numFanInSubObject * numSubObjectInExt;
        numFanOutObject = numFanOutSubObject * numSubObjectOutExt;

        /* get & update digital unit info */
        // numUnit
        int numDU = MIN(maxNumDU, 
                        (int)ceil((double) numFanOutSubObject / numColPerSynapse));
        // adderTree
        designHObj.push_back(numDU); // numUnit -- need to know local output extension size
        designHObj.push_back(bitSubObjectOut); // numAdderBit
        designHObj.push_back(numSubObjectInExt); // numAdd -- need to know local input extension size
        bitObjectOut = bitSubObjectOut + (int)ceil(log2(numSubObjectInExt));
        // reLu & MaxPool (Single layer processing requires these units for last level only)
        designHObj.push_back(0); // reLu, numUnit
        designHObj.push_back(0); // reLu, numBit
        designHObj.push_back(0); // MaxPool, numUnit
        designHObj.push_back(0); // MaxPool, numBit
        designHObj.push_back(0); // MaxPool, Window

        /* update interconnect info */
        designHObj.push_back(0);    // delaytolerance if fixed to 0
        designHObj.push_back(hObjParams[2]);     // outType
        designHObj.push_back(hObjParams[3]);     // outBW / flitSize(Mesh)
        designHObj.push_back(hObjParams[4]);     // inType / numPorti(Mesh)
        designHObj.push_back(hObjParams[5]);     // InBW

        /* get & update buffer unit info */
        // bet numBit buffering data
        int inBUSize, outBUSize;
        // get #bit of output stored in the buffer
        outBUSize = ceil( param->numColCIMArray / numColPerSynapse ) * bitObjectOut
                    * numSubObjectOutExt;

        // check if the next hObj has Linear array IC for Row Systolic operation
        vector<int> nextHObjParams = archParams[h+1];
        bool nextRowSystolic = (nextHObjParams[4]==1/*LinearArray*/)? true : false;
        // get #bit of input stored in the buffer
        if ( nextRowSystolic ) {
            inBUSize = numFanInObject * numBitInput;
        } else {
            inBUSize = param->numRowCIMArray * numBitInput * numSubObjectInExt;
        }

        // get buffer info
        int outBUCoreBW, numOutBUCore, inBUCoreBW, numInBUCore;
        outBUCoreBW = hObjParams[7]; // coreBW is the size of buffer core row/col size 
        inBUCoreBW = hObjParams[8];
        // add inBUSize ato outBUSize if no seperated input buffer
        if ( !hasInputBuffer ) {
            outBUSize = outBUSize + inBUSize;
            inBUSize = 0;
        }
        // get number of BU core
        if ( buType > 0 ) {
            numOutBUCore = (int)ceil((double) outBUSize / (outBUCoreBW * outBUCoreBW) );
            if ( inBUCoreBW > 0 ) {
                numInBUCore = (int)ceil((double) inBUSize / (inBUCoreBW * inBUCoreBW) );
            } else {
                numInBUCore = 0;
            }
        } else {
            outBUCoreBW = outBUSize;
            inBUCoreBW = inBUSize;
            numOutBUCore = 1;
            numInBUCore = 1;
        }
        designHObj.push_back(buType); // buType
        designHObj.push_back(outBUSize); // outBUSize
        designHObj.push_back(outBUCoreBW); // outBUCoreBW 
        designHObj.push_back(numOutBUCore); // numOutBUCore
        designHObj.push_back(inBUSize); // inBUSize
        designHObj.push_back(inBUCoreBW); // inBUCoreBW
        designHObj.push_back(numInBUCore); // numInBUCore

        /* update designArch */
        designArch.push_back(designHObj);
       
        /* update current hObject info as next subObject info */
        numFanInSubObject = numFanInObject;
        numFanOutSubObject = numFanOutObject;
        bitSubObjectOut = bitObjectOut;

    }

    /* update the Top hierarchy design */
    vector<int> designHObj;
    /* get hObject param & imporant design properties */
    vector<int> hObjParams = archParams.back();
    int buType = hObjParams[6];

    /* get & update number of subObject */
    int numSubObjectRow = hObjParams[0];
    int numSubObjectCol = hObjParams[1];
    int numSubObject = numSubObjectRow * numSubObjectCol;

    /* Check minimum #subObject to finish the network operation on the accelerator*/
    // initialize the architecture
    vector<vector<int>> designArch_tmp;
    designArch_tmp.assign(designArch.begin(), designArch.end());
    //designArch_tmp.push_back({1, 1, 
    designArch_tmp.push_back({hObjParams[0], hObjParams[1], 
                    0,0,0, 0,0, 0,0,0, 
                    0,hObjParams[2],hObjParams[3],hObjParams[4],hObjParams[5],
                    hObjParams[6],hObjParams[7],hObjParams[7],hObjParams[7],hObjParams[8],hObjParams[8],hObjParams[8]});
    HierarchyObject *prevObject;
    HierarchyObject *lastObject;
    for ( int h=0; h < numHierarchy; h++ ) {
        vector<double> design_db;
        for ( int i=0; i < designArch_tmp[h].size(); i++ ) {
            design_db.push_back((double)designArch_tmp[h][i]);
        }
        lastObject = new HierarchyObject(inputParameter, tech, cell, h+1, hRoot, prevObject, design_db);
        lastObject->Initialize(param->clkFreq);
        prevObject = lastObject;
        design_db.clear();
    }
    designArch_tmp.clear();

    // calculate the minimum number of subObject required for top
    vector<vector<double>> networkStructure_db; // network structure with double format
    for ( int i=0; i < networkStructure.size(); i++ ) {
        vector<double> layerStructure_db;
        for ( int j=0; j < networkStructure[i].size(); j++ ) {
            layerStructure_db.push_back((double)networkStructure[i][j]);
        }
        networkStructure_db.push_back(layerStructure_db);
        layerStructure_db.clear();
    }
    NetworkScheduler *networkScheduler = new NetworkScheduler();
    networkScheduler->Initialize(networkStructure_db, lastObject); // lastObject ->hTop
    networkScheduler->Scheduling(scheduler_type);
    int numMinSubObject_top = networkScheduler->numUsedSubObject_top;
    networkStructure_db.clear();

    //printf("numMinSubObject_top: %d\n", numMinSubObject_top);

    // update the number of subObject
    numSubObject = numMinSubObject_top;
    numSubObjectRow = (int)ceil( sqrt((double) numSubObject) ); 
    numSubObjectCol = (int)ceil((double) numSubObject / numSubObjectRow );
    //if ( numSubObject < numMinSubObject_top ) {
    //    numSubObject = numMinSubObject_top;
    //    numSubObjectRow = (int)ceil( sqrt((double) numSubObject) ); 
    //    numSubObjectCol = (int)ceil((double) numSubObject / numSubObjectRow );
    //}

    designHObj.push_back(numSubObjectRow);
    designHObj.push_back(numSubObjectCol);
    //designHObj.push_back(hObjParams[0]);
    //designHObj.push_back(hObjParams[1]);

    /* get network max info */
    vector<int> network_max = get_network_max(networkStructure);
    int max_fanIn = network_max[0];
    int max_fanOut = network_max[1];
    int max_fmap = network_max[2];
    int mpWindow = network_max[3];
    network_max.clear();

    /* get the number of arrays used for in/out extension */
    numSubObjectInExt = (int)ceil( (double) max_fanIn / numFanInSubObject );
    numSubObjectOutExt = (int)floor( (double) max_fanOut / numFanOutSubObject );
    numFanInObject = numFanInSubObject * numSubObjectInExt;
    numFanOutObject = numFanOutSubObject * numSubObjectOutExt;

    /* get & update digital unit info */
    // numUnit
    int numDU = MIN(maxNumDU, 
                    (int)ceil((double) numFanOutSubObject / numColPerSynapse));
    // adderTree
    designHObj.push_back(numDU); // numUnit -- need to know local output extension size
    designHObj.push_back(bitSubObjectOut); // numAdderBit
    designHObj.push_back(numSubObjectInExt); // numAdd -- need to know local input extension size
    bitObjectOut = bitSubObjectOut + (int)ceil(log2(numSubObjectInExt));
    // reLu & MaxPool (Single layer processing requires these units for last level only)
    // reLu
    designHObj.push_back(numDU); // numUnit
    designHObj.push_back(numBitInput); // numBit
    // MaxPool
    if ( mpWindow > 0 ) {
        designHObj.push_back((int)ceil((double)numDU / mpWindow)); // numUnit
        designHObj.push_back(numBitInput); // numBit
        designHObj.push_back(mpWindow); // Window
    } else {
        designHObj.push_back(0); // numUnit
        designHObj.push_back(0); // numBit
        designHObj.push_back(0); // Window
    }

    /* update interconnect info */
    designHObj.push_back(0);    // delaytolerance if fixed to 0
    designHObj.push_back(hObjParams[2]);     // outType
    designHObj.push_back(hObjParams[3]);     // outBW / flitSize(Mesh)
    designHObj.push_back(hObjParams[4]);     // inType / numPorti(Mesh)
    designHObj.push_back(hObjParams[5]);     // InBW

    /* get & update buffer unit info */
    // bet numBit buffering data
    // top store all the input feature map & output feature map on the buffer
    int outBUSize = max_fmap;
    int inBUSize = 0;
    // get buffer info
    int outBUCoreBW, numOutBUCore, inBUCoreBW, numInBUCore;
    outBUCoreBW = hObjParams[7]; // coreBW is the size of buffer core row/col size 
    inBUCoreBW = 0;
    // get number of BU core
    if ( buType > 0 ) {
        numOutBUCore = (int)ceil((double) outBUSize / (outBUCoreBW * outBUCoreBW) );
        numInBUCore = 0;
    } else {
        outBUCoreBW = outBUSize;
        inBUCoreBW = inBUSize;
        numOutBUCore = 1;
        numInBUCore = 1;
    }
    designHObj.push_back(buType); // buType
    designHObj.push_back(outBUSize); // outBUSize
    designHObj.push_back(outBUCoreBW); // outBUCoreBW 
    designHObj.push_back(numOutBUCore); // numOutBUCore
    designHObj.push_back(inBUSize); // inBUSize
    designHObj.push_back(inBUCoreBW); // inBUCoreBW
    designHObj.push_back(numInBUCore); // numInBUCore

    /* update designArch */
    designArch.push_back(designHObj);
    

    ///* print designArch */
    //printIntVector2(&designArch);

    /* save designArch */
    saveIntVector2(filename, &designArch);

}
