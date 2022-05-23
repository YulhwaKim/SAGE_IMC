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
    
    /* set const info */
    const int maxNumDU = 128;

    /* get input information */
    vector<int> layerStructure;
    vector<vector<int>> archParams;

    archParams = readCSVint(argv[1]);
    layerStructure = readCSVint(argv[2])[0];

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

    // define filename for saving architecture design
    string filename = argv[5];

    /* get layerSturcture
    * (0-inW, 1-inH, 2-inC, 3-kW, 4-kH, 5-outC, 6-maxPool, 7-padding, 8-stride) */
    int inW = layerStructure[0];
    int inH = layerStructure[1];
    int inC = layerStructure[2];
    int kW = layerStructure[3];
    int kH = layerStructure[4];
    int outC = layerStructure[5];
    int mpWindow = (layerStructure[6]==1)? 4 : 0;
    int fanIn = kW * kH * inC;   // fanIn/Out - number of in/out processed in array
    int fanOut = outC * numColPerSynapse;
    int outW = (int)ceil( ( inW + 2.0*layerStructure[7] - kW ) / layerStructure[8] ) + 1;
    int outH = (int)ceil( ( inH + 2.0*layerStructure[7] - kH ) / layerStructure[8] ) + 1;
    int numConv = outW * outH;

    /* calculate the number of arrays required for
    * input/output extension and total number of arrays for the extension */
    int numInputExt = (int)ceil( (double)fanIn / param->numRowCIMArray );
    int numOutputExt = (int)ceil( (double)fanOut / param->numColCIMArray);
    int totalExt = numInputExt * numOutputExt;

    /* get properties of CIM array output */
    HierarchyRoot *hRoot;
    hRoot = new HierarchyRoot(inputParameter, tech, cell);
    hRoot->Initialize();
    int bitArrayOut = (int)hRoot->numOutBit;


    /* generate architecture design with given info */
    vector<vector<int>> designArch;
    int bitSubObjectOut = bitArrayOut;
    int bitObjectOut;
    int numArraysSubObject = 1;
    int numArraysObject;
    int numFanInSubObject = param->numRowCIMArray;
    int numFanOutSubObject = param->numColCIMArray;
    int numFanInObject, numFanOutObject;
    int numSubObjectInExt, numSubObjectOutExt;
    for (int h=0; h < numHierarchy; h++) {
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
        numArraysObject = numArraysSubObject * numSubObject;
        if ( h == numHierarchy - 1 ) { // last level->meet the minimum #array requirement
            if ( numArraysObject < totalExt ) {
                numSubObject = (int)ceil((double) totalExt / numArraysSubObject );
                numSubObjectRow = (int)ceil( sqrt((double) numSubObject) ); 
                numSubObjectCol = (int)ceil((double) numSubObject / numSubObjectRow );
                numArraysObject = numArraysSubObject * numSubObject;
            }
        }
        designHObj.push_back(numSubObjectRow);
        designHObj.push_back(numSubObjectCol);

        /* get the number of arrays used for in/out extension */
        if ( idxType == 0 ) { // row/col for input/output extension
            numSubObjectInExt = numSubObjectRow;
            numSubObjectOutExt = numSubObjectCol;
        } else {
            if ( h == numHierarchy - 1 ) {
                numSubObjectInExt = (int)ceil( (double) fanIn / numFanInSubObject );
                numSubObjectOutExt = (int)ceil( (double) fanOut / numFanOutSubObject );
            } else {
                if ( numFanInSubObject * numSubObjectRow <= fanIn ) {
                    numSubObjectInExt = numSubObjectRow;
                    numSubObjectOutExt = numSubObjectCol;
                } else {
                    numSubObjectInExt = (int)ceil( (double) fanIn / numFanInSubObject );
                    numSubObjectOutExt = (int)floor( (double) numSubObject / numSubObjectInExt );
                }
            }
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
        if ( h == numHierarchy - 1 ) {
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
        } else {
            designHObj.push_back(0); // reLu, numUnit
            designHObj.push_back(0); // reLu, numBit
            designHObj.push_back(0); // MaxPool, numUnit
            designHObj.push_back(0); // MaxPool, numBit
            designHObj.push_back(0); // MaxPool, Window
        }

        /* update interconnect info */
        designHObj.push_back(0);    // delaytolerance if fixed to 0
        designHObj.push_back(hObjParams[2]);     // outType
        designHObj.push_back(hObjParams[3]);     // outBW / flitSize(Mesh)
        designHObj.push_back(hObjParams[4]);     // inType / numPorti(Mesh)
        designHObj.push_back(hObjParams[5]);     // InBW

        /* get & update buffer unit info */
        // bet numBit buffering data
        int inBUSize, outBUSize;
        if ( h == numHierarchy - 1 ) { // top hiearchy
            // top store all the input feature map & output feature map on the buffer
            outBUSize = outW * outH * outC * numBitInput;
            inBUSize = inW * inH * inC * numBitInput;
        } else {
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
        numArraysSubObject = numArraysObject;
        numFanInSubObject = numFanInObject;
        numFanOutSubObject = numFanOutObject;
        bitSubObjectOut = bitObjectOut;

    }

    /* print designArch */
    //printIntVector2(&designArch);

    /* save designArch */
    saveIntVector2(filename, &designArch);

}
