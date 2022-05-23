#include <cmath>
#include <iostream>
#include <stdlib.h>
#include "LayerScheduler.h"
#include "Param.h"

extern Param *param;

LayerScheduler::LayerScheduler() {
    // set scheduling parameter
    numBitInput = param->numBitInput;
    numCellPerSynapse = param->numColPerSynapse;
    numColMuxed = param->numColMuxed;
    lengthInfoRead = param->lengthInfoRead;
}

void LayerScheduler::Initialize(int _layerIdx, vector<double> _layerStructure,
                                int _hlevelMappingUnit, vector<vector<int>> _idxOffsetVector,
                                const HierarchyObject* hTop){
   
    layerIdx = _layerIdx;
    layerStructure = _layerStructure;
    //hlevelMappingUnit = _hlevelMappingUnit;
    hlevelMappingUnit = hTop->hlevel - 1;
    idxOffsetVector = _idxOffsetVector;

    /* layerSturcture
    * (0-inW, 1-inH, 2-inC, 3-kW, 4-kH, 5-outC, 6-maxPool, 7-padding, 8-stride) */
    inW = layerStructure[0];
    inH = layerStructure[1];
    inC = layerStructure[2];
    kW = layerStructure[3];
    kH = layerStructure[4];
    outC = layerStructure[5];
    fanIn = kW * kH * inC;   
    fanOut = outC;
    outW = ceil( ( inW + 2*layerStructure[7] - kW ) / layerStructure[8] ) + 1;
    outH = ceil( ( inH + 2*layerStructure[7] - kH ) / layerStructure[8] ) + 1;
    numConv = outW * outH;

    // initialize parameters
    doingAct = false;
    doingMaxPool = false;
    doneAct = false;
    doneMaxPool = false;
    hlevelRowSystolic = -1;
    hlevelColSystolic = -1;
    hlevelTop = hTop->hlevel;
    numUsedSubObject_top = 0;

    // check Linear Array
    CheckLinearArray(hTop);

}

/* Recursively check if any HiearchyObject has Linear Array & get hlevel of HierarchyObject with Systolic array */
void LayerScheduler::CheckLinearArray(const HierarchyObject* hObject) {
    // Check if the HierarchyObject has Linear Array
    if ( hObject->interConnect->outType == 1 ) {
        hlevelRowSystolic = hObject->hlevel; 
        numHObjectRowSAExt = ceil((double)kH / hObject->numSubObjectRow); // NOTE: Assume ext can be finished in parent hObject
    }
    if ( hObject->interConnect->inType == 1 ) {
        hlevelColSystolic = hObject->hlevel; 
        numHObjectColSAExt = ceil((double)kW / hObject->numSubObjectCol); // NOTE: Assume ext can be finished in parent hObject
    }
    // move to next
    if ( hObject->hlevel > 1 ) {
        CheckLinearArray(hObject->subObject);
    }
}

/* HierarchyRoot Scheduling */
vector<double> LayerScheduler::HRootScheduling(const HierarchyRoot* hRoot, 
                    double idxRow, double idxCol, 
                    double weightMatrixRow, double weightMatrixCol) {

    // Check if we received weight matrix with proper size
    if ( hRoot->numRow < weightMatrixRow ) {
        cerr << "[Error] weightMatrixRow is larger than hRoot numRow" << endl;
        exit(-1);
    } else if ( hRoot->numCol < weightMatrixCol ) {
        cerr << "[Error] weightMatrixCol is larger than hRoot numCol" << endl;
        exit(-1);
    }

    // clear vector
    vector<double> infoRead;
    infoRead.assign(lengthInfoRead, 0);

    // generate inforRead
    infoRead[0] = (double)layerIdx;
    infoRead[1] = 0; // hlevel
    infoRead[2] = idxRow;
    infoRead[3] = idxCol;
    infoRead[4] = weightMatrixRow;
    infoRead[5] = weightMatrixCol;
    infoRead[6] = numBitInput;
    infoRead[7] = numCellPerSynapse;
    infoRead[8] = numConv; // numReadArray

    return infoRead;

}

/* HierarchyObject Scheduling 
* Scheduler type00 - base (tile-wise mapping) */
vector<vector<double>> LayerScheduler::HObjectScheduling_00(const HierarchyObject *hObject,
                                double idxRow, double idxCol, 
                                double wRow, double wCol, double wInC, double wOutC,
                                bool offsetObject, bool lastObject) {

    // Check if the kernel is split in spatial dim or Channel dim
    // Row-dim
    bool spatialRowMapping, spatialColMapping;
    if ( hlevelRowSystolic < 0 ) { // no spatial row-dim mapping for the entire mapping
        spatialRowMapping = false; 
        wInC = wRow * wInC;
        wRow = 1;
    } else if ( hObject->hlevel == hlevelRowSystolic ) {
        spatialRowMapping = true;
    } else if ( (hObject->hlevel == hlevelRowSystolic + 1) && (numHObjectRowSAExt > 1) ) {
        /* NOTE: Systolic spill is allowed for one more higher hObect */
        spatialRowMapping = true;
    } else {
        spatialRowMapping = false;
    }
    // Col-dim
    if ( hlevelColSystolic < 0 ) { // no spatial col-dim mapping for the entire mapping
        spatialColMapping = false; 
        wInC = wCol * wInC;
        wCol = 1;
    } else if ( hObject->hlevel == hlevelColSystolic ) {
        spatialColMapping = true;
    } else if ( (hObject->hlevel == hlevelColSystolic + 1) && (numHObjectColSAExt > 1) ) {
        /* NOTE: Systolic spill is allowed for one more higher hObect */
        spatialColMapping = true;
    } else {
        spatialColMapping = false;
    }

    // get idxOffset info
    int idxOffsetRow=0;
    int idxOffsetCol = 0;
    if ( offsetObject && ( idxOffsetVector.size() > 0 ) ) {
        vector<int> idxOffset = idxOffsetVector.back();        
        if ( idxOffset[0] == ( hObject->hlevel - 1 ) ) { // subObject offset
            // update offset
            idxOffsetRow = idxOffset[1];
            idxOffsetCol = idxOffset[2];
            idxOffsetVector.pop_back();
        } else {
            offsetObject = false;
            idxOffsetVector.clear();
        }
    } else { 
        offsetObject = false;
    }

    // check the type of subObject address, 0: row/col for input/output extension, 1: linear address
    int idxType = ( hObject->inputBuffer || spatialColMapping || spatialRowMapping )? 0: 1;

    // check if IC info
    const InterConnect *IC = hObject->interConnect;
    int icOutType = IC->outType;
    int icInType = IC->inType;
    BusMode icInBusMode = IC->inBusMode;

    // generate vector
    vector<vector<double>> infoReadVector;

    //*** [STEP1] Scheduling Sub-HierarchyObject ***//

    // Calculate the number of used subObject // NOTE: modify for intra-layer SubObject splitting
    double numSubObjectRow, numSubObjectCol, numSubObject;
    // Row-dim
    if ( spatialRowMapping ) {
        numSubObjectRow = wRow;
    } else {
        numSubObjectRow = ceil( wInC / hObject->numInCSubObject );
    }
    // Col-dim
    if ( spatialColMapping ) {
        numSubObjectCol = wCol;
    } else {
        numSubObjectCol = ceil( wOutC / hObject->numOutCSubObject );
    }
    numSubObject = numSubObjectRow * numSubObjectCol;

    // get & update sub-object scheduling results
    int subObjectCounter = 0;
    int idxSubObjectRow, idxSubObjectCol;
    bool offsetSubObject = offsetObject;
    bool lastSubObject = false;
    double wRowSubObject, wColSubObject, wInCSubObject, wOutCSubObject;
    // for mesh scheduling
    vector<double> infoReadICIn, infoReadICOut;
    vector<vector<double>> infoReadICInVector, infoReadICOutVector;
    double numOutBitSubObject, numInBitSubObject;
    int prevIdxSubObjectRow;
    double icState = 1; /*reset before update*/

    for (int subObjectRow=0; subObjectRow < (int)numSubObjectRow; subObjectRow++) {
        for (int subObjectCol=0; subObjectCol < (int)numSubObjectCol; subObjectCol++) {
            if ( subObjectCounter == ( numSubObject - 1) ) { // last subObject
                lastSubObject = lastObject;
            }
            // update size of weight matrix assigned for each subObject
            // Row-dim
            if ( spatialRowMapping ) {
                wRowSubObject = 1; // assume only 1 pixel is assigned for spatial dim mapping
                wInCSubObject = wInC;
            } else {
                wRowSubObject = wRow;
                wInCSubObject = MIN(hObject->numInCSubObject, wInC - subObjectRow * hObject->numInCSubObject); 
            }

            // Col-dim
            if ( spatialColMapping ) {
                wColSubObject = 1; // assume only 1 pixel is assigned for spatial dim mapping
                wOutCSubObject = wOutC;
            } else {
                wColSubObject = wCol;
                wOutCSubObject = MIN(hObject->numOutCSubObject, wOutC - subObjectCol * hObject->numOutCSubObject); 
            }

            // get index 
            if ( idxType == 0 ) {
                idxSubObjectRow = idxOffsetRow + subObjectRow;
                idxSubObjectCol = idxOffsetCol + subObjectCol;
            } else { // idxType == 1
                div_t divResult = div((subObjectCounter + idxOffsetRow), hObject->numSubObjectCol);
                idxSubObjectRow = divResult.quot;
                idxSubObjectCol = divResult.rem;
            }
            // check if the idx is valid
            if ( (idxSubObjectRow < 0) || (idxSubObjectRow >= hObject->numSubObjectRow) ||
                 (idxSubObjectCol < 0) || (idxSubObjectCol >= hObject->numSubObjectCol) ) {
                //cout << hlevelTop << " " << hObject->hlevel << " " << idxSubObjectRow << " " << idxSubObjectCol << endl;
                //cout << hObject->numSubObjectRow << " " << hObject->numSubObjectCol << endl;
                cerr << "[Error] Invalid idxSubObjectRow/Col!" << endl;
                exit(-1);
            }
        
            // scheduling sub-object 
            if ( hObject->hlevel == 1 ) {
                if ( (wRowSubObject > 1) || (wColSubObject > 1) ) {
                    cerr << "[Error] wRow/Col for hRootObject should be 1!" << endl;
                    exit(-1);
                }
                vector<double> infoReadSub = HRootScheduling( hObject->rootObject, idxSubObjectRow, idxSubObjectCol, 
                                                            wInCSubObject, wOutCSubObject);
                infoReadVector.push_back(infoReadSub); // update scheduling result
                infoReadSub.clear(); // clear updated result
            } else {
                vector<vector<double>> infoReadVectorSub 
                                            = HObjectScheduling_00( hObject->subObject, idxSubObjectRow, idxSubObjectCol, 
                                                            wRowSubObject, wColSubObject, wInCSubObject, wOutCSubObject,
                                                            offsetSubObject, lastSubObject);
                infoReadVector.insert(infoReadVector.end(), infoReadVectorSub.begin(), infoReadVectorSub.end()); // update
                infoReadVectorSub.clear(); // clear updated result
            }

            // update infoReadIC
            // (layerIdx, hlevel, idxRow, idxCol, state, dataType, numRead, numHopsRow, numHopsCol, dummy)
            if ((icOutType == 2 /*2D Mesh*/) || (icOutType == 3 /*HBus*/)) {
                //printf("wOutC: %10d, wOutCSubObject: %10d\n", (int)wOutC, (int)wOutCSubObject);
                numOutBitSubObject = numConv * wOutCSubObject / numCellPerSynapse * hObject->numOutBitSubObject;
                infoReadICOut = {(double)layerIdx, (double)hObject->hlevel, idxRow, idxCol, icState, 0/*dataType*/,
                                numOutBitSubObject / IC->outBusWidth, (double)idxSubObjectRow, 0, 0}; // no col dim move of out
                infoReadICOutVector.push_back(infoReadICOut);
                icState = 2; /*cumulate*/
            }

            if ((icInType == 2 /*2D Mesh*/) || (icInType == 3 /*HBus*/)) {
                if ( subObjectCol == 0 ) { // input switch -> get feeding info of new input
                    if ( hlevelColSystolic > 0 ) {
                        if ( hObject->hlevel >= hlevelColSystolic ) {
                            numInBitSubObject = outH * ( wRowSubObject * inW * wInCSubObject) * numBitInput;
                        } else {
                            numInBitSubObject = outH * ( wRowSubObject * outW * wInCSubObject) * numBitInput;
                        }
                    } else {
                        numInBitSubObject = numConv * ( wRowSubObject * wColSubObject * wInCSubObject) * numBitInput;
                    }
                    infoReadICIn = {(double)layerIdx, (double)hObject->hlevel, idxRow, idxCol, 2/*icState*/, 1/*dataType*/,
                                numInBitSubObject / IC->inBusWidth, (double)idxSubObjectRow, (double)idxSubObjectCol};
                } else { // input sharing
                    if ( prevIdxSubObjectRow == idxSubObjectRow ) {
                        infoReadICIn.back() = MAX(infoReadICIn.back(), idxSubObjectCol); // horizontal & same row -> update the col idx
                    } else { // rowIdx of input changed -> update prev input info & get new row input feeding info
                        if ( (icInType == 3 /*HBus*/) && (prevIdxSubObjectRow % 2 == 0)  ) { // row bus sharing for HBus
                            infoReadICIn.rbegin()[1] = (double)idxSubObjectRow;
                            infoReadICIn.back() = MAX(infoReadICIn.back(), idxSubObjectCol);
                        }
                        else {
                            infoReadICIn.push_back(0/*dummy*/);
                            infoReadICInVector.push_back(infoReadICIn); // NOTE: BUG: Latency should be overlapped
                            double rowMove = (double)(idxSubObjectRow - prevIdxSubObjectRow);
                            if ( icInType == 2 /*2D Mesh*/ ) {
                                infoReadICIn = {(double)layerIdx, (double)hObject->hlevel, idxRow, idxCol, 2, 1/*dataType*/,
                                        numInBitSubObject / IC->inBusWidth, rowMove, (double)idxSubObjectCol};
                            } else {
                                infoReadICIn = {(double)layerIdx, (double)hObject->hlevel, idxRow, idxCol, 2, 1/*dataType*/,
                                        numInBitSubObject / IC->inBusWidth, rowMove + 1, (double)idxSubObjectCol};
                            }
                        }
                    }
                }
            }

            // move to next subObject
            prevIdxSubObjectRow = idxSubObjectRow;
            subObjectCounter += 1;
            offsetSubObject = false;
        }

        // update infoReadICIn
        if ((icInType == 2 /*2D Mesh*/) || (icInType == 3 /*HBus*/)) {
            infoReadICIn.push_back(0/*dummy*/);
            infoReadICInVector.push_back(infoReadICIn);
        }

    }

    if ( hObject->hlevel == hlevelTop ) {
        numUsedSubObject_top = subObjectCounter + idxOffsetRow;
    }


    // check if there was any overflow in the subObject index
    if ( idxType == 0 ) {
        if ( ( (idxSubObjectRow + 1) > hObject->numSubObjectRow ) || ( (idxSubObjectCol + 1) > hObject->numSubObjectCol) ) {
            cerr << "[ERROR] index of subObject row overflow" << endl;
            cerr << "last idxSubObjectRow: " << idxSubObjectRow <<  ", subObjectRow: " << hObject->numSubObjectRow << endl;
            cerr << "last idxSubObjectCol: " << idxSubObjectCol <<  ", subObjectCol: " << hObject->numSubObjectCol << endl;
            exit(-1);
        }
    } else { // idxType == 1
        int numUsedSub = subObjectCounter + idxOffsetRow;
        if ( numUsedSub > hObject->numSubObject ) {
            cerr << "[ERROR] number of subObject overflow" << endl;
            cerr << "used numSubObject: " << numUsedSub << ", available numSubObject: " << hObject->numSubObject << endl;
            exit(-1);
        }
    }

    // update nextIdxOffsetVector
    if ( ( hObject->hlevel > hlevelMappingUnit ) & ( lastSubObject ) ) {
        if ( hObject->inputBuffer ) { // fixed flow
            // NOTE: only row dim split is allowed - Not used now
            if ( ( idxSubObjectRow + 1 ) < ( hObject->numSubObjectRow - 1 ) ) { 
                vector<int> nextIdxOffset = {hObject->hlevel-1, idxSubObjectRow+1, 0};
                nextIdxOffsetVector.push_back(nextIdxOffset);
                nextIdxOffset.clear();
            }
        } else { // non-fixed flow
            idxSubObjectRow = subObjectCounter + idxOffsetRow; // linear index
            idxSubObjectCol = 0;
            vector<int> nextIdxOffset = {hObject->hlevel-1, idxSubObjectRow, idxSubObjectCol};
            nextIdxOffsetVector.push_back(nextIdxOffset);
            nextIdxOffset.clear();
        }
    }

    if ( doingAct ) {
        doneAct = true;
    }
    if ( doingMaxPool ) {
        doneMaxPool = true;
    }

    //*** [STEP1] Scheduling Current HierarchyObject ***//
    vector<double> infoRead; // storage for update scheduling results
    
    // update object info
    infoRead.push_back((double)layerIdx);
    infoRead.push_back((double)hObject->hlevel); // hlevel
    infoRead.push_back(idxRow); // object row index
    infoRead.push_back(idxCol); // object col index
    infoRead.push_back(wRow * wCol * wInC); 
    infoRead.push_back(wOutC); 

    // update infoReadCIM
    infoRead.push_back(0);  // numBitInput
    infoRead.push_back(0);  // numCellPerSynapse
    infoRead.push_back(0);  // numRead Array

    //*** [STEP2] Scheduling DE ***//
    const DigitalElements *DE = hObject->digitalElements;
    double numRead, numUnitAdd;

    // number of total output 
    double numOut = numConv * wOutC / numCellPerSynapse; // NOTE: What for systolic array?  -- seems like the same for adderTree, but different in subObject addition & the size of adderTree in the object

    // adderTree
    if ( DE->placeAdderTree && !doneAct ) {
        numRead = numOut / DE->adderTree->numAdderTree;
        numUnitAdd = numSubObjectRow;
    } else {
        numRead = 0;
        numUnitAdd = 0;
    }
    infoRead.push_back(numRead);
    infoRead.push_back(numUnitAdd);

    // reLu
    if ( DE->placeReLu && (wRow * wCol * wInC == fanIn) && !doneAct) { // has reLu unit & finished MAC operation
        numRead = numOut / DE->reLu->numUnit;
        doingAct = true;
    } else {
        numRead = 0;
    }
    infoRead.push_back(numRead);

    // maxPool
    if ( DE->placeMaxPooling && (doingAct || doneAct) && !doneMaxPool) { // has maxPooling unit & finished MAC operation
        numRead = numOut / DE->maxPooling->window / DE->maxPooling->numMaxPooling;
        doingMaxPool = true;
    } else {
        numRead = 0;
    }
    infoRead.push_back(numRead);
    

    //*** [STEP3] Scheduling IC ***//
    double totalOutBit, totalInBit;
    if ( (hlevelRowSystolic > 0 ) & (hObject->hlevel == hlevelRowSystolic) ) {
        totalOutBit = numOut * hObject->numOutBitSubObject;
    } else {
        totalOutBit = numSubObjectRow * numOut * hObject->numOutBitSubObject;
    }

    if ( hlevelColSystolic > 0 ){ // NOTE:BUG: with stride, add additional FF access count for just passing!!
        if ( hObject->hlevel >= hlevelColSystolic ) {
            totalInBit = outH * (wRow * inW * wInC) * numBitInput;
        } else {
            totalInBit = outH * (wRow * outW * wInC) * numBitInput;
        }
    } else {
        totalInBit = numConv * (wRow * wCol * wInC) * numBitInput;
    }

    
    // update infoReadICVector for output
    if ((icOutType == 2 /*2D Mesh*/) || (icOutType == 3 /*HBus*/)) {
        infoReadVector.insert(infoReadVector.end(), infoReadICOutVector.begin(), infoReadICOutVector.end());
        infoReadICOutVector.clear();
    } else {
        // numOutRead
        double numOutRead = totalOutBit / IC->outBusWidth;
        
        // update IC infoRead // NOTE:BUG: assume single layer on hObject w/ Linear Array
        vector<double> tmpInfoReadIC; 
        if ( IC->outType == 1 /*LinearArray*/) {
            tmpInfoReadIC = {(double)layerIdx, (double)hObject->hlevel, idxRow, idxCol, 1/*icState*/, 0/*dataType*/,
                            numOutRead, (double)numSubObjectRow, (double)numSubObjectCol, 0/*dummy*/};
        } else{
            tmpInfoReadIC = {(double)layerIdx, (double)hObject->hlevel, idxRow, idxCol, 1/*icState*/, 0/*dataType*/,
                            numOutRead, 0, 0, 0/*dummy*/};
        }
        infoReadVector.push_back(tmpInfoReadIC);
        tmpInfoReadIC.clear();
    }

    // update infoReadICVector for input
    if ((icInType == 2 /*2D Mesh*/) || (icInType == 3 /*HBus*/)) {
        infoReadVector.insert(infoReadVector.end(), infoReadICInVector.begin(), infoReadICInVector.end());
        infoReadICInVector.clear();
    } else {
        double numInRead;

        // numInRead
        if ( IC->inType == 0 /*Bus*/) {
            numInRead = totalInBit / IC->inBusWidth;
        } else { /*Linear Array*/
            numInRead = ( totalInBit / wRow ) / IC->inBusWidth;
        }

        // update IC infoRead // NOTE:BUG: assume single layer on hObject w/ Linear Array
        vector<double> tmpInfoReadIC; 
        if ( IC->inType == 1 /*LinearArray*/) {
            tmpInfoReadIC = {(double)layerIdx, (double)hObject->hlevel, idxRow, idxCol, 2/*icState*/, 1/*dataType*/,
                            numInRead, (double)numSubObjectRow, (double)numSubObjectCol, 0/*dummy*/};
        } else{
            tmpInfoReadIC = {(double)layerIdx, (double)hObject->hlevel, idxRow, idxCol, 2/*icState*/, 1/*dataType*/,
                            numInRead, 0, 0, 0/*dummy*/};
        }
        infoReadVector.push_back(tmpInfoReadIC);
        tmpInfoReadIC.clear();
    }

    // empty IC infoRead
    vector<double> tmpInfoReadIC = {0, 0, 0, 0, 0, 0};
    infoRead.insert(infoRead.end(), tmpInfoReadIC.begin(), tmpInfoReadIC.end());
    tmpInfoReadIC.clear();


    //*** [STEP4] Scheduling BU ***//
    const BufferUnit *BU = hObject->bufferUnit;
    double numOutRead, numOutWrite, numInRead, numInWrite, outParallelism, inParallelism;

    // get output data size stored in buffer unit (BU)
    double numOutBitBU, numOutBU;
    if ( doneMaxPool || doingMaxPool ) {
        numOutBU = numOut / DE->maxPooling->window;
        numOutBitBU = numBitInput;
    } else {
        numOutBU = numOut;
        if ( doneAct || doingAct ) {
            numOutBitBU = numBitInput;
        } else {
            if ( DE->placeAdderTree ) {
                numOutBitBU = DE->adderTree->numOutBit;
            } else {
                numOutBitBU = hObject->numOutBitSubObject;
            }
        }
    }
    double totalOutBitBU = numOutBU * numOutBitBU;

    // get number of read 
    if ( BU->inBUSize > 0 ) {
        numOutRead = totalOutBitBU / BU->outBUCoreBW;
        numInRead = totalInBit / BU->inBUCoreBW;
    } else {
        numOutRead = (totalOutBitBU + totalInBit) / BU->outBUCoreBW;
        numInRead = 0;
    }

    // get number of write & parallelism
    if ( BU->buType == 0 /*D Flip-Flop*/ ) {
        numOutWrite = 0;
        numInWrite = 0;
        outParallelism = 1;
        inParallelism = 1;
    } else {
        numOutWrite = numOutRead;
        numInWrite = numInRead;
        // NOTE: separate IC if separate buffer unit for in/out
        if ( BU->inBUSize > 0 ) { // have input buffer
            outParallelism = MIN(BU->numOutBUCore, IC->outBusWidth/BU->outBUCoreBW);
            inParallelism = MIN(BU->numInBUCore, IC->inBusWidth/BU->inBUCoreBW);
        } else { // only output buffer
            outParallelism = MIN(BU->numOutBUCore, (IC->outBusWidth + IC->inBusWidth)/BU->outBUCoreBW);
            inParallelism = 0;
        }
    }
    infoRead.push_back(numOutRead);
    infoRead.push_back(numOutWrite);
    infoRead.push_back(outParallelism);
    infoRead.push_back(numInRead);
    infoRead.push_back(numInWrite);
    infoRead.push_back(inParallelism);
    
    //*** [STEP5] Update scheduling results ***//
    infoReadVector.push_back(infoRead);

    return infoReadVector;
}
