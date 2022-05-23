#ifndef MESH_H_
#define MESH_H_

#include "typedef.h"
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"
#include "FunctionUnit.h"
#include "DFF.h"
#include "DigitalMux.h"

class Mesh: public FunctionUnit {
public:
	Mesh(const InputParameter& _inputParameter, const Technology& _tech);
	virtual ~Mesh() {}
	const InputParameter& inputParameter;
	const Technology& tech;

	/* Functions */
	void PrintProperty(const char* str);
    void Initialize(int _numPort, int _flitSize, int _numRow, int _numCol, double _delaytolerance,
                    double _unitHeight, double _unitWidth, double _foldedratio, double _clkFreq);
	void CalculateArea();
	void CalculateLatency(int numHopsRow, int numHopsCol, double numRead);
	void CalculatePower(int numHopsRow, int numHopsCol, double numBitAccess, double numRead);

	/* Properties */
	bool initialized;	/* Initialization flag */
    int numPort, flitSize;
    int numRow, numCol, numDff, numMux, numRepeaterH, numRepeaterV;
    int numRouterRow, numRouterCol, numRouter;
    double delaytolerance, wireWidthV, wireWidthH, wireLengthV, wireLengthH;
    double unitHeight, unitWidth;
    double foldedratio;
	double clkFreq;
    double switchingRatio;

    double unitLengthWireResistance, unitLengthWireCap, unitLatencyRep, unitLatencyWire;
    double unitLengthEnergyRep, unitLengthEnergyWire;
    double repeaterSize, minDist, capMinInvInput, capMinInvOutput, hMinInv, wMinInv, widthMinInvN, widthMinInvP;
    double hInv, wInv, capInvInput, capInvOutput, widthInvN, widthInvP;
    double hRep, wRep, capRepInput, capRepOutput;
   
    DigitalMux mux;
    DFF dff;
};

#endif /* MESH_H_ */
