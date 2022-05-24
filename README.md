# SAGE_IMC
An open-source framework for the IMC architecture design.

SAGE-IMC consists of 
1. Architecture Template (./ArchitectureTemplate)
2. DNN Mapper (./Scheduler)
3. Performance Analysis Engine (./CircuitModule)

Note taht the performance analysis engine (./CircuitModule) is based on [NeuroSim](https://github.com/neurosim/MLP_NeuroSim_V3.0).

## Requirements
Conda environment required for the experiment is updated in 'environment.yml'.

## Makefile & Program Usage
First, compile the C++ source code.

    cmake CMakeLikst.txt
    make
    
Second, run python script to evaluate previous IMC accelerator design.

    pythons script_evaluate_previous_arch.py
    
Third, run python script to explore design space of IMC accelerators.

    pythons script_dse.py
