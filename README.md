# EIT Simulation Framework

Integration between EIDORS and PSPICE for EIT simulation.

 - A/D and D/A simulation;
 - Automatic generation of PWL stimulus for excitation signal and multiplexing commands;
 - Example Testbench code;
 - Instructions for PSPICE implementation.
 
Requirements: 
- EIDORS and standard MATLAB libraries;
- PSPICE or other alternative SPICE simulator with PWL functions.

Simulation Steps (from examples testbench_stimulation.m and testbench_acquisition.m):

- EIDORS SIDE (testbench_stimulation.m)
  - Set the path variables on the testbench; 
  - Set the desired FEM model, stimulation pattern and homogeneous/inhomogeneous structures using EIDORS;
  - Convert the image structures into spice netlists to the .lib files at netlist_path;
  - Configure stimulation signal using DAC_MODEL class and pwl_write() function;
  - Configure control and trigger signals using MUX_CONTROL class and pwl_write() function;
  - Run testbench_stimulation.m to generate PWL and netlist files. 
  
- PSPICE SIDE
  - Create CAPTURE parts for the .lib FEM netlists using MODEL EDITOR; 
  - Use VPWL(forever) source for the stimulation signal. Associate the part with the 'DA_output.txt' file;
  - Use VPWL sources for the switching control signal bus. Associate the parts to the MUX files, each source being a control bit for each mux;
  - Obs: 'PWL_paths.txt' can be used to facilitate the two previous steps;
  - Build your dedicated EIT electronic circuit, simulating each stage, electrode and sample (FEM netlists). Run a transient analysis (Remember to simulate until mux.time(end) seconds, choosing "maximum step size" to be at least 100 times smaller than your ADC sampling period);
  - Export the voltage measurement from the measuring circuit, using .csv files for the homogeneous and inhomogenous cases;
  
- EIDORS SIDE (testbench_acquisition.m)
  - Read PSPICE files using READ_CURVES class;
  - Configure the analog-to-digital converter using ADC_MODEL;
  - Sample and discretize PSPICE data using ADC_MODEL methods;
  - Pre-process the data using ADC_MODEL.avg(), ADC_MODEL.avg_norm() or another dedicated pre-processing method that is compatible with the implemented reconstruction algorithm;
  - Create EIDORS data structures using the pre-processed data;
  - Configure the reconstruction method using the EIDORS functions;
  - Run testbench_acquisition.m to reconstruct the ideal and experimental images (ideal images were created by testbench_stimulation.m). 
