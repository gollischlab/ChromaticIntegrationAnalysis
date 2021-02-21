# ChromaticIntegrationAnalysis

Code for analyzing chromatic integration properties of recorded retinal ganglion cells in the data, accompanying the paper "Linear and nonlinear chromatic integration in the mouse retina" by Khani and Gollisch, 2021.

The Matlab file "Analyze_Chromatic_Integration_Stimulus.m" shows how the recorded spike times of a ganglion cell under the chromatic integration stimulus are analyzed, as described in the paper.

The code reads the spike times from file, recreates the sequence of chromatic contrast steps according to a random-number generator and plots the firing rates and PSTHs for the different chromatic contrast combinations.

To run the sript, use the command "Analyze_Chromatic_Integration_Stimulus" and select one of the provided data folders in the ensuing file dialogue, e.g. "example cell 1".

For more data, see the corresponding data repository at https://gin.g-node.org/gollischlab/Khani_and_Gollisch_2021_RGC_spike_trains_chromatic_integration.

The code relies on the mex file "ran1.mexw64", which contains the random-number generator. It should run as is under Windows 10, 64 bit. If not, you may have to recompile the provided .cpp file ("mex ran1.cpp" in Matlab).
