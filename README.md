AXLdiffEQ
=========

Gas6 signaling model for AXL


External Library

The external MatLab library may be accessed by the following function interfaces. Each can be found in the attached source code.

## pyEntry

Description | Fits parameter set to data, returning chi-squared of data fit. Single threaded.
----------- | --------------------------------------------------------------------------
Parameters  | Double*: Output vector the length of one parameter set.
Returns     | Chi squared error. Integration failures will return error of 10^6.


## pyEntryVec

Description | Fits parameter set to data, returning chi-squared of data fit. Single threaded but vectorized.
----------- | --------------------------------------------------------------------------
Parameters  | Double*: Output vector the length of one parameter set.
            | Double*: 
            | Int: Number of parameter sets passed.
Returns     | Void


## calcProfileMatlab

Description | Fits parameter set to data, returning chi-squared of data fit. Single threaded but vectorized.
----------- | --------------------------------------------------------------------------
Parameters  | Double*: Output data array, length equal to the number of time points.
            | Double*: Array of a parameter set.
            | Double*: Array of time points to calculate. Should monotonically increase.
            | Int: Number of time points (length of previous array).
            | Double: Gas6 stimulation concentration in nanomolar.
            | Int: Specifier for pY, total, or surface receptor.
Returns     | Int: 0 for success. 1 for failure.





