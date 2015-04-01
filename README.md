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
Parameters  | Double*: 
            | Double*: 
            | Double*: 
            | Int: 
            | Double: Gas6 stimulation concentration in nanomolar.
            | Int: 
Returns     | Int: 





