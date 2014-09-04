AXLdiffEQ
=========

Gas6 signaling model for TAM receptors


External Library

The external MatLab library may be accessed by the following function interfaces. Each can be found in the attached source code.

## matlabEntry

Description | Fits parameter set to data, returning chi-squared of data fit. Multithreaded.
----------- | --------------------------------------------------------------------------
Parameters  | Double*: Output vector as long as the number of parameter sets passed.
            | Double*: 
            | Int: Number of parameter sets passed in.
Returns     | Zero always. Integration failures will return error of 10^6.




\begin{tabular}{r| p{4.7in}}
\hline
\multicolumn{2}{ l }{  \texttt{rEntry} } \\
\hline
Description  & Fits parameter set to data, returning $\chi^2$ of data fit.\\
Parameters & Double*: Output value. \\
 & Const double*: Parameter set. \\
Returns & Zero always. Integration failures will return error of $10^6$. \\
\end{tabular}
\\
\begin{tabular}{r| p{4.7in}}
\hline
\multicolumn{2}{ l }{  \texttt{calcProfileMatlab} } \\
\hline
Description  & Integrates the homogeneous model and returns resulting profile. \\
Parameters & Double*: Output vector the length of the number of time points. \\
 & Double*: Vector containing the parameters. \\
 & Double*: Vector containing the time points to return (must be in order). \\
 & Int: Integer containing the number of time points. \\
 & Double: Autocrine abundance of Gas6, in nM. \\
 & Double: AXL expression in units of receptors per minute. \\
 & Double: Concentration of Gas6 stimulation, in nM.  \\
 & Int: Switch whether output should be (1) pY, (2) fractional pY, or (3) total AXL. \\
Returns & Int: Zero on success, one on failure.\\
\end{tabular}
\\
\begin{tabular}{r| p{4.7in}}
\hline
\multicolumn{2}{ l }{  \texttt{matlabDiffTPS} } \\
\hline
Description  & Integrates the spatial model and returns species abundance at specified times. \\
Parameters & Double*: Output vector length of which is the product of time point number, species, and radial discretization points. \\
 & Double: AXL expression in units of receptors per minute. \\
 & Double*: Input vector of Gas6 concentration at each spatial point in nM, length equal to the number of spatial points. \\
 & Double: Number of spatial grid points. \\
 & Double: Autocrine abundance of Gas6, in nM. \\
 & Double*: Vector containing the parameters. \\
 & Double*: Vector containing the time points to return (must be in order). \\
 & Double*: Diffusivity of each species, in units of $L^2~\textrm{min}^{-1}$. \\
 & Double: Fractional endocytosis activity of Gas6-bound species. Should usually be 1. \\
 & Double: Fractional impaired degradation of Gas6-bound species. Should usually be 1. \\
Returns & Int: Zero on success, one on failure.\\
\end{tabular}
\\
\begin{tabular}{r| p{4.7in}}
\hline
\multicolumn{2}{ l }{  \texttt{matlabDiffTPS\_pY} } \\
\hline
Description  & Integrates the spatial model and returns species abundance at specified times. \\
Parameters & Double*: Output vector the length of the number of time points times the number of spatial points. \\
 & Double: AXL expression in units of receptors per minute. \\
 & Double*: Input vector of Gas6 concentration at each spatial point in nM, length equal to the number of spatial points. \\
 & Double: Number of spatial grid points. \\
 & Double: Autocrine abundance of Gas6, in nM. \\
 & Double*: Vector containing the parameters. \\
 & Double*: Vector containing the time points to return (must be in order). \\
 & Double*: Diffusivity of each species, in units of $L^2~\textrm{min}^{-1}$. \\
 & Double: Fractional endocytosis activity of Gas6-bound species. Should usually be 1. \\
 & Double: Fractional impaired degradation of Gas6-bound species. Should usually be 1. \\
 & Int: Switch whether output should be (1) pY, (2) fractional pY, or (3) total AXL. \\
Returns & Int: Zero on success, one on failure.\\
\end{tabular}
\\
\begin{tabular}{r| p{4.7in}}
\hline
\multicolumn{2}{ l }{  \texttt{matlabDiffTPS\_pYavg} } \\
\hline
Description  & Integrates the spatial model and returns the average pY output at specified time points. \\
Parameters & Double*: Output vector the length of the number of time points. \\
 & Double: AXL expression in units of receptors per minute. \\
 & Double*: Input vector of Gas6 concentration at each spatial point in nM, length equal to the number of spatial points. \\
 & Double: Number of spatial grid points. \\
 & Double: Autocrine abundance of Gas6, in nM. \\
 & Double*: Vector containing the parameters. \\
 & Double*: Vector containing the time points to return (must be in order). \\
 & Double*: Diffusivity of each species, in units of $L^2~\textrm{min}^{-1}$. \\
 & Double: Fractional endocytosis activity of Gas6-bound species. Should usually be 1. \\
 & Double: Fractional impaired degradation of Gas6-bound species. Should usually be 1. \\
 & Int: Switch whether output should be (1) pY, (2) fractional pY, or (3) total AXL. \\
Returns & Int: Zero on success, one on failure.\\

\end{tabular}
