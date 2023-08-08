# carbondate-cpp

This depends on the R standalone library. See [here](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#The-standalone-Rmath-library) for more information.
The library source code has been included here and is built and compiled into the executable as a static library.

It expects input files as OxCal code, and looks for the command `NP_model(“name”)`, where the name for the model is required. 
Dates can then be entered in as radiocarbon age BP using the `R_Date()` command e.g.
```
Plot()
  {
  NP_Model("my_model_1")
    {
    R_Date("my_point_1",1421,32);
    R_Date("my_point_2",1551,33);
  };
};
```

The dates can also be entered as F14C concentrations using `R_F14C()`. The `R_Date()` and `R_F14C()` commands are identical 
to those used in other models in OxCal, where the first argument is optional, and assigns a label to the data point, and the last two arguments are the radiocarbon determination and the error respectively.

There are various options that are recognised by the program (note this does not include all options that OxCal recognises).

| OPTION        | DEFAULT     | DESCRIPTION                                                                                   |
|---------------|-------------|-----------------------------------------------------------------------------------------------|
| `Curve`       | `intcal20.14c` | The calibration curve to use. Only the IntCal curves are currently supported                  |
| `Floruit`     | `FALSE`	      | Whether quantile ranges are calculated instead of the default highest posterior density (hpd) |
| `kIterations` | `100`	        | The default number of MCMC passes                                                             | 
| `Resolution`  | `5`	          | The default bin size for probability distributions                                            |
| `SD1`	        | `TRUE`	       | Whether 68.2% (1 σ) ranges are given in the log and tab delimited files                       |
| `SD2`	        | `TRUE`	       | Whether 95.4% (2 σ) ranges are given in the log and tab delimited files                       |
| `SD3`	        | `FALSE`	      | Whether 99.7% (3 σ) ranges are given in the log and tab delimited files                       |
| `UseF14C`	    | `TRUE`	       | Whether all calibrations take place in F14C space (rather than BP space)                      |


This program is designed to be run after the OxCal program, therefore it expects the OxCal output `*.js`, `*.txt` and `*.log` files to be present.
When compiled with the `OXCAL_RELEASE` preprocessor definition present, it expects all input and output files to be in the working directory of the executable.
When compiled without it, it expects the directory structure in the code here.
