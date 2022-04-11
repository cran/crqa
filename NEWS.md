# crqa (cross-recurrence quantification analysis)

# crqa 2.0.3

* Fixed minor bug on drpfromts() when method was `mdcrqa`. The argument `datatype` was forcing the input data to be in vector form, and should instead be left as a matrix. Added if statement to check (line 25)

# crqa 2.0.1

* Fixed bug on wincrqa() and windowdrp() to run windowed RQA on multidimensional data.
    * Simplified output of wincrqa() now returning a dataframe
    * Included new import from FSA() package to use diags() function. 

# crqa 2.0.0

# Major features

* Extension of recurrence analysis to multidimensional time-series data, and significant update of computational procedures in `crqa()`
    * Added arguments: `method`, `metric` and `datatype`
    * Improved method for phase-space reconstruction (line 207:237)
    * Included computation of categorical entropy (line 411:420)
    * Removal of argument `checkl`
    
* Simplified structure of functions and better division between core and ancillary functions:
    * `crqa_helpers` contains several functions previously exported (e.g., `theiler` or `tt`) that are now only accessed internally by the `crqa()` package.
    
* New functions `mdDelay()` and `mdFnn()` to estimate Average Mutual Information and False Nearest Neighbours of multidimensional time-series.

* Experimental `piecewiseRQA()` function created to better handle the computational load of large time-series.

* `optimizeParam()` now works also with multidimensional time-series. 
    
# Minor features

* Deprecated functions: `CTcrqa()`, `runcrqa()`, `calcphi()`, `takephi()`

* On `crqa()` 
    * improved error checking and warning messages (line 133:170)
    * added two more ways of rescaling the distance matrix (line 252:272)
    
* On `plotRP()` added arguments to improve the plotting of Recurrence Plots 

* `drpdfromts()` is now called `drpfromts()` and it has been rewritten to align with the new version of `crqa()`.

* A convenience function called `numerify` in `crqa_helpers` is automatically called when a user inputs categorical series (i.e., it contains either characters or factors) and this function is used to recode the levels of such time-series into numerical codes (to run crqa). A warning is send to the user when `crqa()`.

* `windowdrp()` has been rewritten to align with the new version of `crqa()`.  

* `wincrqa()` has been rewritten to align with the new version of `crqa()` and better names for the output were provided.  

# crqa 1.0.9

* On `drpdfromts()` removed a left over constant used for testing (line: 42)

# crqa 1.0.8

## Improvements and bug fixes

* On `drpdfromts()` fixed initialisation of dimensions for empty RP (line: 52)

# crqa 1.0.7

## New functions

* `plotRP()` convenience function based on the standard `plot()` to visualize a Recurrence Plot

## Improvements and bug fixes

* On `crqa()` added a few more checks (`stop`) if the data inputted did not comply with the function, and send a warning message.

* `drpdfromts()` entirely rewritten around `crqa()` to better deal with continuous valued time-series.

* `runcrqa()` fixed to fit with the revised functions: `drpdfromts()` and `windowdrp()`, 

* `tt()` fixed `rBind` (line 23 and 91) which was deprecated from the `Matrix()` package.

* `windowdrp()`entirely rewritten around the new version of `drpdfromts()` 

* `wincrqa()` adjusted indexing of windows (line: 40:41)

# crqa 1.0.6

## New functions

* `ami()` externalised from `optimizeParam()`

* `lorenzattractor()` simulates and plots 3D data from a Lorenz Attractor

## Improvements and bug fixes

* On `crqa()` include a `stop` (line: 95:100) if time-series were shorter than their phase space reconstructed portaits

* On `optimizeParam()` included argument `typeami` to set the type of `ami()` desired (either, minimum dip or maximum lag)

# crqa 1.0.5

## Improvements and bug fixes

* On `checkts()`. added argument `pad` (line: 38:67), which gives the option, in case of series of different length to extend the shortest sequence either with mean value (if the variable is in a continuous scale) or with a random label not present in either series (if the variable is categorical).

* On `crqa()`: 
    * Added default values for arguments in the call of the function. 
    * Included the argument `side` to select the region of the Recurrence Plot to extract measures on
    * Included the argument `checkl`, a wrapper to call the function `checkts()` directly inside this function. 

* On `optimizeParam()`:
    * Added arguments `min.rec` and `max.rec` to the call. 
    * Added argument within `par` (`fnnpercent`) to estimate False Nearest Neighbours based on a percentage reduction with respect to first dimension (line: 199:241)
    * Simplified code throughout and markedly improved the estimation of the radius (line: 299:360)
    
* On `runcrqa()` included argument `pad` in the call to work with the revised version of `checkts()`

* On `wincrqa()` added calculation of TREND (line: 64:87) 

# crqa 1.0.4

* crqa 1.0.3 did not pass CRAN check because of `DESCRIPTION` (`Depends` field) and was resubmitted as 1.0.4.

# crqa 1.0.3

## New functionx 

* Added `theiler()` function to choose the separation between values on the time series when specifying a delay reconstruction vector, i.e., the Theiler window.

## Improvements and bug fixes

* On `crqa()` added theiler window (`theiler()`, line: 170:175)

* On `optimizeParam()` improved calculation of average mutual information (line 56:119), and added estimation of radius within user specified expected recurrence values (line 237:284)

* On `runcrqa()` added missing arguments when calling `wincrqa()` (line 104-110:112) 

* On `tt()` simplified calculation of laminarity (line 65)

* On `wincrqa()` added missing arguments in the function call to exploit better functionality in `crqa()`

# crqa 1.0.2

* On `tt()` added in line comments as header of function.

# crqa 1.0.1

## Improvements and bug fixes

* On `crqa()`: Implementation of embedding dimensions and phase space recostruction, added in line comments in the code.

# crqa 1.0

First version of the package featuring the following original functions:

* `calcphi`
* `checkts`
* `crqa`
* `CTcrqa`
* `drpdfromts`
* `optimizeParam`
* `runcrqa`
* `simts`
* `spdiags`
* `takephi`
* `tt`
* `wincrqa`
* `windowdrp`



