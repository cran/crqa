# Unidimensional and Multidimensional Methods for Recurrence Quantification Analysis with crqa.

The crqa, R package, allows users to conduct a wide range of recurrence-based analyses on single (e.g., auto-recurrence) and multivariate time series (e.g., multidimensional cross-recurrence quantification), examine coupling properties underlying leader-follower relationships (i.e., diagonal-profile methods), as well as, track the evolution of recurrence rate over the time course (i.e., windowed methods). 

## Installation

``` r
# You can install the latest version of crqa on CRAN by running:
install.packages("crqa")

# Or for the development version from GitHub:
# install.packages("devtools")
devtools::install_github("morenococo/crqa")
```

# Usage

crqa comes with some data that can be used to test and study the different functions therein. 

``` r
data(crqa) # load the data 
```

## RQA on a categorical time-series (auto-recurrence)

First, specify the arguments that will be used in the crqa, core function.

``` r
## parameter setting 
delay = 1; embed = 1; rescale = 0; radius = 0.0001;
normalize = 0; mindiagline = 2; minvertline = 2;
tw = 1; whiteline = FALSE; recpt = FALSE; 
side = "both"; method = 'rqa'; metric = 'euclidean';  
datatype = "categorical"
```

Then, run crqa on a nursery rhyme “The wheels on the bus” by Verna Hills: a vector of 120 strings (i.e., the words of the song),

``` r
ans = crqa(text, text, delay, embed, rescale, radius, normalize, 
           mindiagline, minvertline, tw, whiteline, recpt, side, method, metric, 
           datatype)
```

Have a look at the output, which contains different measures extracted from the recurrence plot (RP), and the RP itself, which can be plotted using the plotRP function. 

``` r
str(ans)
```
## CRQA on a categorical time-series (cross-recurrence)

Cross-recurrence extends univariate recurrence analysis into a bivariate analysis that allows quantification of the temporal coupling properties of two time series. We use eye-tracking data, 2,000 observations of six possible screen locations that are looked at by a dyad engaged into a joint task.  

``` r
listener = eyemovement$listener
narrator = eyemovement$narrator

delay = 1; embed = 1; rescale = 0; radius = .01;
normalize = 0; mindiagline = 2; minvertline = 2;
tw = 0; whiteline = FALSE; recpt = FALSE; side = "both"
method = 'crqa'; metric = 'euclidean';  
datatype = "categorical"

ans = crqa(narrator, listener, delay, embed, rescale, radius, normalize, 
           mindiagline, minvertline, tw, whiteline, recpt, side, method, metric, 
           datatype)
```

### Diagonal cross-recurrence profile. 

From cross-recurrence plots is possible to extract the diagonal cross-recurrence profiles (DCRPs) and use them to capture leader-follower-relationships. 

``` r
timecourse = round( seq(-3300,3300,33)/1000, digit = 2)  ## construct the time-course for the diagonal profile

res = drpfromts(narrator, listener, windowsize = 100,
                 radius = 0.001, delay = 1, embed = 1, rescale = 0,
                 normalize = 0, mindiagline = 2, minvertline = 2,
                 tw = 0, whiteline = F, recpt = F, side = 'both', 
                 method = 'crqa', metric = 'euclidean', 
                 datatype = 'continuous')
                 
 ## visualise the diagonal-profile
 profile = res$profile*100 ## extract it from the res object
 
plot(timecourse, profile, type = "l", lwd = 2.5, xlab = "Lag (seconds)",
     ylab = "Recurrence Rate %")                 
```

## Multidimensional cross-recurrence quantification analysis

Multidimensional cross-recurrence quantification analysis allows for the computation of cross-recurrences between two multidimensional time-series. We use hand-movement data from a complex LEGO joint construction task. The dataframe comprises of 5,799 observations.

``` r
# reduce the dimensionality of the time series to make the computation faster
# handset = handmovement[1:3000, ]
handset = handmovement[1:1000, ]

P1 = cbind(handset$P1_TT_d, handset$P1_TT_n) 
P2 = cbind(handset$P2_TT_d, handset$P2_TT_n)

delay = 5; embed = 2; rescale = 0; radius = .1;
normalize = 0; mindiagline = 10; minvertline = 10;
tw = 0; whiteline = FALSE; recpt = FALSE; side = "both"
method = 'mdcrqa'; metric = 'euclidean';  
datatype = "continuous"

ans = crqa(P1, P2, delay, embed, rescale, radius, normalize, 
           mindiagline, minvertline, tw, whiteline, recpt, side, method, metric, 
           datatype)

RP = ans$RP
results = unlist(ans[1:10])
print(results)
```

## Authors

* **Moreno I Coco** - *role = [cre, aut]* - (moreno.cocoi@gmail.com)
* **Dan Mønster** - *role = [aut]* - (danm@econ.au.dk)
* **Giuseppe Leonardi** - *role = [aut]* - (g.leonardi@vizja.pl)
* **Rick Dale** - *role = [aut]* - (rdale@ucla.edu)
* **Sebastian Wallot** - *role = [aut] - (sebastian.wallot@ae.mpg.de)

## Acknowledgments

* **James D. Dixon** - *role = [ctb]* - (james.dixon@uconn.edu)
* **John  C. Nash** - *role = [ctb]* -  (nashjc@uottawa.ca)