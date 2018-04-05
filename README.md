#### "Born to run? Quantifying the balance of prior bias and new information in prey escape decisions"

Sutton, NM and JP O'Dwyer

5 April 2018

This readme file describes the data files and R script accompanying the above publication. For further questions please contact nmsutto2@illinois.edu

## Data file: "BornToRun_DATA.txt"

Tab delimited text file containing angle of approach (in both degrees and radians), flight-initiation distances (meters), alert distances (meters), and site for human approaches with white-tailed deer. For use with our R script, header names must not be changed (i.e. when using your own data, please use the same headings we have).

Key for headings:

**FullANGLE**: angle with which human approached deer in degrees

**FullRADIANS**: angle with which human approached deer in radians

**FullFID**: flight inition distance of deer in meters (i.e. distance between approacher and prey when prey fled)

**FullAD**: alert distance of deer in meters (i.e. distance between approacher and prey when prey became aware of approacher)

**Site**: site/population where deer were encountered

## R script: "BornToRun_CODE.r"

R script for inferring prey prior distirbutions, analyzing goodness of fit of models, and comparing models based off of different decision mechanisms/risk factors. Steps for running the script and explanation of outputs are as follows:

*Step 1*: The following R packages must be installed:

* rootSolve
* hypergeo
* numDeriv
* optimx
* gsl

*Step 2*: Specify the path to a tab delimited data file, such as a text file, by copying the path in the quotation marks on line 15

*Step 3*: Lines 21-36 set species-specific energetic parameters that are held constant in our models. If not analyzing deer behavior,
	please adjust these parameters as needed depending on the species being considered (note that within species, results are typically
	qualitatively insensitive to the choice of these parameters).

*Step 4*: On line 42, specify the name of the site/population you wish to analyze (multiple sites cannot be run simultaneously,
	and must be analyzed one at a time).

*Step 5*: On line 50, specify the number of data sets to be generated during the goodness of fit test (data sets are randomnly
	generated by drawing from inferred distributions as part of an exact test for model goodness of fit). Smaller values will run
	much faster, but will sacrifice accuracy. Default is 1000, which will take a few minutes to run.

*Step 6*: After setting a value in step 5, the code is ready to run, and the entire script may be sourced.

## Description of outputs

Following is a copy of the outputs of one run of this script, with annotations added for explanation:

The script will first infer maximum likelihood estimates for shape parameters of the beta distribution for both the one and two risk factor models, and will output these parameters as follows

`Maximizing -- use negfn and neggr`
`Inferred shape parameters for two risk factor model:  1.621205877 1.932032681`

`Maximizing -- use negfn and neggr`
`Inferred shape parameters for one risk factor model 0.6671999957 11.70804459 `


After inferring shape parameters, the script will perform the goodness of fit test and report a p-value. Models pass this exact test when the log likelihood for observed data lies within the bulk of the distribution of log likelihoods of randomnly generated data sets (formed by random sampling from inferred beta distributions). In other words, if 0.05<p<0.95, the model passes the gof test

`Two risk factor goodness of fit: p = 0.6`

`One risk factor goodness of fit: p = 0.2`


Finally, the script will output the result of Vuong's closeness test, which tests whether the two risk factor model is better than the one risk factor model. If p<0.05, the two risk factor model is significantly closer to the true data generating process than the one risk factor model

`Vuong's closeness test (two versus one risk factor model): p = 0.009543313238`
