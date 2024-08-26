# Improving convergence of the joint model estimation procedure


Sometimes the joint model for stepped wedge cluster-randomized trials
with non-ignorable dropout fails to converge due to the increased
computational complexity. In this document, we illustrate several steps
and approaches that can be used to get the joint model to converge.

We start by setting the version of the Stata interpreter to version 18;
this means that the following code will continue to work in all future
versions of Stata, even if Stata were to change.

``` stata
version 18
```

<style>div.jp-Notebook .datagrid-container {min-height: 448px; }</style>

# Data

We use the following dataset for illustration purposes:

``` stata
use "04-data.dta", clear
codebook, compact
```


    Variable       Obs Unique      Mean        Min       Max  Label
    -------------------------------------------------------------------------------
    repn          8000      1         4          4         4  
    i             8000     32      16.5          1        32  
    id            8000   1600     800.5          1      1600  
    j             8000      5         3          1         5  
    jid           8000    160      80.5          1       160  
    x             8000      2        .5          0         1  
    cumx          8000      5         1          0         4  
    y             8000   8000  32.93547  -6.646588  75.60858  
    yobs          4360   4360  35.44544   4.923843  75.60858  
    t0            5335      5  1.585942          0         4  
    t             5335    980  2.489409   .0009573         5  
    d             5335      2  .1827554          0         1  
    eventtime     8000    976  3.012497   .0009573         5  
    actual_eve~e  8000   1600   7.52727   .0009573  117.5387  
    -------------------------------------------------------------------------------

Note that this dataset is one of the datasets that did not converge in
the simulation study with our default settings.

The relevant columns are:

- `i`, the cluster indicator variable;
- `id`, the participant indicator variable;
- `j`, the period indicator variable;
- `x`, the binary treatment indicator variable;
- `yobs`, the observed longitudinal outcome values;
- `t0` and `t`, the time to dropout process in start-stop notation;
- `d`, the dropout indicator variable.

# Default settings

Our experience suggests that convergence of the joint model improves
when using the `startgrid` and `technique(bfgs)` options. `datagrid`
instructs `gsem` to test a grid of possible starting values for the
random effects, while `technique(bfgs)` instructs `gsem` to use the
Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm instead of Stata’s
modified Newton-Raphson (NR) algorithm, the default.

The `gsem` code is as follows:

``` stata
capture noisily gsem ///
    (yobs <- ibn.j i.x M1[i]@1 M2[i>id]@1, noconstant family(gaussian)) ///
    (t <- i.x M1[i] M2[i>id], family(weibull, failure(d) lt(t0))) ///
    , startgrid technique(bfgs) nolog
```

    cannot compute an improvement -- flat region encountered

Note that we prefix the call to `gsem` with `capture noisily` to capture
the error (and continue compiling the document), and we use the `nolog`
option to omit the iterations log.

Using these settings, the BFGS algorithm fails to converge to a
solution.

# Using `gsem`’s default settings

The first step we can try consists of using the default `gsem` settings:

``` stata
capture noisily gsem ///
    (yobs <- ibn.j i.x M1[i]@1 M2[i>id]@1, noconstant family(gaussian)) ///
    (t <- i.x M1[i] M2[i>id], family(weibull, failure(d) lt(t0))) ///
    , nolog
```

    convergence not achieved

    Generalized structural equation model               Number of obs   =    5,335

    Response: yobs                                      Number of obs   =    4,360
    Family:   Gaussian            
    Link:     Identity            

    Response: t                                         Number of obs   =    5,335
    Family:   Weibull                                   No. of failures =      975
    Form:     Proportional hazards                      Time at risk    = 4,819.99
    Link:     Log                 

    Log likelihood = -17586.699

     ( 1)  [yobs]_cons = 0
     ( 2)  [yobs]M1[i] = 1
     ( 3)  [yobs]M2[i>id] = 1
    ------------------------------------------------------------------------------
                 | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
    -------------+----------------------------------------------------------------
    yobs         |
               j |
              1  |   30.21271   .3095227    97.61   0.000     29.60606    30.81936
              2  |   30.64101   .3576861    85.66   0.000     29.93996    31.34206
              3  |   30.28756     .41302    73.33   0.000     29.47805    31.09706
              4  |   30.22942   .4731771    63.89   0.000     29.30201    31.15683
              5  |   30.32292   .5287668    57.35   0.000     29.28655    31.35928
                 |
             1.x |   4.788127   .3469723    13.80   0.000     4.108074     5.46818
                 |
           M1[i] |          1  (constrained)
        M2[i>id] |          1  (constrained)
                 |
           _cons |          0  (omitted)
    -------------+----------------------------------------------------------------
    t            |
             1.x |  -.3384015   .0865106    -3.91   0.000    -.5079591   -.1688438
                 |
           M1[i] |  -1.91e+32          .        .       .            .           .
        M2[i>id] |  -.1187545   .0091306   -13.01   0.000    -.1366502   -.1008588
                 |
           _cons |  -1.569263   .0605129   -25.93   0.000    -1.687866   -1.450659
    -------------+----------------------------------------------------------------
    /t           |
            ln_p |   .0497505   .0372434                     -.0232452    .1227461
    -------------+----------------------------------------------------------------
       var(M1[i])|   4.10e-97   3.44e-82                             .           .
    var(M2[i>id])|   56.41505   3.528683                      49.90608    63.77294
    -------------+----------------------------------------------------------------
      var(e.yobs)|   40.51348   1.017213                      38.56804    42.55705
    ------------------------------------------------------------------------------
    convergence not achieved

This did not work either.

# Using the `difficult` option

The second step we can try is to use the `difficult` option, with both
BFGS and Stata’s modified Newton-Raphson. This instructs Stata to use a
different stepping algorithm in non-concave regions, which may help with
convergence.

``` stata
capture noisily gsem ///
    (yobs <- ibn.j i.x M1[i]@1 M2[i>id]@1, noconstant family(gaussian)) ///
    (t <- i.x M1[i] M2[i>id], family(weibull, failure(d) lt(t0))) ///
    , difficult nolog
```

    convergence not achieved

    Generalized structural equation model               Number of obs   =    5,335

    Response: yobs                                      Number of obs   =    4,360
    Family:   Gaussian            
    Link:     Identity            

    Response: t                                         Number of obs   =    5,335
    Family:   Weibull                                   No. of failures =      975
    Form:     Proportional hazards                      Time at risk    = 4,819.99
    Link:     Log                 

    Log likelihood = -17586.699

     ( 1)  [yobs]_cons = 0
     ( 2)  [yobs]M1[i] = 1
     ( 3)  [yobs]M2[i>id] = 1
    ------------------------------------------------------------------------------
                 | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
    -------------+----------------------------------------------------------------
    yobs         |
               j |
              1  |   30.21271   .3095228    97.61   0.000     29.60606    30.81936
              2  |   30.64102   .3576862    85.66   0.000     29.93996    31.34207
              3  |   30.28757   .4130201    73.33   0.000     29.47806    31.09707
              4  |   30.22943   .4731772    63.89   0.000     29.30202    31.15684
              5  |   30.32292   .5287669    57.35   0.000     29.28656    31.35929
                 |
             1.x |    4.78813   .3469723    13.80   0.000     4.108077    5.468183
                 |
           M1[i] |          1  (constrained)
        M2[i>id] |          1  (constrained)
                 |
           _cons |          0  (omitted)
    -------------+----------------------------------------------------------------
    t            |
             1.x |  -.3384061   .0865104    -3.91   0.000    -.5079635   -.1688488
                 |
           M1[i] |  -5.659763          .        .       .            .           .
        M2[i>id] |  -.1187531   .0091305   -13.01   0.000    -.1366486   -.1008576
                 |
           _cons |  -1.569246    .060512   -25.93   0.000    -1.687847   -1.450644
    -------------+----------------------------------------------------------------
    /t           |
            ln_p |   .0497447   .0372432                     -.0232506    .1227399
    -------------+----------------------------------------------------------------
       var(M1[i])|   4.89e-24          .                             .           .
    var(M2[i>id])|   56.41482   3.528658                       49.9059    63.77266
    -------------+----------------------------------------------------------------
      var(e.yobs)|   40.51347   1.017212                      38.56803    42.55704
    ------------------------------------------------------------------------------
    convergence not achieved

``` stata
capture noisily gsem ///
    (yobs <- ibn.j i.x M1[i]@1 M2[i>id]@1, noconstant family(gaussian)) ///
    (t <- i.x M1[i] M2[i>id], family(weibull, failure(d) lt(t0))) ///
    , difficult technique(bfgs) nolog
```

    cannot compute an improvement -- discontinuous region encountered

Using a different stepping algorithm did not lead to convergence, with
either of the methods.

# Custom starting values

The third step we can try is passing custom starting values to `gsem`.
To this end, we first fit the corresponding linear mixed effects model
(ignoring the dropout process):

``` stata
gsem (yobs <- ibn.j i.x M1[i]@1 M2[i>id]@1, noconstant family(gaussian))
```


    Fitting fixed-effects model:

    Iteration 0:  Log likelihood = -15825.816  
    Iteration 1:  Log likelihood = -15825.816  

    Refining starting values:

    Grid node 0:  Log likelihood = -15573.441

    Fitting full model:

    Iteration 0:  Log likelihood = -15573.441  (not concave)
    Iteration 1:  Log likelihood = -15288.349  
    Iteration 2:  Log likelihood = -15277.955  
    Iteration 3:  Log likelihood = -15248.326  
    Iteration 4:  Log likelihood = -15199.172  
    Iteration 5:  Log likelihood = -15195.574  
    Iteration 6:  Log likelihood = -15195.544  
    Iteration 7:  Log likelihood = -15195.544  

    Generalized structural equation model                    Number of obs = 4,360
    Response: yobs    
    Family:   Gaussian
    Link:     Identity
    Log likelihood = -15195.544

     ( 1)  [yobs]_cons = 0
     ( 2)  [yobs]M1[i] = 1
     ( 3)  [yobs]M2[i>id] = 1
    ------------------------------------------------------------------------------
                 | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
    -------------+----------------------------------------------------------------
    yobs         |
               j |
              1  |   31.76479   .2734336   116.17   0.000     31.22887    32.30071
              2  |   32.75502   .3099614   105.67   0.000     32.14751    33.36253
              3  |   32.67481   .3668142    89.08   0.000     31.95587    33.39375
              4  |   32.78311   .4317538    75.93   0.000     31.93689    33.62934
              5  |   32.97506   .4920062    67.02   0.000     32.01074    33.93937
                 |
             1.x |    4.75271   .3524938    13.48   0.000     4.061835    5.443585
                 |
           M1[i] |          1  (constrained)
        M2[i>id] |          1  (constrained)
                 |
           _cons |          0  (omitted)
    -------------+----------------------------------------------------------------
       var(M1[i])|   .1795762   .4252089                      .0017327    18.61071
    var(M2[i>id])|   44.17984   2.454326                      39.62208    49.26188
    -------------+----------------------------------------------------------------
      var(e.yobs)|   40.96164   1.037622                       38.9776    43.04667
    ------------------------------------------------------------------------------

This model converged without issues. Then, we save the parameter
estimates in a matrix named `bfrom`:

``` stata
matrix bfrom = e(b)
```

Finally, we can pass these starting values to `gsem` using the `from`
option:

``` stata
gsem ///
    (yobs <- ibn.j i.x M1[i]@1 M2[i>id]@1, noconstant family(gaussian)) ///
    (t <- i.x M1[i] M2[i>id], family(weibull, failure(d) lt(t0))) ///
    , from(bfrom) startgrid technique(bfgs)
```


    Fitting fixed-effects model:

    Iteration 0:  Log likelihood = -20425.271  
    Iteration 1:  Log likelihood = -18341.449  
    Iteration 2:  Log likelihood = -18324.004  
    Iteration 3:  Log likelihood = -18323.957  
    Iteration 4:  Log likelihood = -18323.957  

    Refining starting values:

    Grid node 0:  Log likelihood = -17899.457
    Grid node 1:  Log likelihood = -19054.862
    Grid node 2:  Log likelihood = -18882.525
    Grid node 3:  Log likelihood = -18158.564
    Grid node 4:  Log likelihood = -19018.855
    Grid node 5:  Log likelihood = -18852.216
    Grid node 6:  Log likelihood = -18152.002
    Grid node 7:  Log likelihood = -19032.516
    Grid node 8:  Log likelihood = -18866.722
    Grid node 9:  Log likelihood = -18170.503

    Fitting full model:

    Iteration 0:  Log likelihood = -17899.457  
    Iteration 1:  Log likelihood =  -17697.12  (backed up)
    Iteration 2:  Log likelihood =  -17689.09  (backed up)
    Iteration 3:  Log likelihood = -17672.976  (backed up)
    Iteration 4:  Log likelihood =  -17654.24  
    Iteration 5:  Log likelihood = -17649.427  
    Iteration 6:  Log likelihood =  -17641.68  
    Iteration 7:  Log likelihood = -17634.953  
    Iteration 8:  Log likelihood = -17632.377  (backed up)
    Iteration 9:  Log likelihood = -17625.577  
    Iteration 10: Log likelihood = -17623.274  
    Iteration 11: Log likelihood = -17615.286  
    Iteration 12: Log likelihood = -17608.875  
    Iteration 13: Log likelihood = -17605.247  
    Iteration 14: Log likelihood = -17602.121  
    Iteration 15: Log likelihood = -17591.085  
    Iteration 16: Log likelihood = -17587.817  
    Iteration 17: Log likelihood = -17586.995  
    Iteration 18: Log likelihood = -17586.382  
    Iteration 19: Log likelihood = -17585.699  
    Iteration 20: Log likelihood = -17585.607  
    Iteration 21: Log likelihood = -17585.605  
    Iteration 22: Log likelihood = -17585.605  
    Iteration 23: Log likelihood = -17585.605  
    Iteration 24: Log likelihood = -17585.604  
    Iteration 25: Log likelihood = -17585.604  
    Iteration 26: Log likelihood = -17585.604  

    Generalized structural equation model               Number of obs   =    5,335

    Response: yobs                                      Number of obs   =    4,360
    Family:   Gaussian            
    Link:     Identity            

    Response: t                                         Number of obs   =    5,335
    Family:   Weibull                                   No. of failures =      975
    Form:     Proportional hazards                      Time at risk    = 4,819.99
    Link:     Log                 

    Log likelihood = -17585.604

     ( 1)  [yobs]_cons = 0
     ( 2)  [yobs]M1[i] = 1
     ( 3)  [yobs]M2[i>id] = 1
    ------------------------------------------------------------------------------
                 | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
    -------------+----------------------------------------------------------------
    yobs         |
               j |
              1  |    30.2072   .3383947    89.27   0.000     29.54395    30.87044
              2  |   30.63811   .3829617    80.00   0.000     29.88752     31.3887
              3  |   30.28871   .4353314    69.58   0.000     29.43548    31.14194
              4  |   30.23615   .4931733    61.31   0.000     29.26955    31.20276
              5  |   30.33463   .5471961    55.44   0.000     29.26215    31.40712
                 |
             1.x |   4.767812   .3491602    13.66   0.000      4.08347    5.452153
                 |
           M1[i] |          1  (constrained)
        M2[i>id] |          1  (constrained)
                 |
           _cons |          0  (omitted)
    -------------+----------------------------------------------------------------
    t            |
             1.x |  -.3339431   .0872295    -3.83   0.000    -.5049097   -.1629764
                 |
           M1[i] |  -.1618064   .0696937    -2.32   0.020    -.2984036   -.0252093
        M2[i>id] |   -.118595   .0091703   -12.93   0.000    -.1365685   -.1006216
                 |
           _cons |  -1.571052   .0645307   -24.35   0.000     -1.69753   -1.444575
    -------------+----------------------------------------------------------------
    /t           |
            ln_p |   .0500247   .0372978                     -.0230776     .123127
    -------------+----------------------------------------------------------------
       var(M1[i])|   .6116154   .5639303                      .1003787    3.726622
    var(M2[i>id])|   55.86785   3.515298                       49.3859    63.20057
    -------------+----------------------------------------------------------------
      var(e.yobs)|    40.5172   1.017378                      38.57145    42.56111
    ------------------------------------------------------------------------------

Finally, this joint model converged to a solution.

# Other options

If the joint model did not converge after trying the steps above, there
is still hope. Specifically, we suggest:

1.  Testing different maximisation algorithms (e.g., passing different
    options to the `technique` argument). Please consult the `gsem`
    manual to find all available options;

2.  Testing different combinations of the options introduced above;

3.  Testing a different software implementation, e.g., using the
    user-written `merlin` command (more details
    [here](https://arxiv.org/abs/1806.01615v1)).

If the joint model still does not converge after trying all these steps,
then one may need to simplify the model, e.g., because some coefficients
are not identifiable using the available data.
