Code for the Manuscript: ‘An Understanding of Principal Differential
Analysis’
================

<!-- ![](outputs/SHM/paper-plots/3d-phase-plane.pdf) -->

<img src="outputs/SHM/paper-plots/3d-phase-plane.png" width="1889" />

------------------------------------------------------------------------

# Welcome

This repository contains code and data for the manuscript ***‘An
Understanding of Principal Differential Analysis’*** by [Edward
Gunning](https://edwardgunning.github.io/) and [Giles
Hooker](http://www.gileshooker.com/) ([arXiv
link](https://arxiv.org/abs/2406.18484)).

------------------------------------------------------------------------

# Repository Structure

- 📂 [**code**](code/)
  - 📄 contains R scripts with analysis and functions for the paper –
    the simple harmonic motion (SHM) model, the van der Pol (VdP) model
    and the real data analysis of the runner’s centre of mass (COM).
  - 📂 [**paper-figures**](code/paper-figures/) R scripts to generate
    figures from the paper including the simple harmonic motion (SHM)
    model, the van der Pol (VdP) model and the real data analysis of the
    runner’s centre of mass (COM).
- 📂 [**data**](data/)
  - data for the analysis of running kinematics. We are very grateful to
    Prof. Kieran Moran for providing this dataset.

------------------------------------------------------------------------

# References

- Ramsay, J. O. (1996). Principal Differential Analysis: Data Reduction
  by Differential Operators. Journal of the Royal Statistical Society.
  Series B (Methodological), 58(3), 495–508.

------------------------------------------------------------------------

# Computing Information

``` r
R.version
```

    ##                _                           
    ## platform       aarch64-apple-darwin20      
    ## arch           aarch64                     
    ## os             darwin20                    
    ## system         aarch64, darwin20           
    ## status                                     
    ## major          4                           
    ## minor          4.1                         
    ## year           2024                        
    ## month          06                          
    ## day            14                          
    ## svn rev        86737                       
    ## language       R                           
    ## version.string R version 4.4.1 (2024-06-14)
    ## nickname       Race for Your Life

``` r
sessionInfo()
```

    ## R version 4.4.1 (2024-06-14)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Sonoma 14.5
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: Europe/Dublin
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_4.4.1    fastmap_1.2.0     cli_3.6.3         tools_4.4.1      
    ##  [5] htmltools_0.5.8.1 rstudioapi_0.16.0 yaml_2.3.8        rmarkdown_2.27   
    ##  [9] highr_0.11        knitr_1.47        xfun_0.45         digest_0.6.36    
    ## [13] rlang_1.1.4       png_0.1-8         evaluate_0.24.0

``` r
packageVersion(pkg = "fda")
```

    ## [1] '6.1.8'

``` r
packageVersion(pkg = "deSolve")
```

    ## [1] '1.40'
