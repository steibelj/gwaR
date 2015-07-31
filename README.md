# gwaR
This package implements GWA from GBLUP model following Gualdron Duarte et al (2014) - http://www.biomedcentral.com/1471-2105/15/246/abstract

The packagte is under developement, but it is functional.

To install use: `devtools::install_github("steibelj/gwaR")`


A quick introduction: `vignette("intro-to-gblup")`

###Future commits:
1) ONGOING: improve documentation to expand to other capabilities

2) DONE: improve plot functions: better graphics

3) include predict method - DONE for training set (summary.gblup) PENDING: validation set 

4) implement multi-core computation

5)  DONE: LRT functions for forward, backward comparison

6) implement high-throughput functions for massively paralell phenotypes

7) include functions for GWA peak annotation from genomic dBs

8) make memory efficient for extra-large systems of equations (seveal 10000s animals)

9) implement META-ANALYSIS according to Yeni's paper

10) implement shrinkage tests according to Misztal's groups papers
