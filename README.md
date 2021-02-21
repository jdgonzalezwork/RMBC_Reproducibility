Summary
-------

This repository is for reproducing the analysis on real data of the
article *Robust Model-Based* Clustering whose authors are Juan D.
Gonzalez, [Ricardo
Maronna](https://scholar.google.com/citations?user=Cqa77SYAAAAJ&hl=en),
[Victor J.
Yohai](https://scholar.google.com/citations?user=fBUfRqcAAAAJ&hl=en),
and [Ruben H.
Zamar](https://scholar.google.ca/citations?user=XgRp4HUAAAAJ&hl=en).

The only thing you should do in order to reproduce the results is to
download the repository and run in R

    source("Fig2_and_Table3.R")

This command will do the following steps

-   Install the package “RMBC” from github. This packages implements the
    estimator developed in \[1\] and contains the real data set used to
    asses our procedure.
-   Install the necessary packages from CRAN, namely
    -   [`tclust`](https://cran.r-project.org/web/packages/tclust/index.html)
    -   [`RSKC`](https://cran.r-project.org/web/packages/RSKC/index.html)
    -   [`GSE`](https://cran.r-project.org/web/packages/GSE/index.html)
    -   [`otrimle`](https://cran.r-project.org/web/packages/otrimle/index.html)
    -   [`mclust`](https://cran.r-project.org/web/packages/mclust/index.html)
    -   [`mvtnorm`](https://cran.r-project.org/web/packages/mvtnorm/index.html)
    -   [`ktaucenters`](https://cran.r-project.org/web/packages/ktaucenters/index.html)
    -   [`combinat`](https://cran.r-project.org/web/packages/combinat/index.html)
-   Run the auxiliary routines that aid to plot the results, compute
    performance measures and give to the different estimators the same
    format.
    -   `auxComputeClusters.R`
    -   `auxPlotGroupsModelBasedCluster.R`
    -   `performance_measures.R`

**Reference**:  
\[1\] Gonzalez, J. D., Maronna, R., Yohai, V. J., & Zamar, R. H. (2021).
Robust Model-Based Clustering. arXiv preprint arXiv:2102.06851.
