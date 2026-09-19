# QsRutils 0.3.0

# QsRutils 0.3.0 Adds the following functions:

* extract_adapter_content() - Extracts information on adapter content from FastQC zip files allowing plotting of the proportion of each adapter versus the position in the read using ggplot2.

* find_truncation_parameters() - Finds the position in amplicon reads where the 75th percentile of the Q scores first falls below a Q value of 20.

* avg_dist() - Modification of vegan::avgdist() to include option for parallelization.  

# QsRutils 0.2.1

* Addressed CRAN comments re version 0.2.0.

# QsRutils 0.2.0

* `avg_alpha()` now uses a compiled C++ rarefaction routine (`rrarefy_cpp`) in place of `vegan::rrarefy()`, yielding roughly a 34% end-to-end speedup.
