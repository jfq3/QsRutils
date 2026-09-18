## Test environments
* local Windows 11, R 4.6.1
* Via GitHub Actions:
  - macos-latest,   r: 'release'
  - windows-latest, r: 'release'
  - ubuntu-latest,  r: 'devel'
  - ubuntu-latest,  r: 'release'
  - ubuntu-latest,  r: 'oldrel-1'
          
## R CMD check results
There were 0 ERRORS, 0 WARNINGS, and 0 NOTES.


## Submission summary
This is a update release that adds the following functions:
* extract_adapter_content() - Extracts information on adapter content from FastQC zip files allowing plotting of the proportion of each adapter versus the position in the read using ggplot2.

* find_truncation_parameters() - Finds the position in amplicon reads where the 75th percentile of the Q scores first falls below a Q value of 20.

* avg_dist() - Modification of vegan::avgdist() that includes an option for parallelization.  

## Downstream dependencies

There are no downstream dependencies on CRAN.
