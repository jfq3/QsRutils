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

## Response to CRAN check issues

Submitting QsRutils version 0.3 resulted in 1 ERROR and 2 NOTES.

The previous submission (v0.2.1) produced an ERROR on r-oldrel-macos-arm64 because
the `dada2` Bioconductor package was not available on that platform. This has been
addressed by moving `dada2` from `Imports` to `Suggests` and adding a
`requireNamespace("dada2")` check inside `check_primer_hits()`, which is the only
function that calls `dada2`. Users on platforms where `dada2` is unavailable can
install the package normally; they will receive an informative error only if they
call `check_primer_hits()`.

One NOTE was the result of the Additional-repositories line in the DESCRIPTION. This is an informational note only.

The second NOTE was because the avg_dist() example took more than 5 seconds elapsed time. 
I changed the example parameters so that the example completes in < 5 seconds.

## Downstream dependencies

There are no downstream dependencies on CRAN.
