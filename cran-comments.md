## Resubmission

This is a resubmission. I have addressed the folloiwtn cCRAN comments:

* **Warning: Unexecutable code**: I have removeed the phrase "See vignette" from the example blocks for the functions comp_assemble, comp_comparisons, comp_make_f_tests. comp_means, comp_prepare_otu_tables, comp_prepare_phyloseq and make_comparisons.

* **Missing Rd-tags**: I have added the following for the return value of the function clear_warnings: "No return value, called for side effects. Clears the stored warning list so that stale warnings are no longer reported in the console."

* **Unwrap examples**: All functions now have examples that run in less than 5 seconds. None are wrapped in \dontrun.

* **Do not write information to the console that cannot be easily suppressed**: 
check_variance now returns a data frame of class check_var which has its own print.check_var function.

* **Do not set a seed**: I removed the line setting a seed from the function generate_passwords.

## Original submission

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

* NOTE: "Availability using Additional_repositories specification"
  Several imported packages (phyloseq, Biostrings, dada2, ShortRead) are
  hosted on Bioconductor rather than CRAN. The Bioconductor repository is
  listed in `Additional_repositories` in DESCRIPTION. All packages are
  available at https://bioconductor.org/packages/release/bioc.

## Test environments

* Windows 11, R 4.5.3 (local)
* R-hub (ubuntu-latest, R-devel) — via `rhub::rhub_check()`

## Downstream dependencies

There are no downstream dependencies on CRAN.
