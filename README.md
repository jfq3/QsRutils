<!-- badges: start -->
[![R-CMD-check](https://github.com/jfq3/QsRutils/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jfq3/QsRutils/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->


QsRutils - R Functions Useful for Community Ecology
===================================================

Well, at least most of the functions are. There are some other things here, too, including a password generator. For several years I collected functions I wrote and used often into a package I called MyRutils. But I really can't post a package with that name, can I? It just doesn't make sense if you are going to use it. So I changed the name to QsRutils. After all, it fits because decades ago a student gave me the nickname Q. Knowing she was a Trekkie, I let it slide.

The functions here are a bit of a mish-mash. Several are related to phyloseq and vegan, packages I use extensively. veganotu extracts a phyloseq otu\_table as a matrix with samples as rows and taxa as columns, the format used by vegan functions. vegansam extracts a phyloseq sample\_data object as a data frame while preserving all data types, again ready to use in vegan. vegan\_stand applies any vegan decostand transformation to the otu\_table in a phyloseq object.

Several other functions are related to making tables comparing taxa abundances among treatments, as in the supplementary tables in Zhang, *et.al.*, 2017. A vignette for creating such table is included in this package.

Not all functions from MyRutils have made it here yet, so you can expect that more will be added as time goes by, especially any I use in my workshops, tutorials and blogs at john-quensen.com.

Version 0.1.3 adds a function to root trees in phyloseq objects by their longest terminal branch.

Version 0.1.4 adds a function to get ggplot2 plot limits.

Version 0.1.5 adds function srs_p to normalize the OTU table in a phyloseq object. This is an alternative to "rarefying" to the same number of counts per sample.

Version 0.2.1 adds several functions.  
  -  avg_alpha - Calculates Shannon, Observed, Pielou, Simpson and Inverse Simpson alpha-diversity metrics as the mean or median of repetitive samplings of the OTU table.  
  - check_primer_hits - Determines hits of all orientations of the primers to paired sequence files. Used to determine if merging will create overhangs.  
 - format_taxon - Adds *'s around proper parts of taxon names so that the names can be rendered in italics by Rmarkdown. Useful in making ggplots.  
  - hash_dna_seqs - Converts DNA seequences into hashes encoding the sequences. Useful in shortening the taxa names especially when using the R verision of DADA2 which outputs the OTU column names as the sequences themselves.  
  - plot_transition_stats - Makes a plot of DADA2 transition rates from the transition stats qza file output by QIIME2 DADA2 beginning with QIIME 2 version 2025.7. Useful in determing how well DADA2 corrected sequence errors.  
  - se - calculates the standard eror.  
  - %wo% - given vectors x and y, returns elements of x that are not in y.  

Installation
------------
Install from CRAN with install.packages("QsRutils")  

You can install the development version of QsRutils from github with:

``` r
install.packages("pak")
pak::pak("jfq3/QsRutils")
```

References
----------

Zhang, B., Penton, C.R., Xue, C., Quensen, J.F., Roley, S.S., Guo, J., Garoutte, A., Zheng, T., Tiedje, J.M., 2017. Soil depth and crop determinants of bacterial communities under ten biofuel cropping systems. Soil Biology & Biochemistry 112, 140-152.
