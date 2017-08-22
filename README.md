
<!-- README.md is generated from README.Rmd. Please edit that file -->
QsRutils - R Functions Useful for Community Ecology
===================================================

Well, at least most of the functions are. There are some other things here, too, including a password generator. For several years I collected functions I wrote and used often into a package I called MyRutils. But I really can't post a package with that name, can I? It just doesn't make sense if you are going to use it. So I changed the name to QsRutils. After all, it fits because decades ago a student gave me the nickname Q. Knowing she was a Trekkie, I let it slide.

The functions here a a bit of a mish-mash. Several are related to phyloseq and vegan, packages I use extensively. veganotu extracts a phyloseq otu\_table as a matrix with samples as rows and taxa as columns, the format used by vegan functions. vegansam extracts a phyloseq sample\_data object as a data frame while preserving all data types, again ready to use in vegan. vegan\_stand applies any vegan decostand transformation to the otu\_table in a phyloseq object.

Several other functions are related to making tables comparing taxa abundances among treatments, as in the supplementary tables in Zhang, *et.al.*, 2017. A vignette for creating such table is included in this package.

Not all functions from MyRutils have made it here yet, so you ca expect that more will be added as time goes by, especially any I use in my workshops, tutorials and blogs at john-quensen.com.

Installation
------------

You can install QsRutils from github with:

``` r
# install.packages("devtools")
devtools::install_github("jfq3/QsRutils")
```

References
----------

Zhang, B., Penton, C.R., Xue, C., Quensen, J.F., Roley, S.S., Guo, J., Garoutte, A., Zheng, T., Tiedje, J.M., 2017. Soil depth and crop determinants of bacterial communities under ten biofuel cropping systems. Soil Biology & Biochemistry 112, 140-152.
