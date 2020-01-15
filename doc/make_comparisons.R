## ---- tidy=TRUE---------------------------------------------------------------
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(QsRutils))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
data("its.root")
its.root

## ---- tidy=TRUE---------------------------------------------------------------
comp.phyla <- make_comparisons(its.root, taxrank = "Phylum", grps = "Label", transformation = "sqrt_arc_sine", pc.filter = 0.01, p.adjust.method ="BH", pool.sd = TRUE)
comp <- comp.phyla$comparison.table
comp <- comp[order(comp$mean, decreasing = TRUE), ]
comp

## -----------------------------------------------------------------------------
check_var(comp.phyla$taxa.pc.transformed, comp.phyla$groups)

## ---- tidy=TRUE---------------------------------------------------------------
asco <- subset_taxa(its.root, Phylum == "Ascomycota")
asco <- prune_taxa(taxa_sums(asco)>0, asco)
comp.asco <- make_comparisons(asco, taxrank = "Class", grps = "Label", transformation = "sqrt_arc_sine", pc.filter = 0.01, p.adjust.method ="BH", pool.sd = TRUE)
comp.asco$comparison.table

## ---- tidy=TRUE---------------------------------------------------------------
expt <- its.root
expt.pc <- transform_sample_counts(expt, function(x) 100*(x/sum(x)))

# Subset to higher rank
expt.taxon <- subset_taxa(expt, Phylum == "Ascomycota")
expt.taxon.pc <- subset_taxa(expt.pc, Phylum == "Ascomycota")

# Agglomerate to desired rank.
taxrank <- "Class"
expt.taxon <- tax_glom(expt.taxon, taxrank)
expt.taxon.pc <- tax_glom(expt.taxon.pc, taxrank)

# Remove ranks other than taxrank
tax_table(expt.taxon) <- tax_table(expt.taxon)[ , taxrank]
tax_table(expt.taxon.pc) <- tax_table(expt.taxon.pc)[ , taxrank]

# Filter out taxa that are < 0.1% of the total sequences in expt.
pc.filter <- 0.001
n <- sum(taxa_sums(expt)) * pc.filter
expt.taxon <- prune_taxa(taxa_sums(expt.taxon)>=n, expt.taxon)
expt.taxon.pc <- prune_taxa(taxa_names(expt.taxon), expt.taxon.pc)

# Make comparisons
temp2 <- comp_prepare_otu_table(expt.taxon.pc, grps = "Label", transformation = "sqrt_arc_sine")
temp3 <- comp_means_sd(temp2$otu.pc)
temp4 <- comp_make_f_tests(temp2$otu.pc.trans, grps = temp2$groups, var.equal = TRUE)
temp5 <- comp_comparisons(otu.pc = temp2$otu.pc, otu.pc.trans = temp2$otu.pc.trans, grps = temp2$groups, p.adjust.method = "BH",  pool.sd = TRUE)
comp <- comp_assemble(temp3, temp4, temp5)
comp <- comp[order(comp$mean, decreasing = TRUE), ]
comp

## ---- tidy=TRUE, fig.width=5, fig.align='center'------------------------------
data("plot_df")
ggplot(data = plot_df, aes(x = Treatment, y = Percent, fill = Family)) +
  stat_summary(fun.y = "mean", geom = "col", position = position_stack())

## ---- tidy=TRUE---------------------------------------------------------------
head(plot_df)

## ---- tidy=TRUE---------------------------------------------------------------
plot_df$seq <- with(plot_df, ave(Percent, Treatment, Family, FUN = seq_along))
my.data <- dcast(seq + Treatment ~ Family, data = plot_df, value.var = "Percent")
head(my.data)

## ---- tidy=TRUE---------------------------------------------------------------
my.grps <- my.data$Treatment
my.pc <- my.data[ , -c(1,2)]

## ---- tidy=TRUE---------------------------------------------------------------
my.pc.trans <- apply(my.pc, 2, sqrt_arc_sine)

## ---- tidy=TRUE---------------------------------------------------------------
my.pc <- t(my.pc)
my.pc.trans <- t(my.pc.trans)

## ---- tidy=TRUE---------------------------------------------------------------
temp3 <- comp_means_sd(my.pc)
temp4 <- comp_make_f_tests(my.pc.trans, grps = my.grps, var.equal = TRUE)
temp5 <- comp_comparisons(otu.pc = my.pc, otu.pc.trans = my.pc.trans, grps = my.grps, p.adjust.method = "BH",  pool.sd = TRUE)
comp <- comp_assemble(temp3, temp4, temp5)
comp <- comp[order(comp$mean, decreasing = TRUE), ]
comp

