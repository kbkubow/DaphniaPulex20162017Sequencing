Analyses of Hardy-Weinberg deviations and attempt to cluster based on LD. They need to be run in roughly this order for one file to generate the samples for the next script.

> `HWE_simulations.*.*.R`
> These scripts perform analysis showing there is an excess of SNPs that are AA or Aa, consistent with a sex-chromosome like model

> `LD_clump.*.*.R`
> These scripts try to clump SNPs into linkage blocks focusing on those SNPs that are close to the peak of the ternary histograms. There is a significant overlap of ZW-like regions between D8 and DBunk. I wonder if this is due to shared history?

> `A_vs_B.analysis.workstation.R`
> Where do A & B superclones fall in their genotype distribution at these ZW like regions? A is over-represented for Aa, B is over-represented for AA.

> `association_test.workstation.R`
> What does the old phenotyping data look like?
