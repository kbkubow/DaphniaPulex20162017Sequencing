Analyses of Hardy-Weinberg deviations and attempt to cluster based on LD. They need to be run in roughly this order for one file to generate the samples for the next script.

> `HWE_simulations.*.*.R`
> These scripts perform analysis showing there is an excess of SNPs that are AA or Aa, consistent with a sex-chromosome like model

> `LD_clump.*.*.R`
> These scripts try to clump SNPs into linkage blocks focusing on those SNPs that are close to the peak of the ternary histograms. There is a significant overlap of ZW-like regions between D8 and DBunk. I wonder if this is due to shared history?

> `getTargetSNPs.analysis.workstation.R`
> generates list of ZW-like SNPs to use for downstream analyses; based on single pops, intersections, and annotated with linkage groups, as defined above. Takes input from LD_clump, makes snpSet.Rdata object

> `TwoCloneTest.analysis.workstation.R`
> Generalization of A vs B test (and it's derivative, Alike-Blike); generalized a bit to work for D8 & DCat (or any two clones), different input SNP lists, and works along chromosomes as well as across genome.

> `TwoCloneTest_pairwise.analysis.workstation.R`
> Even further generalization of A vs B test (and it's derivative, Alike-Blike); generalized a bit to work for D8 & DCat (or any two clones), different input SNP lists, and works along chromosomes as well as across genome.


> `association_test.workstation.R`
> What does the old phenotyping data look like? Look first at ZW-like SNPs in the


TO DO:
  • DCat test of AxB script
  • Redo Dcat & D8 with Dbunk SNP list
  • Test reference genome, Karen will send bam file path


### Defunct scripts

  > `DEFUNCT_D8_A_vs_B.analysis.workstation.R`
  > Where do A & B superclones fall in their genotype distribution at these ZW like regions? A is over-represented for Aa, B is over-represented for AA.

  > `DEFUNCT_Dcat_Alike_vs_Blike.analysis.workstation.R`
  > Where do the two superclones fro DCat which are A-like and-Blike fall in their genotype distribution at these ZW like regions? A is over-represented for Aa, B is over-represented for AA.
