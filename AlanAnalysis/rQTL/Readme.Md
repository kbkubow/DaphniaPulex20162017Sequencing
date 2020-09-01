### Run RABBIT to reconstruct phase of F1s
  #git clone https://github.com/chaozhi/RABBIT.git
  Apparently does not need phased parents, orininally built for only phenotyped F1s. Now it does all F1s

  ### Run RABBIT
  `rabbit.runRabbit.sh` runs RABBIT
      (requires `template.m`, runs `rabbit.formatData.R` & `rabbit.parseHaplotypes.py`, and files `chrs.csv`)

  `rabbit.plotRabbit.R` makes some diagnostic analysis

  ### convert RABBIT to
  `rabbit.convert_output.R` makes files for import into rQTL
    outputs to `/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm_onlyPheno`

### Run rQTL
`rqtl.R` Runs rQTL on phenotype data (needs to be worked on a bit)













### Phase A & C parents using quartet phasing with Whatshap
`whatshap.prepare_vcf.R` makes the input VCF file
`whatshap.run.sh` runs Whatshap to make output
