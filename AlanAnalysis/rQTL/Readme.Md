### Phase A & C parents using quartet phasing with Whatshap
`whatshap.prepare_vcf.R` makes the input VCF file
`whatshap.run.sh` runs Whatshap to make output

### Run RABBIT to reconstruct phase of F1s
#git clone https://github.com/chaozhi/RABBIT.git

`rabbit.formatData.R` takes input of A&C haplotypes plus GDS to make rabbit data files
`rabbit.runRabbit.sh` runs RABBIT (requires `template.m`)
`rabbit.plotRabbit.R` makes some diagnostic analysis
`rabbit.convert_output.R` makes output files for import into rQTL

### Run rQTL
`rqtl.R` Runs rQTL on phenotype data
