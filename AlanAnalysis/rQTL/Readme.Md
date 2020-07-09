### first phase A & C parents using quartet phasing
`prepare_vcf_rQTL.R` makes the input VCF file
`whatshap_rQTL.sh` runs Whatshap to make output

### prepare input data for RABBIT
`formatData.R` takes input of A&C haplotypes plus GDS to make rabbit file

### run RABBIT
#git clone https://github.com/chaozhi/RABBIT.git


eps="0.005"
epsF="0.005"
RABBITmodel="jointModel"
RABBITestfun="origViterbiDecoding"
RABBITpackageLocation="/scratch/aob2x/daphnia_hwe_sims/RABBIT/RABBIT_Packages/"
outdir="/scratch/aob2x/daphnia_hwe_sims/rabbit_out"
generations_for_RABBIT="1"

python - <<EOF > /scratch/aob2x/daphnia_hwe_sims/trioPhase/AxC_Rabbit.m
print """SetDirectory["%s"]""" % "${RABBITpackageLocation}"
print """Needs["MagicReconstruct\`"]"""
print """SetDirectory["%s"]""" % "${outdir}"
print """popScheme = Table["RM1-E", {%s}]""" % "${generations_for_RABBIT}"
print 'epsF = %s' % "${epsF}"
print 'eps = %s' % "${eps}"
print 'model = "%s"' % "${RABBITmodel}"
print 'estfun = "%s"' % "${RABBITestfun}"
print 'inputfile = "%s"' % "/scratch/aob2x/daphnia_hwe_sims/trioPhase/AxC_Rabbit.in"
print 'resultFile = "%s.txt"' % "/scratch/aob2x/daphnia_hwe_sims/trioPhase/AxC_Rabbit.out"
print """magicReconstruct[inputfile, model, epsF, eps, popScheme, resultFile, HMMMethod -> estfun, PrintTimeElapsed -> True]"""
print 'summaryFile = StringDrop[resultFile, -4] <> ".csv"'
print 'saveAsSummaryMR[resultFile, summaryFile]'
print 'Exit'
EOF

module load mathematica/11.1.1

math -noprompt -script /scratch/aob2x/daphnia_hwe_sims/trioPhase/AxC_Rabbit.m


SetDirectory["/scratch/aob2x/daphnia_hwe_sims/RABBIT/RABBIT_Packages/"]
Needs["MagicReconstruct`"]
SetDirectory["/scratch/aob2x/daphnia_hwe_sims/trioPhase/"]
popScheme = Table["RM1-E", {1}]
epsF = 0.005
eps = 0.005
model = "jointModel"
estfun = "origViterbiDecoding"
inputfile = "AxC_Rabbit.in"
resultFile = "AxC_Rabbit.out.txt"
magicReconstruct[inputfile, model, epsF, eps, popScheme, resultFile, HMMMethod -> estfun, PrintTimeElapsed -> True]
summaryFile = StringDrop[resultFile, -4] <> ".csv"
saveAsSummaryMR[resultFile, summaryFile]

Exit







wd=/scratch/aob2x/

git clone https://github.com/chaozhi/RABBIT.git
