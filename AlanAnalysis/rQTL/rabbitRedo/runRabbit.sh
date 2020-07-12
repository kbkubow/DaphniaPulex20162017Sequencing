#!/usr/bin/env bash
#
#SBATCH -J trioPhase_whatshapp # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH --cpus-per-task=1 ### standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 1:00:00 # Running time of 1 hours
#SBATCH --mem 12G # Memory request of 8 GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/trioPhase_whatshapp.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/trioPhase_whatshapp.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### run as
# sbatch --array=1-12 ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/rabbitRedo/runRabbit.sh


module load intel/18.0 intelmpi/18.0 R/3.6.3
module load mathematica/11.1.1
module load parallel

#SLURM_ARRAY_TASK_ID=2

wd="/scratch/aob2x/daphnia_hwe_sims"
datadir="/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase"
chr=$( grep "^${SLURM_ARRAY_TASK_ID}" ${datadir}/chrs.csv | cut -f2 -d',' )
echo $chr

### generate RABBIT input data
Rscript ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/rabbitRedo/formatData.R ${chr}


### format RABBIT script file
RABBITpackageLocation="/scratch/aob2x/daphnia_hwe_sims/RABBIT/RABBIT_Packages/"
RABBITmodel="jointModel"
RABBITestfun="origViterbiDecoding"
generations_for_RABBIT=1
topDirectory="/scratch/aob2x/daphnia_hwe_sims/trioPhase"

# Generate mathematica script for RABBIT
python - <<EOF > ${datadir}/${chr}.RABBIT.m
print """SetDirectory["%s"]""" % "${RABBITpackageLocation}"
print """Needs["MagicReconstruct\`"]"""
print """Needs["MagicMap\`"]"""
print """SetDirectory["%s"]""" % "${datadir}"
print """popScheme = popScheme="ped.ped" """
print 'model = "%s"' % "${RABBITmodel}"
print 'estfun = "%s"' % "${RABBITestfun}"
print 'inputfile = "%s"' % "${chr}.all.in"
print 'resultFile = "%s.txt"' % "${chr}.out"
print """magicImpute[inputfile, model, popScheme, isFounderInbred -> False, outputFileID -> resultFile, isPrintTimeElapsed -> True]"""
print 'imputed = "%s"' % "${chr}.out_ImputedGenotype.csv"
print """magicReconstruct[imputed, model, popScheme, isFounderInbred -> False, outputFileID -> resultFile, reconstructAlgorithm -> estfun, isPrintTimeElapsed -> True]"""
print """summaryFile = StringDrop[resultFile, -4] <> ".csv" """
print """saveAsSummaryMR[resultFile<>"_magicReconstruct.txt", summaryFile]"""
print 'Exit'
EOF


### Run RABBIT impute & reconstruct
math -script ${datadir}/${chr}.RABBIT.m

### convert paths
python /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/parseHaplotypes.py \
${chr}.all.csv \
/scratch/aob2x/ >> ${chr}.all.haps
