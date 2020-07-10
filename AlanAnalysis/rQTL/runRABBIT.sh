#!/usr/bin/env bash
#
#SBATCH -J trioPhase_whatshapp # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH --cpus-per-task=1 ### standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 12:00:00 # Running time of 2 hours
#SBATCH --mem 12G # Memory request of 8 GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/trioPhase_whatshapp.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/trioPhase_whatshapp.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


module load gcc
module load R/3.6.3
module load mathematica

chromosome="${1}"
clone="${2}"

#clone="AxB_R4_P17_B"; chromosome="Scaffold_2373_HRSCAF_2879"

eps="0.005"
epsF="0.005"
RABBITmodel="jointModel"
RABBITestfun="origViterbiDecoding"
RABBITpackageLocation="/scratch/aob2x/daphnia_hwe_sims/RABBIT/RABBIT_Packages/"
generations_for_RABBIT=1
topDirectory="/scratch/aob2x/daphnia_hwe_sims/trioPhase"

# Generate mathematica script for RABBIT
python - <<EOF > /scratch/aob2x/daphnia_hwe_sims/trioPhase/rabbit_m/$chromosome.$clone.RABBIT.m
print """SetDirectory["%s"]""" % "${RABBITpackageLocation}"
print """Needs["MagicReconstruct\`"]"""
print """SetDirectory["%s"]""" % "${topDirectory}"
print """popScheme = Table["RM1-E", {%s}]""" % "${generations_for_RABBIT}"
print 'epsF = %s' % "${epsF}"
print 'eps = %s' % "${eps}"
print 'model = "%s"' % "${RABBITmodel}"
print 'estfun = "%s"' % "${RABBITestfun}"
print 'inputfile = "%s"' % "${topDirectory}/rabbitIn/${chromosome}.${clone}.in"
print 'resultFile = "%s.txt"' % "${topDirectory}/rabbitOut/${chromosome}.${clone}.out"
print """magicReconstruct[inputfile, model, popScheme, outputFileID -> resultFile, reconstructAlgorithm -> estfun, PrintTimeElapsed -> True]"""
print 'summaryFile = StringDrop[resultFile, -4] <> ".csv"'
print 'saveAsSummaryMR[resultFile, summaryFile]'
print 'Exit'
EOF

# Run RABBIT
  math -noprompt -script /scratch/aob2x/daphnia_hwe_sims/trioPhase/rabbit_m/${chromosome}.${clone}.RABBIT.m
