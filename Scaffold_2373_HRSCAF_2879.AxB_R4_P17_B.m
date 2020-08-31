SetDirectory["/scratch/aob2x/daphnia_hwe_sims/RABBIT/RABBIT_Packages/"]
Needs["MagicReconstruct`"]
SetDirectory["/scratch/aob2x/daphnia_hwe_sims/trioPhase/rabbitIn"]
popScheme = Table["RM1-E", {1}]
epsF = 0.005
eps = 0.005
model = "indepModel"
estfun = "origViterbiDecoding"
inputfile = "Scaffold_2373_HRSCAF_2879.AxB_R4_P17_B.in"
resultFile = "/scratch/aob2x/daphnia_hwe_sims/trioPhase/rabbitOut/Scaffold_2373_HRSCAF_2879.AxB_R4_P17_B.out"
magicReconstruct[inputfile, model, popScheme, isFounderInbred -> True, outputFileID -> resultFile, reconstructAlgorithm -> estfun, isPrintTimeElapsed -> True]
summaryFile = StringDrop[resultFile, -4] <> ".csv"
saveAsSummaryMR[resultFile<>"_magicReconstruct.txt", summaryFile]
Exit
