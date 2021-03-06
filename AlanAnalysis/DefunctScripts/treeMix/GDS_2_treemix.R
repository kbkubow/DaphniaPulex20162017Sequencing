### Run Treemix on Daphnia samples
#module load intel/18.0  intelmpi/18.0 gsl/2.4 R/3.6.3 boost/1.60.0; R
 # module load

### libraries
	library(data.table)
	library(foreach)
	library(SeqArray)

setwd("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/")

### open GDS file
	genofile <- seqOpen("MapJune2020_ann.seq.gds", readonly=TRUE)

### snpFilter file
  #snpFilter <- fread("snpsvarpulexpresentinhalf_table_20200623")
	load("dpfiltsnps_20200623.Rdata")

### load metadata file
  samps <- fread("Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")
  samps[,set:=paste(Species, population, sep="_")]

### functions
	### write treemix data
		writeTreeMix <- function(samps, filestem) {
			o.dat <- foreach(s=na.omit(unique(samps$set)), .combine="cbind")%do%{
        #s <- unique(samps$set)[1]
				print(s)
        subSamp <- samps[set==s,list(clone=clone[which.max(medrd)]), list(SC)]$clone


				seqSetFilter(genofile,
							sample.id=samps[set==s,list(clone=clone[which.max(medrd)]), list(SC)]$clone,
              variant.id=dpfiltsnps$variant.ids)


				dosage <- seqGetData(genofile, "$dosage")

				refCount <- apply(dosage, 2, FUN=function(x) 2*sum(x==0L, na.rm=TRUE) + sum(x==1L, na.rm=TRUE))
				altCount <- apply(dosage, 2, FUN=function(x) 2*sum(x==2L, na.rm=TRUE) + sum(x==1L, na.rm=TRUE))


				o.temp <- data.table(V1=paste(refCount, altCount, sep=","))

				setnames(o.temp, "V1", s)

				o.temp
			}

			seqResetFilter(genofile)

  		### write file
  			gz1 <- gzfile(paste("/scratch/aob2x/treemixIn.", filestem, ".gz", sep=""), "w")
  			write.table(o.dat, gz1, sep=" ", quote=F, row.names=F)
  			close(gz1)
	}

### What is relationship between species/populations

	### write data
		writeTreeMix(samps=samps, filestem="pond_species")

	### run threepop
		#ret <- system("/home/aob2x/treemix/src/threepop -i /scratch/aob2x/treemixIn.pond_species.gz -k 500", intern=T)



  #./treemix \
  #-i /scratch/aob2x/treemixIn.pond_species.gz \
  #-k 500 \
  #-o /scratch/aob2x/treemixIn.pond_species \
  #-root pulicaria_Pond21



### plot

  source("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanAnalysis/treeMix/plotting_funcs.R")

  plot_tree(stem="~/scratch/aob2x/treemixOut.pond_species")
  plot_resid(stem="~/scratch/aob2x/treemixOut.pond_species")




			ret <- system("threepop -i /mnt/spicy_2/treemix_daps/A_B_DCat.gz -k 5000", intern=T)
			ret3 <- system("threepop -i /mnt/spicy_2/treemix_daps/A_B_DBunk_D10.gz -k 500", intern=T)



	### f4 test
		### load metadata file
			samps <- fread("/mnt/spicy_3/Karen/Sp2017/filteredsnps/Superclones20181001.csv")

		### define sets
			samps[Superclone=="A", set:="A"]
			samps[Superclone=="B", set:="B"]
			samps[pond=="DBunk", set:="DBunk"]
			samps[pond=="D10",  set:="D10"]
			samps[grepl("Dcat", pond), set:="W"]

		### write data
			writeTreeMix(samps=samps, filestem="A_B_DBunk_D10_W")

		### run threepop & fourpop
			ret3 <- system("~threepop -i /mnt/spicy_2/treemix_daps/A_B_DBunk_D10_W.gz -k 500", intern=T)
			ret4 <- system("fourpop -i /mnt/spicy_2/treemix_daps/A_B_DBunk_D10_W.gz -k 500", intern=T)
