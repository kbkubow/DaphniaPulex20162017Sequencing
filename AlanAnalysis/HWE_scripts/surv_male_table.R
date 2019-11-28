### libraries
  library(data.table)
  library(foreach)

### data
  surv <- fread("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanAnalysis/DataFiles/D82017ClonesMortality.csv")
  male <- fread("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanAnalysis/DataFiles/males_Dorset_forAlan_withClones.csv")

### merge
  surv[,Clone:=gsub(" ", "_", Clone)]
  male[,Clone:=paste(tstrsplit(clone, "_")[[3]], tstrsplit(clone, "_")[[4]], sep="_")]

  setkey(surv, Clone)
  setkey(male, Clone)

  o <- merge(surv, male)

### aggregate
  o.ag <- o[,list(frac=mean(total_males/total_individuals), anyMale=any(total_males>0),
                  alive=mean(AliveMarch72018=="Y"), days=mean(DaystoMortality, na.rm=T)),
              list(Clone, pond)]

  chisq.test(table(male=o.ag$anyMale, alive=o.ag$alive))

  summary(lm(frac~alive + pond, o.ag))
