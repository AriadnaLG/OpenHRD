#OpenHRD runing script##

#### Run one function from loading the .CEL files until obtaining the HRDscore) ####

if (dir.exists("/data/RScripts/")) {
  script.dir = "/data/RScripts/"
}

if ( !exists("script.dir") ||  is.na(script.dir) || script.dir == '' || !file.exists( paste0(script.dir, "HRDAnalysis.R") )) {
  print("script folder not set or found")
  quit(status=1)
}


run.HRD.analysis <- function (AT.Cel = NULL, GC.Cel = NUFL, FileName = "",
                              Gender="XX", Platform= "AffyOncoScan", GammaVal= 0.55,
                              Removechr=FALSE, Reference="hg19", AberrantCF=NA, AvgCN=NA
) 
{
  source(paste0(script.dir,"HRDAnalysis.R"))
  
  parallel::detectCores()
  n.cores <- parallel::detectCores() - 1
  my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
  doParallel::registerDoParallel(cl = my.cluster)
  
  cat(paste0("Using ", n.cores, " CPU cores\n"))
  
  LogandBAF <- OS.Process(ATChannelCel = AT.Cel, GCChannelCel = GC.Cel, samplename = FileName)
  print("BAF and L2R ready")
  
  qc<-qcdata(FileName = FileName)
  print(qc)
  
  data.load<- ascat.loadData(Tumor_LogR_file = LogandBAF$TumorLOG, Tumor_BAF_file =LogandBAF$TumorBAF , gender = Gender)
  print("loading ready")
  
  data.gg<- ascat.predictGermlineGenotypes(ASCATobj = data.load, platform ="AffyOncoScan") #is needed to run this function and save the data otherwise ascat.gg doesn't work
  print("gg ready")
  
  Chrominfo<- chrominfoHRD(reference = Reference, chr.in.names=TRUE)
  
  HRDp<-foreach (Penalty = 70, .combine=rbind, .export = c("script.dir")) %dopar% 
    {
      source(paste0(script.dir,"HRDAnalysis.R"))
      
      data.segments<- ascat.aspcf(ASCATobj = data.load, ascat.gg = data.gg, penalty = Penalty, out.prefix = paste0(FileName))
      #print("segments ready")
      
      final.ascat<-ascat.runAscat(ASCATobj = data.segments, filename = paste0(FileName), gamma = GammaVal , rho_manual = AberrantCF, psi_manual = AvgCN,
                                  writetable = F, img.prefix = paste0(FileName))
      #print("ASCAT ready")
      
      table.scarHRD<-tmp.ascat.scarHRD(dataRDS= final.ascat, name=paste0(FileName), writetable=F, removechr=F)
      #print("ready to do HRDscore")
      
      table.TD <- tmp.TD(data = table.scarHRD, Name = paste0(FileName))
      
      HRD.score<-scar_score(seg = table.scarHRD, chr.in.names=TRUE, m, seqz=FALSE, ploidy=NULL, sizelimitLOH=15e6, outputdir=NULL, chrominfo = Chrominfo)
      #print("HRD ready")
      
      scores<- final.scores.openHRD(hrd = HRD.score, qc = qc, Name = paste0(FileName))
      
    }
  
  parallel::stopCluster(cl = my.cluster)
  return(scores)
  
} 