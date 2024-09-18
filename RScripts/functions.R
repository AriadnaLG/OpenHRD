#### Call packages ####

require(foreach)
require(dplyr)
require(tibble)
require(ggplot2)
require(tidyr)
require(doParallel)
        
#### Main functions ####

## For normalization (EaCoN)  ####
OS.Process <- function(ATChannelCel = NULL, GCChannelCel = NULL, samplename = NULL, dual.norm = TRUE,
                       l2r.level = "weighted", gc.renorm = FALSE, gc.rda = NULL, wave.renorm = FALSE,
                       wave.rda = NULL, mingap = 1E+06, out.dir = getwd(), oschp.keep = FALSE, force.OS = NULL,
                       apt.version = "2.4.0", apt.build = "na33.r2", genome.pkg = "BSgenome.Hsapiens.UCSC.hg19",
                       return.data = FALSE, write.data = TRUE, plot = FALSE, force = FALSE)
  {
  
  ## TEMP
  
  #ATChannelCel = NULL
  #GCChannelCel = NULL
  #samplename = NULL
  #dual.norm = TRUE
  #out.dir = getwd()
  #temp.files.keep = TRUE
  #force.OS = NULL
  #wave.renorm = TRUE
  #wave.rda <- NULL
  #gc.renorm = TRUE
  #gc.rda <- NULL
  #l2r.level <- "normal"
  #apt.build = "na33.r2"
  #apt.version = "2.4.0"
  #mingap = 1E+06
  #genome.pkg = "BSgenome.Hsapiens.UCSC.hg19"
  #BAF.filter = .75
  #force = FALSE
  #write.data = TRUE
  #plot = TRUE
  #return.data = FALSE
  #require(foreach)
  
  `%do%` <- foreach::"%do%"
  
  ## Early checks
  if (is.null(ATChannelCel)) stop(tmsg("An ATChannel CEL file is required !"), call. = FALSE)
  if (is.null(GCChannelCel)) stop(tmsg("A GCChannel CEL file is required !"), call. = FALSE)
  if (is.null(samplename)) stop(tmsg("A samplename is required !"), call. = FALSE)
  if (gc.renorm) { if (!is.null(gc.rda)) { if (!file.exists(gc.rda)) stop(tmsg(paste0("Could not find gc.rda file ", gc.rda)), call. = FALSE) } }
  if (wave.renorm) { if (!is.null(wave.rda)) { if (!file.exists(wave.rda)) stop(tmsg(paste0("Could not find wave.rda file ", wave.rda)), call. = FALSE) } }
  if(!file.exists(ATChannelCel)) stop(tmsg(paste0("Could not find ATChannelCel file ", ATChannelCel, " !")), call. = FALSE)
  if(!file.exists(GCChannelCel)) stop(paste0("Could not find GCChannelCel file ", GCChannelCel, " !"), call. = FALSE)
  if (ATChannelCel == GCChannelCel) stop(tmsg("ATChannelCel and GCChannelCel files are identical !"), call. = FALSE)
  if (!genome.pkg %in% BSgenome::installed.genomes()) {
    if (genome.pkg %in% BSgenome::available.genomes()) {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")), call. = FALSE)
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")), call. = FALSE)
    }
  }
  if (dir.exists(samplename)) { if (!force) stop(tmsg(paste0("A [", samplename, '] dir already exists !')), call. = FALSE) else unlink(samplename, recursive = TRUE, force = FALSE) }
  
  l2r.lev.conv <- list("normal" = "Log2Ratio", "weighted" = "WeightedLog2Ratio")
  if (!(l2r.level %in% names(l2r.lev.conv))) stop(tmsg("Option 'l2r.level' should be 'normal' or 'weighted' !"), call. = FALSE)
  
  ## Handling compressed files
  CEL.OS <- compressed_handler(c(ATChannelCel, GCChannelCel))
  ATChannelCel <- CEL.OS[1]
  GCChannelCel <- CEL.OS[2]
  
  ## Secondary checks
  sup.array <- c("OncoScan", "OncoScan_CNV")
  arraytype.celA = affxparser::readCelHeader(filename = ATChannelCel)$chiptype
  arraytype.celC = affxparser::readCelHeader(filename = GCChannelCel)$chiptype
  if (!arraytype.celA %in% sup.array) stop(tmsg(paste0("Identified array type for ATChannelCel file '", arraytype.celA, "' is not supported by this function !")), call. = FALSE)
  if (!arraytype.celC %in% sup.array) stop(tmsg(paste0("Identified array type for GCChannelCel file '", arraytype.celC, "' is not supported by this function !")), call. = FALSE)
  if (arraytype.celA != arraytype.celC) stop(tmsg(paste0("ATChannelCel and GCChannelCel files are not of the same design : ", arraytype.celA, " != ", arraytype.celC, " !")), call. = FALSE)
  
  ## Checking APT version compatibility
  valid.apt.versions <- c("2.4.0")
  if (!(apt.version %in% valid.apt.versions)) tmsg(paste0("APT version ", apt.version, " is not supported. Program may fail !"))
  
  ## Checking build compatibility
  valid.builds <- c("na33.r2", "na33.r4", "na36.r1")
  if (!(tolower(apt.build) %in% valid.builds)) tmsg(paste0("Build ", apt.build, " is not supported. Program may fail !"))
  
  ## Checking apt-copynumber-onco-ssa package loc
  apt.onco.pkg.name <- paste0("apt.oncoscan.", apt.version)
  if (!(apt.onco.pkg.name %in% installed.packages())) stop(tmsg(paste0("Package ", apt.onco.pkg.name, " not found !")), call. = FALSE)
  suppressPackageStartupMessages(require(package = apt.onco.pkg.name, character.only = TRUE))
  
  ## Processing CELs to an OSCHP file
  if (dir.exists(samplename) && force) unlink(samplename, recursive = TRUE, force = FALSE)
  oscf <- apt.oncoscan.process(ATChannelCel = ATChannelCel, GCChannelCel = GCChannelCel, samplename = samplename, dual.norm = dual.norm, out.dir = out.dir, temp.files.keep = FALSE, force.OS = force.OS, apt.build = apt.build)
  
  ## Reading OSCHP
  my.oschp <- oschp.load(file = oscf)
  if (length(names(my.oschp$Meta)) == 3) names(my.oschp$Meta) <- c("ATChannelCel", "GCChannelCel", "analysis")
  sex.chr <- c("chrX", "chrY")
  
  ## Processing : meta (and checks)
  if (!("affymetrix-chipsummary-snp-qc" %in% names(my.oschp$Meta$analysis))) my.oschp$Meta$analysis[["affymetrix-chipsummary-snp-qc"]] <- NA
  
  ### Loading genome info
  tmsg(paste0("Loading ", genome.pkg, " ..."))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  genome2 <- BSgenome::providerVersion(BSg.obj)
  cs <- chromobjector(BSg.obj)
  
  ### Getting basic meta
  genome <- getmeta("affymetrix-algorithm-param-genome-version", my.oschp$Meta$analysis)
  if (genome != genome2) stop(tmsg(paste0("Genome build name given with BSgenome package '", genome.pkg, "', (", genome2, ") is different from the genome build specified by provided APT build version '", apt.build, "' (", genome, ") !")), call. = FALSE)
  arraytype <- getmeta("affymetrix-array-type", my.oschp$Meta$analysis)
  manufacturer <- getmeta("program-company", my.oschp$Meta$analysis)
  species <- getmeta("affymetrix-algorithm-param-genome-species", my.oschp$Meta$analysis)
  
  gender.conv <- list("female" = "XX", "male" = "XY", "NA" = "NA")
  pgender <- gender.conv[[(getmeta("affymetrix-chipsummary-Y-gender-call", my.oschp$Meta$analysis))]]
  
  if (!(arraytype %in% sup.array)) stop(tmsg(paste0("Unsupported array : '", arraytype, "' !")), call. = FALSE)
  
  meta.b <- list(
    samplename = samplename,
    source = "microarray",
    source.file = list(ATChannelCel = ATChannelCel, GCChannelCel = GCChannelCel),
    type = arraytype,
    manufacturer = manufacturer,
    species = species,
    genome = genome,
    genome.pkg = genome.pkg,
    predicted.gender = pgender
  )
  
  ## Extracting data : L2R
  ao.df <- data.frame(ProbeSetName = my.oschp$ProbeSets$CopyNumber$ProbeSetName, chrs = as.vector(my.oschp$ProbeSets$CopyNumber$Chromosome), pos = as.vector(my.oschp$ProbeSets$CopyNumber$Position), L2R.ori = as.vector(my.oschp$ProbeSets$CopyNumber[[l2r.lev.conv[[l2r.level]]]]), L2R = as.vector(my.oschp$ProbeSets$CopyNumber[[l2r.lev.conv[[l2r.level]]]]), BAF = as.vector(my.oschp$ProbeSets$AllelicData$BAF), AD = as.vector(my.oschp$ProbeSets$AllelicData$AllelicDifference), CallF = as.vector(my.oschp$Genotyping$Calls$ForcedCall), ASignal = my.oschp$Genotyping$Calls$ASignal, BSignal = my.oschp$Genotyping$Calls$BSignal, stringsAsFactors = FALSE)
  affy.chrom <- my.oschp$Chromosomes$Summary
  ak <- affy.chrom$Display
  names(ak) <- affy.chrom$Chromosome
  ao.df$chrA <- as.vector(ak[as.character(ao.df$chrs)])
  ao.df$chr <- paste0("chr", ao.df$chrA)
  ao.df$chrN <- unlist(cs$chrom2chr[ao.df$chr])
  ao.df <- ao.df[order(ao.df$chrN, ao.df$pos, ao.df$ProbeSetName),]
  ao.df <- ao.df[!(is.na(ao.df$L2R) & is.na(ao.df$BAF)),]
  
  ## Adding specific data for FACETS
  rcmat <- round(cbind(ao.df$BAF * (ao.df$ASignal + ao.df$BSignal), (1-ao.df$BAF)*(ao.df$ASignal + ao.df$BSignal)))
  ao.df$LOR <- log(rcmat[,1]+1/6) - log(rcmat[,2]+1/6)
  ao.df$LORvar <- 1/(rcmat[,1]+1/6) + 1/(rcmat[,2]+1/6)
  rm(rcmat)
  gc()
  
  ## L2R renormalizations
  smo <- round(nrow(ao.df) / 550)
  if(smo%%2 == 0) smo <- smo+1
  
  ### Wave  
  if (wave.renorm) {
    tmsg("Wave re-normalization ...")
    
    ren.res <- renorm.go(input.data = ao.df, renorm.rda = wave.rda, track.type = "Wave", smo = smo, arraytype = arraytype, genome = genome)
    
    ao.df <- ren.res$data
    fitted.l2r <- ren.res$renorm$l2r$l2r
    
    if(is.null(ren.res$renorm$pos)) {
      meta.b <- setmeta("wave.renorm", "None", meta.b)
      tmsg(" No positive fit.")
    } else {
      ## Tweaking sex chromosomes
      sex.idx <- ao.df$chr %in% sex.chr
      auto.ori.med <- median(ao.df$L2R[!sex.idx], na.rm = TRUE)
      auto.rn.med <- median(fitted.l2r[!sex.idx], na.rm = TRUE)
      if (any(sex.idx)) {
        for (k in sex.chr) {
          k.idx <- ao.df$chr == k
          if (any(k.idx)) {
            k.ori.diffmed <- median(ao.df$L2R.ori[k.idx], na.rm = TRUE) - auto.ori.med
            k.rn.diffmed <- median(fitted.l2r[k.idx], na.rm = TRUE) - auto.rn.med
            fitted.l2r[k.idx] <- fitted.l2r[k.idx] - k.rn.diffmed + k.ori.diffmed
          }
        }
      }
      meta.b <- setmeta("wave.renorm", paste0(ren.res$mrenorm$pos, collapse = ","), meta.b)
    }
    rm(ren.res)
    ao.df[["L2R.Wave"]] <- fitted.l2r - median(fitted.l2r, na.rm = TRUE)
    ao.df$L2R <- ao.df[["L2R.Wave"]]
  } else {
    meta.b <- setmeta("wave.renorm", "FALSE", meta.b)
  }
  
  ### GC
  if (gc.renorm) {
    tmsg("GC renormalization ...")
    
    ren.res <- renorm.go(input.data = ao.df, renorm.rda = gc.rda, track.type = "GC", smo = smo, arraytype = arraytype, genome = genome)
    ao.df <- ren.res$data
    fitted.l2r <- ren.res$renorm$l2r$l2r
    
    if(is.null(ren.res$renorm$pos)) {
      meta.b <- setmeta("gc.renorm", "None", meta.b)
      message(tmsg(" No positive fit."))
    } else {
      ## Tweaking sex chromosomes
      sex.idx <- ao.df$chr %in% sex.chr
      auto.ori.med <- median(ao.df$L2R[!sex.idx], na.rm = TRUE)
      auto.rn.med <- median(fitted.l2r[!sex.idx], na.rm = TRUE)
      if (any(sex.idx)) {
        for (k in sex.chr) {
          k.idx <- ao.df$chr == k
          if (any(k.idx)) {
            k.ori.diffmed <- median(ao.df$L2R.ori[k.idx], na.rm = TRUE) - auto.ori.med
            k.rn.diffmed <- median(fitted.l2r[k.idx], na.rm = TRUE) - auto.rn.med
            fitted.l2r[k.idx] <- fitted.l2r[k.idx] - k.rn.diffmed + k.ori.diffmed
          }
        }
      }
      meta.b <- setmeta("gc.renorm", paste0(ren.res$renorm$pos, collapse = ","), meta.b)
    }
    rm(ren.res)
    ao.df[["L2R.GC"]] <- fitted.l2r - median(fitted.l2r, na.rm = TRUE)
    ao.df$L2R <- ao.df[["L2R.GC"]]
  } else {
    meta.b <- setmeta("gc.renorm", "FALSE", meta.b)
  }
  
  ## Rough median-centering of L2R
  ao.df$L2R <- ao.df$L2R - median(ao.df$L2R, na.rm = TRUE)
  ao.df$AD[is.nan(ao.df$AD)] <- NA
  
  ## Identifying gaps and clustering chromosomal portions
  gaps <- which(diff(ao.df$pos) >= mingap)
  kends <- vapply(unique(ao.df$chrs), function(k) { max(which(ao.df$chrs == k)) }, 1L)
  kbreaks <- sort(unique(c(gaps, kends)))
  ao.df$chrgap <- rep(seq_along(kbreaks), times = c(kbreaks[1], diff(kbreaks)))
  
  ## Preparing germline
  germ <- ao.df$CallF
  
  
  if(as.numeric(version$major) >= 4) {
    germ[germ %in% as.raw(c(8,11))] <- as.raw(0)
    germ[germ != 0 ] <- as.raw(1)
  } else {
    germ[germ %in% c(8,11)] <- 0
    germ[germ != 0 ] <- 1
  }
  
  ## Building ASCAT-like object
  tmsg("Building normalized object ...")
  my.ch <- sapply(unique(ao.df$chrs), function(x) { which(ao.df$chrs == x) })
  my.ascat.obj <- list(
    data = list(
      Tumor_LogR.ori = data.frame(sample = ao.df$L2R.ori, row.names = ao.df$ProbeSetName),
      Tumor_LogR = data.frame(sample = ao.df$L2R, row.names = ao.df$ProbeSetName),
      Tumor_BAF = data.frame(sample = ao.df$BAF, row.names = ao.df$ProbeSetName),
      Tumor_AD = data.frame(sample = ao.df$AD, row.names = ao.df$ProbeSetName),
      #Tumor_LogR_segmented = NULL,
      #Tumor_BAF_segmented = NULL,
      #Germline_LogR = NULL,
      #Germline_BAF = NULL,
      # SNPpos = data.frame(chrs = ao.df$chrA, pos = ao.df$pos, row.names = rownames(ao.df)),
      SNPpos = data.frame(chrs = ao.df$chr, pos = ao.df$pos, row.names = ao.df$ProbeSetName),
      ch = sapply(unique(ao.df$chr), function(x) { which(ao.df$chr == x) }),
      chr = sapply(unique(ao.df$chrgap), function(x) { which(ao.df$chrgap == x) }),
      chrs = unique(ao.df$chr),
      samples = samplename,
      gender = as.vector(meta.b$predicted.gender),
      sexchromosomes = sex.chr,
      failedarrays = NULL,
      additional = data.frame(RD.test = ao.df$ASignal, RD.ref = ao.df$BSignal, LOR = ao.df$AD, LORvar = ao.df$LORvar, stringsAsFactors = FALSE)
    ),
    meta = list(
      basic = meta.b,
      affy = my.oschp$Meta
    ),
    germline = list(
      germlinegenotypes = matrix(as.logical(germ), ncol = 1),
      failedarrays = NULL
    ),
    CEL = list(
      ATChannelCel = affxparser::readCel(filename = ATChannelCel),
      GCChannelCel = affxparser::readCel(filename = GCChannelCel)
    )
  )
  colnames(my.ascat.obj$germline$germlinegenotypes) <- colnames(my.ascat.obj$data$Tumor_LogR) <- colnames(my.ascat.obj$data$Tumor_LogR.ori) <- colnames(my.ascat.obj$data$Tumor_BAF) <- samplename
  rownames(my.ascat.obj$germline$germlinegenotypes) <- ao.df$ProbeSetName
  my.ascat.obj$CEL$ATChannelCel$intensities <- as.integer(my.ascat.obj$CEL$ATChannelCel$intensities)
  my.ascat.obj$CEL$GCChannelCel$intensities <- as.integer(my.ascat.obj$CEL$GCChannelCel$intensities)
  gc()
  
  
  #To obtain the files for next steps
  #require(dplyr)
  #require(tibble)
  
  LogRdata<- bind_cols(my.ascat.obj$data$SNPpos, my.ascat.obj$data$Tumor_LogR)
  #print(head(LogRdata))
  LogRdata.ori<- bind_cols(my.ascat.obj$data$SNPpos, my.ascat.obj$data$Tumor_LogR.ori)
  
  
  BAFdata <- bind_cols(my.ascat.obj$data$SNPpos, my.ascat.obj$data$Tumor_BAF)
  #print(head(LogBAF))
  
  #if (write.data)
  #  write.table(LogRdata, paste0(out.dir, "/", samplename, "/", samplename, "_LogR.txt"), sep = "\t", row.names = T)
  #  write.table(BAFdata, paste0(out.dir, "/", samplename, "/", samplename, "_BAF.txt"), sep = "\t", row.names = T)
  
  
  
  genopos <- ao.df$pos + cs$chromosomes$chr.length.toadd[ao.df$chrN]
  rm(ao.df)
  gc()
  if (plot) {
    tmsg("Plotting ...")
    kend <- genopos[vapply(my.ascat.obj$data$ch, max, 1L)]
    l2r.notna <- which(!is.na(my.ascat.obj$data$Tumor_LogR[,1]))
    l2r.rm <- runmed(my.ascat.obj$data$Tumor_LogR[,1][l2r.notna], smo)
    l2r.ori.rm <- runmed(my.ascat.obj$data$Tumor_LogR.ori[,1][l2r.notna], smo)
    png(paste0(out.dir, "/", samplename, "/", samplename, "_", arraytype, "_", genome, "_rawplot.png"), 1600, 1050)
    par(mfrow = c(3,1))
    plot(genopos, my.ascat.obj$data$Tumor_LogR.ori[,1], pch = ".", cex = 3, col = "grey70", xaxs = "i", yaxs = "i", ylim = c(-2,2), main = paste0(samplename, " ", arraytype, " raw L2R profile (median-centered) / ", round(sum(abs(diff(l2r.ori.rm))), digits = 3)), xlab = "Genomic position", ylab = "L2R")
    lines(genopos[l2r.notna], l2r.ori.rm, col = 1)
    abline(v = kend, col = 4, lty = 3, lwd = 2)
    abline(h = 0, col = 2, lty = 2, lwd = 2)
    plot(genopos, my.ascat.obj$data$Tumor_LogR[,1], pch = ".", cex = 3, col = "grey70", xaxs = "i", yaxs = "i", ylim = c(-2,2), main = paste0(samplename, " ", arraytype, " L2R profile (", l2r.level, ", median-centered)) / ", round(sum(abs(diff(l2r.rm))), digits = 3)), xlab = "Genomic position", ylab = "L2R")
    lines(genopos[l2r.notna], l2r.rm, col = 1)
    abline(v = kend, col = 4, lty = 3, lwd = 2)
    abline(h = 0, col = 2, lty = 2, lwd = 2)
    plot(genopos, my.ascat.obj$data$Tumor_BAF[,1], pch = ".", cex = 3, col = "grey75", xaxs = "i", yaxs = "i", ylim = c(0,1), main = paste0(samplename, " ", arraytype, " BAF profile"), xlab = "Genomic position", ylab = "BAF")
    points(genopos[germ == 0], my.ascat.obj$data$Tumor_BAF[germ == 0,1], pch = ".", cex = 3, col = 2)
    # plot(ao.df$genopos, my.ascat.obj$data$Tumor_BAF[,1], pch = ".", cex = 3, col = ao.df$ForcedCall-5, xaxs = "i", yaxs = "i", ylim = c(0,1), main = paste0(samplename, " ", arraytype, " BAF profile"), xlab = "Genomic position", ylab = "BAF")
    abline(v = kend, col = 4, lty = 3, lwd = 2)
    abline(h = .5, col = 2, lty = 2, lwd = 2)
    dev.off()
  }
  #New position return data 
  LogandBAF<- return(list(TumorLOG= LogRdata, TumorBAF=BAFdata))
  ## Cleaning
  if(!oschp.keep) {
    tmsg("Removing temporary OSCHP file ...")
    file.remove(oscf)
  }
  
  message(tmsg("Done."))
  gc()
  if(return.data) return(my.ascat.obj)
}

## Quality control ####
qcdata<- function(FileName, apt.version="2.4.0", apt.build="na33.r2")
  {
  qc<-read.table(file = paste0(getwd(),"/", FileName, "/",FileName,"_",apt.version, "_",apt.build,".qc.txt"), header = T)
  qc1<- qc[,c("cel_files", "Y.Gender", "MAPD", "ndSNPQC")]
  
  write.table(x = qc1, file=paste0("qc.",FileName, ".txt"),sep = "\t", row.names = F, col.names = T, quote = T) 
  return(qc1)
  
}

## ASCAT ####
ascat.loadData = function(Tumor_LogR_file= NULL, Tumor_BAF_file=NULL, Germline_LogR_file = NULL, Germline_BAF_file = NULL,
                                                    chrs = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14", #chrs has been change to include "chr"
                                   "chr15","chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22","chrX","chrY"),
                          gender = "XX", sexchromosomes = c("X","Y"))
  {
  
  Tumor_LogR <- Tumor_LogR_file
  Tumor_BAF  <- Tumor_BAF_file
  
  
  #infinite values are a problem - change those
  Tumor_LogR[Tumor_LogR==-Inf]=NA
  Tumor_LogR[Tumor_LogR==Inf]=NA
  
  Germline_LogR = NULL
  Germline_BAF = NULL
  if(!is.null(Germline_LogR_file)) {
    print.noquote("Reading Germline LogR data...")
    Germline_LogR <- read.table(Germline_LogR_file, header=T, row.names=1, comment.char="", sep = "\t", check.names=F)
    print.noquote("Reading Germline BAF data...")
    Germline_BAF <- read.table(Germline_BAF_file, header=T, row.names=1, comment.char="", sep = "\t", check.names=F)
    
    
    #infinite values are a problem - change those
    Germline_LogR[Germline_LogR==-Inf]=NA
    Germline_LogR[Germline_LogR==Inf]=NA
  }
  
  # make SNPpos vector that contains genomic position for all SNPs and remove all data not on chromosome 1-22,X,Y (or whatever is given in the input value of chrs)
  print.noquote("Registering SNP locations...")
  SNPpos <- Tumor_LogR[,1:2]
  SNPpos = SNPpos[SNPpos[,1]%in%chrs,]
  print(head(SNPpos))
  
  # if some chromosomes have no data, just remove them
  chrs = intersect(chrs,unique(SNPpos[,1]))
  
  Tumor_LogR = Tumor_LogR[rownames(SNPpos),c(-1,-2),drop=F]
  Tumor_BAF = Tumor_BAF[rownames(SNPpos),c(-1,-2),drop=F]
  # make sure it is all converted to numerical values
  for (cc in 1:dim(Tumor_LogR)[2]) {
    Tumor_LogR[,cc]=as.numeric(as.vector(Tumor_LogR[,cc]))
    Tumor_BAF[,cc]=as.numeric(as.vector(Tumor_BAF[,cc]))
  }
  
  if(!is.null(Germline_LogR_file)) {
    Germline_LogR = Germline_LogR[rownames(SNPpos),c(-1,-2),drop=F]
    Germline_BAF = Germline_BAF[rownames(SNPpos),c(-1,-2),drop=F]
    for (cc in 1:dim(Germline_LogR)[2]) {
      Germline_LogR[,cc]=as.numeric(as.vector(Germline_LogR[,cc]))
      Germline_BAF[,cc]=as.numeric(as.vector(Germline_BAF[,cc]))
    }
  }
  
  # sort all data by genomic position
  last = 0;
  ch = list();
  SNPorder = vector(length=dim(SNPpos)[1])
  for (i in 1:length(chrs)) {
    chrke = SNPpos[SNPpos[,1]==chrs[i],]
    chrpos = chrke[,2]
    names(chrpos) = rownames(chrke)
    chrpos = sort(chrpos)
    ch[[i]] = (last+1):(last+length(chrpos))
    SNPorder[ch[[i]]] = names(chrpos)
    last = last+length(chrpos)
  }
  SNPpos = SNPpos[SNPorder,]
  Tumor_LogR=Tumor_LogR[SNPorder,,drop=F]
  Tumor_BAF=Tumor_BAF[SNPorder,,drop=F]
  
  if(!is.null(Germline_LogR_file)) {
    Germline_LogR = Germline_LogR[SNPorder,,drop=F]
    Germline_BAF = Germline_BAF[SNPorder,,drop=F]
  }
  
  # split the genome into distinct parts to be used for segmentation (e.g. chromosome arms, parts of genome between gaps in array design)
  print.noquote("Splitting genome in distinct chunks...")
  chr = split_genome(SNPpos)
  
  if (is.null(gender)) {
    gender = rep("XX",dim(Tumor_LogR)[2])
  }
  return(list(Tumor_LogR = Tumor_LogR, Tumor_BAF = Tumor_BAF,
              Tumor_LogR_segmented = NULL, Tumor_BAF_segmented = NULL,
              Germline_LogR = Germline_LogR, Germline_BAF = Germline_BAF,
              SNPpos = SNPpos, ch = ch, chr = chr, chrs = chrs,
              samples = colnames(Tumor_LogR), gender = gender,
              sexchromosomes = sexchromosomes,
              failedarrays = NULL))
  
}
#Dont use this one for the analysis
ascat.GCcorrect = function(ASCATobj, GCcontentfile = NULL)
  {
  if(is.null(GCcontentfile)) {
    print.noquote("Error: no GC content file given!")
  }
  else {
    GC_newlist<-read.table(file=GCcontentfile,header=TRUE,as.is=TRUE)
    colnames(GC_newlist)[c(1,2)] = c("Chr","Position")
    GC_newlist$Chr<-as.character(GC_newlist$Chr)
    GC_newlist$Position<-as.numeric(as.character(GC_newlist$Position))
    
    ovl = intersect(row.names(ASCATobj$Tumor_LogR),row.names(GC_newlist))
    
    GC_newlist<-GC_newlist[ovl,]
    
    SNPpos = ASCATobj$SNPpos[ovl,]
    Tumor_LogR = ASCATobj$Tumor_LogR[ovl,,drop=F]
    Tumor_BAF = ASCATobj$Tumor_BAF[ovl,,drop=F]
    
    chrs = intersect(ASCATobj$chrs,unique(SNPpos[,1]))
    
    Germline_LogR = NULL
    Germline_BAF = NULL
    if(!is.null(ASCATobj$Germline_LogR)) {
      Germline_LogR = ASCATobj$Germline_LogR[ovl,,drop=F]
      Germline_BAF = ASCATobj$Germline_BAF[ovl,,drop=F]
    }
    
    last = 0;
    ch = list();
    for (i in 1:length(ASCATobj$chrs)) {
      chrke = SNPpos[SNPpos[,1]==ASCATobj$chrs[i],]
      chrpos = chrke[,2]
      names(chrpos) = rownames(chrke)
      chrpos = sort(chrpos)
      ch[[i]] = (last+1):(last+length(chrpos))
      last = last+length(chrpos)
    }
    
    for (s in 1:length(ASCATobj$samples)) {
      print.noquote(paste("Sample ", ASCATobj$samples[s], " (",s,"/",length(ASCATobj$samples),")",sep=""))
      Tumordata = Tumor_LogR[,s]
      names(Tumordata) = rownames(Tumor_LogR)
      
      # Calculate weighted correlation
      length_tot<-NULL
      corr_tot<-NULL
      for(chrindex in unique(SNPpos[,1])) {
        GC_newlist_chr<-GC_newlist[GC_newlist$Chr==chrindex,]
        td_chr<-Tumordata[GC_newlist$Chr==chrindex]
        
        flag_nona<-(complete.cases(td_chr) & complete.cases(GC_newlist_chr))
        
        #only work with chromosomes that have variance
        chr_var=var(td_chr[flag_nona])#Will be NA if there is exactly one element.
        if(length(td_chr[flag_nona])>0 && !is.na(chr_var) && chr_var>0){
          corr<-cor(GC_newlist_chr[flag_nona,3:ncol(GC_newlist_chr)],td_chr[flag_nona])
          corr_tot<-cbind(corr_tot,corr)
          length_tot<-c(length_tot,length(td_chr))
        }
      }
      corr<-apply(corr_tot,1,function(x) sum(abs(x*length_tot))/sum(length_tot))
      index_1M<-c(which(names(corr)=="X1M"),which(names(corr)=="X1Mb"))
      maxGCcol_short<-which(corr[1:(index_1M-1)]==max(corr[1:(index_1M-1)]))
      maxGCcol_long<-which(corr[index_1M:length(corr)]==max(corr[index_1M:length(corr)]))
      maxGCcol_long<-(maxGCcol_long+(index_1M-1))
      
      cat("weighted correlation: ",paste(names(corr),format(corr,digits=2), ";"),"\n")
      cat("Short window size: ",names(GC_newlist)[maxGCcol_short+2],"\n")
      cat("Long window size: ",names(GC_newlist)[maxGCcol_long+2],"\n")
      
      # Multiple regression
      flag_NA<-(is.na(Tumordata))|(is.na(GC_newlist[,2+maxGCcol_short]))|(is.na(GC_newlist[,2+maxGCcol_long]))
      td_select<-Tumordata[!flag_NA]
      GC_newlist_select <- GC_newlist[!flag_NA,]
      x1<-GC_newlist_select[,2+maxGCcol_short]
      x2<-GC_newlist_select[,2+maxGCcol_long]
      x3<-(x1)^2
      x4<-(x2)^2
      model<-lm(td_select~x1+x2+x3+x4,y=TRUE)
      
      GCcorrected<-Tumordata
      GCcorrected[]<-NA
      GCcorrected[!flag_NA] <- model$residuals
      
      Tumor_LogR[,s] = GCcorrected
      
      chr = split_genome(SNPpos)
    }
    
    # add some plotting code for each sample while it is generated!!!!
    
    return(list(Tumor_LogR = Tumor_LogR, Tumor_BAF = Tumor_BAF,
                Tumor_LogR_segmented = NULL, Tumor_BAF_segmented = NULL,
                Germline_LogR = Germline_LogR, Germline_BAF = Germline_BAF,
                SNPpos = SNPpos, ch = ch, chr = chr, chrs = chrs,
                samples = colnames(Tumor_LogR), gender = ASCATobj$gender,
                sexchromosomes = ASCATobj$sexchromosomes))
  }
} 

ascat.predictGermlineGenotypes = function(ASCATobj, platform ="AffyOncoScan",  img.dir=".", img.prefix="")
  {
  Homozygous = matrix(nrow = dim(ASCATobj$Tumor_LogR)[1], ncol = dim(ASCATobj$Tumor_LogR)[2])
  colnames(Homozygous) = colnames(ASCATobj$Tumor_LogR)
  rownames(Homozygous) = rownames(ASCATobj$Tumor_LogR)
  
  
  if (platform=="Custom10k") {
    maxHomozygous = 0.05
    proportionHetero = 0.59
    proportionHomo = 0.38
    proportionOpen = 0.02
    segmentLength = 20
  }
  else if (platform=="IlluminaASA") {
    maxHomozygous = 0.05
    proportionHetero = 0.15
    proportionHomo = 0.82
    proportionOpen = 0.01
    segmentLength = 100
  }
  else if (platform=="IlluminaGSAv3") {
    maxHomozygous = 0.05
    proportionHetero = 0.16
    proportionHomo = 0.80
    proportionOpen = 0.01
    segmentLength = 100
  }
  else if (platform=="Illumina109k") {
    maxHomozygous = 0.05
    proportionHetero = 0.35
    proportionHomo = 0.60
    proportionOpen = 0.02
    segmentLength = 20
  }
  else if (platform=="IlluminaCytoSNP") {
    maxHomozygous = 0.05
    proportionHetero = 0.28
    proportionHomo = 0.62
    proportionOpen = 0.03
    segmentLength = 100
  }
  else if (platform=="IlluminaCytoSNP850k") {
    maxHomozygous = 0.05
    proportionHetero = 0.23
    proportionHomo = 0.72
    proportionOpen = 0.01
    segmentLength = 60
  }
  else if (platform=="Illumina610k") {
    maxHomozygous = 0.05
    proportionHetero = 0.295
    proportionHomo = 0.67
    proportionOpen = 0.015
    segmentLength = 30
  }
  else if (platform=="Illumina660k") {
    maxHomozygous = 0.05
    proportionHetero = 0.295
    proportionHomo = 0.67
    proportionOpen = 0.015
    segmentLength = 30
  }
  else if (platform=="Illumina700k") {
    maxHomozygous = 0.05
    proportionHetero = 0.295
    proportionHomo = 0.67
    proportionOpen = 0.015
    segmentLength = 30
  }
  else if (platform=="Illumina1M") {
    maxHomozygous = 0.05
    proportionHetero = 0.22
    proportionHomo = 0.74
    proportionOpen = 0.02
    segmentLength = 100
    #previousvalues:
    #proportionHetero = 0.24
    #proportionOpen = 0.01
    #segmentLength = 60
  }
  else if (platform=="Illumina2.5M") {
    maxHomozygous = 0.05
    proportionHetero = 0.21
    proportionHomo = 0.745
    proportionOpen = 0.03
    segmentLength = 100
  }
  else if (platform=="IlluminaOmni5") {
    maxHomozygous = 0.05
    proportionHetero = 0.13
    proportionHomo = 0.855
    proportionOpen = 0.01
    segmentLength = 100
  }
  else if (platform=="Affy10k") {
    maxHomozygous = 0.04
    proportionHetero = 0.355
    proportionHomo = 0.605
    proportionOpen = 0.025
    segmentLength = 20
  }
  else if (platform=="Affy100k") {
    maxHomozygous = 0.05
    proportionHetero = 0.27
    proportionHomo = 0.62
    proportionOpen = 0.09
    segmentLength = 30
  }
  else if (platform=="Affy250k_sty") {
    maxHomozygous = 0.05
    proportionHetero = 0.26
    proportionHomo = 0.66
    proportionOpen = 0.05
    segmentLength = 50
  }
  else if (platform=="Affy250k_nsp") {
    maxHomozygous = 0.05
    proportionHetero = 0.26
    proportionHomo = 0.66
    proportionOpen = 0.05
    segmentLength = 50
  }
  else if (platform=="AffySNP6") {
    maxHomozygous = 0.05
    proportionHetero = 0.25
    proportionHomo = 0.67
    proportionOpen = 0.04
    segmentLength = 100
  }
  else if (platform=="AffyOncoScan") {
    maxHomozygous = 0.04
    proportionHetero = 0.355
    proportionHomo = 0.605
    proportionOpen = 0.025
    segmentLength = 30
    # maxHomozygous = 0.05
    # proportionHetero = 0.24
    # proportionHomo = 0.69
    # proportionOpen = 0.04
    # segmentLength = 30
  }
  else if (platform=="AffyCytoScanHD") {
    # maxHomozygous = 0.05
    # proportionHetero = 0.26
    # proportionHomo = 0.69
    # proportionOpen = 0.03
    # segmentLength = 30
    maxHomozygous = 0.04
    proportionHetero = 0.32
    proportionHomo = 0.60
    proportionOpen = 0.03
    segmentLength = 100
  }
  else {
    print("Error: platform unknown")
  }
  
  failedarrays = NULL
  
  for (i in 1:dim(ASCATobj$Tumor_LogR)[2]) {
    
    Tumor_BAF_noNA = ASCATobj$Tumor_BAF[!is.na(ASCATobj$Tumor_BAF[,i]),i]
    names(Tumor_BAF_noNA) = rownames(ASCATobj$Tumor_BAF)[!is.na(ASCATobj$Tumor_BAF[,i])]
    Tumor_LogR_noNA = ASCATobj$Tumor_LogR[names(Tumor_BAF_noNA),i]
    names(Tumor_LogR_noNA) = names(Tumor_BAF_noNA)
    
    chr_noNA = list()
    prev = 0
    for(j in 1:length(ASCATobj$chr)) {
      chrke = ASCATobj$chr[[j]]
      next2 = prev + sum(!is.na(ASCATobj$Tumor_BAF[chrke,i]))
      chr_noNA[[j]] = (prev+1):next2
      prev = next2
    }
    
    ch_noNA = list()
    prev = 0
    for(j in 1:length(ASCATobj$ch)) {
      chrke = ASCATobj$ch[[j]]
      next2 = prev + sum(!is.na(ASCATobj$Tumor_BAF[chrke,i]))
      ch_noNA[[j]] = (prev+1):next2
      prev = next2
    }
    
    tbsam = Tumor_BAF_noNA
    # sample, mirrored
    bsm = ifelse(tbsam<0.5, tbsam, 1-tbsam)
    
    homoLimit = max(sort(bsm)[round(length(bsm)*proportionHomo)],maxHomozygous)
    
    if(homoLimit>0.25) {
      failedarrays = c(failedarrays,ASCATobj$samples[i])
    }
    
    Hom = ifelse(bsm<homoLimit,T,NA)
    
    Homo = sum(Hom==T, na.rm=T)
    Undecided = sum(is.na(Hom))
    
    extraHetero = round(min(proportionHetero * length(Tumor_BAF_noNA),Undecided-proportionOpen*length(Tumor_BAF_noNA)))
    
    Hetero = 0
    
    if(extraHetero>0) {
      
      allProbes=1:length(Tumor_BAF_noNA)
      nonHomoProbes = allProbes[is.na(Hom)|Hom==F]
      
      lowestDist = NULL
      
      # bsm, with homozygous replaced by NA
      bsmHNA=bsm
      bsmHNA[!is.na(Hom)&Hom]=NA
      
      for (chrke in chr_noNA) {
        
        chrNonHomoProbes = intersect(nonHomoProbes,chrke)
        
        # there must be a minimum number of probes on the chromosome, otherwise these are called homozygous anyway
        if (length(chrNonHomoProbes)>5) {
          
          #make sure we're not going over any borders..
          segmentLength2 = min(length(chrNonHomoProbes)-1,segmentLength)
          
          chrNonHomoProbesStartWindowLeft = c(rep(NA,segmentLength2),chrNonHomoProbes[1:(length(chrNonHomoProbes)-segmentLength2)])
          chrNonHomoProbesEndWindowLeft = c(NA,chrNonHomoProbes[1:(length(chrNonHomoProbes)-1)])
          chrNonHomoProbesStartWindowRight = c(chrNonHomoProbes[2:length(chrNonHomoProbes)],NA)
          chrNonHomoProbesEndWindowRight = c(chrNonHomoProbes[(segmentLength2+1):length(chrNonHomoProbes)],rep(NA,segmentLength2))
          chrNonHomoProbesStartWindowMiddle = c(rep(NA,segmentLength2/2),chrNonHomoProbes[1:(length(chrNonHomoProbes)-segmentLength2/2)])
          chrNonHomoProbesEndWindowMiddle = c(chrNonHomoProbes[(segmentLength2/2+1):length(chrNonHomoProbes)],rep(NA,segmentLength2/2))
          
          chrLowestDist = NULL
          
          for (probeNr in 1:length(chrNonHomoProbes)) {
            probe = chrNonHomoProbes[probeNr]
            if(!is.na(chrNonHomoProbesStartWindowLeft[probeNr])&!is.na(chrNonHomoProbesEndWindowLeft[probeNr])) {
              medianLeft = median(bsmHNA[chrNonHomoProbesStartWindowLeft[probeNr]:chrNonHomoProbesEndWindowLeft[probeNr]], na.rm=T)
            }
            else {
              medianLeft = NA
            }
            if(!is.na(chrNonHomoProbesStartWindowRight[probeNr])&!is.na(chrNonHomoProbesEndWindowRight[probeNr])) {
              medianRight = median(bsmHNA[chrNonHomoProbesStartWindowRight[probeNr]:chrNonHomoProbesEndWindowRight[probeNr]], na.rm=T)
            }
            else {
              medianRight = NA
            }
            
            if(!is.na(chrNonHomoProbesStartWindowMiddle[probeNr])&!is.na(chrNonHomoProbesEndWindowMiddle[probeNr])) {
              medianMiddle = median(c(bsmHNA[chrNonHomoProbesStartWindowMiddle[probeNr]:chrNonHomoProbesEndWindowLeft[probeNr]],
                                      bsmHNA[chrNonHomoProbesStartWindowRight[probeNr]:chrNonHomoProbesEndWindowMiddle[probeNr]]), na.rm=T)
            }
            else {
              medianMiddle = NA
            }
            
            chrLowestDist[probeNr] = min(abs(medianLeft-bsm[probe]),abs(medianRight-bsm[probe]),abs(medianMiddle-bsm[probe]),Inf,na.rm=T)
          }
        }
        
        # if too few probes on the chromosome
        else {
          chrLowestDist = NULL
          if (length(chrNonHomoProbes)>0) {
            # 1 is higher than any practical distance
            chrLowestDist[1:length(chrNonHomoProbes)] = 1
          }
        }
        
        lowestDist = c(lowestDist,chrLowestDist)
      }
      
      lowestDistUndecided = lowestDist[is.na(Hom[nonHomoProbes])]
      names(lowestDistUndecided)=names(Tumor_LogR_noNA)[nonHomoProbes[is.na(Hom[nonHomoProbes])]]
      
      sorted = sort(lowestDistUndecided)
      Hom[names(sorted[1:min(length(sorted),extraHetero)])] = F
      
      Hetero = sum(Hom==F, na.rm=T)
      Homo = sum(Hom==T, na.rm=T)
      Undecided = sum(is.na(Hom))
      
    }
    
    png(filename = file.path(img.dir,paste(img.prefix, "tumorSep",colnames(ASCATobj$Tumor_LogR)[i],".png",sep="")), width = 2000, height = 500, res = 200)
    title = paste(paste(colnames(ASCATobj$Tumor_BAF)[i], Hetero), Homo)
    ascat.plotGenotypes(ASCATobj,title,Tumor_BAF_noNA, Hom, ch_noNA)
    dev.off()
    
    # set all Undecided to homozygous
    Hom[is.na(Hom)] = T
    Homozygous[names(Hom),i] = Hom
  }
  
  return(list(germlinegenotypes = Homozygous, failedarrays = failedarrays))
  
}

ascat.aspcf =  function(ASCATobj, selectsamples = 1:length(ASCATobj$samples), ascat.gg, penalty = 70, #Penalty has changed from 25 to 70 since it produces better results 
                        out.dir= ".", out.prefix="")
  {
  #first, set germline genotypes
  gg = NULL
  if(!is.null(ascat.gg)) {
    gg = ascat.gg$germlinegenotypes
  }
  else {
    gg = ASCATobj$Germline_BAF < 0.3 | ASCATobj$Germline_BAF > 0.7
  }
  #print("gg")
  # calculate germline homozygous stretches for later resegmentation
  ghs = predictGermlineHomozygousStretches(ASCATobj$chr, gg)
  #print("ghs")
  
  segmentlengths = unique(c(penalty,25,50,100,200,400,800))
  segmentlengths = segmentlengths[segmentlengths>=penalty]
  #print("seglengths")
  
  Tumor_LogR_segmented = matrix(nrow = dim(ASCATobj$Tumor_LogR)[1], ncol = dim(ASCATobj$Tumor_LogR)[2])
  rownames(Tumor_LogR_segmented) = rownames(ASCATobj$Tumor_LogR)
  colnames(Tumor_LogR_segmented) = colnames(ASCATobj$Tumor_LogR)
  #print("tlogseg")
  
  Tumor_BAF_segmented = list();
  for (sample in selectsamples) {
    print.noquote(paste("Sample ", ASCATobj$samples[sample], " (",sample,"/",length(ASCATobj$samples),")",sep=""))
    logrfilename = file.path(out.dir, paste(out.prefix, ASCATobj$samples[sample],".LogR.PCFed.txt",sep=""))
    baffilename = file.path(out.dir, paste(out.prefix, ASCATobj$samples[sample],".BAF.PCFed.txt",sep=""))
    logRPCFed = numeric(0)
    bafPCFed = numeric(0)
    for (segmentlength in segmentlengths) {
      logRPCFed = numeric(0)
      bafPCFed = numeric(0)
      tbsam = ASCATobj$Tumor_BAF[,sample]
      names(tbsam) = rownames(ASCATobj$Tumor_BAF)
      homosam = gg[,sample]
      for (chrke in 1:length(ASCATobj$chr)) {
        lr = ASCATobj$Tumor_LogR[ASCATobj$chr[[chrke]],sample]
        #winsorize to remove outliers
        #this has a problem with NAs
        lrwins = vector(mode="numeric",length=length(lr))
        lrwins[is.na(lr)] = NA
        lrwins[!is.na(lr)] = madWins(lr[!is.na(lr)],2.5,25)$ywin
        baf = tbsam[ASCATobj$chr[[chrke]]]
        homo = homosam[ASCATobj$chr[[chrke]]]
        Select_het <- !homo & !is.na(homo) & !is.na(baf) & !is.na(lr)
        bafsel = baf[Select_het]
        # winsorize BAF as well (as a safeguard)
        bafselwinsmirrored = madWins(ifelse(bafsel>0.5,bafsel,1-bafsel),2.5,25)$ywin
        bafselwins = ifelse(bafsel>0.5,bafselwinsmirrored,1-bafselwinsmirrored)
        indices = which(Select_het)
        logRaveraged = NULL;
        if(length(indices)!=0) {
          averageIndices = c(1,(indices[1:(length(indices)-1)]+indices[2:length(indices)])/2,length(lr)+0.01)
          startindices = ceiling(averageIndices[1:(length(averageIndices)-1)])
          endindices = floor(averageIndices[2:length(averageIndices)]-0.01)
          if(length(indices)==1) {
            startindices = 1
            endindices = length(lr)
          }
          nrIndices = endindices - startindices + 1
          logRaveraged = vector(mode="numeric",length=length(indices))
          for(i in 1:length(indices)) {
            if(is.na(endindices[i])) {
              endindices[i]=startindices[i]
            }
            logRaveraged[i]=mean(lrwins[startindices[i]:endindices[i]], na.rm=T)
          }
        }
        # if there are no probes in the segment (after germline homozygous removal), don't do anything, except add a LogR segment
        if(length(logRaveraged)>0) {
          logRASPCF = NULL
          bafASPCF = NULL
          if(length(logRaveraged)<6) {
            logRASPCF = rep(mean(logRaveraged),length(logRaveraged))
            bafASPCF = rep(mean(bafselwins),length(logRaveraged))
          }
          else {
            PCFed = fastAspcf(logRaveraged,bafselwins,6,segmentlength)
            logRASPCF = PCFed$yhat1
            bafASPCF = PCFed$yhat2
          }
          names(bafASPCF)=names(indices)
          logRc = numeric(0)
          for(probe in 1:length(logRASPCF)) {
            if(probe == 1) {
              logRc = rep(logRASPCF[probe],indices[probe])
            }
            # if probe is 1, set the beginning, and let the loop go:
            if(probe == length(logRASPCF)) {
              logRc = c(logRc,rep(logRASPCF[probe],length(lr)-indices[probe]))
            }
            else if(logRASPCF[probe]==logRASPCF[probe+1]) {
              logRc = c(logRc, rep(logRASPCF[probe],indices[probe+1]-indices[probe]))
            }
            else {
              #find best breakpoint
              d = numeric(0)
              totall = indices[probe+1]-indices[probe]
              for (bp in 0:(totall-1)) {
                dis = sum(abs(lr[(1:bp)+indices[probe]]-logRASPCF[probe]), na.rm=T)
                if(bp!=totall) {
                  dis = sum(dis, sum(abs(lr[((bp+1):totall)+indices[probe]]-logRASPCF[probe+1]), na.rm=T), na.rm=T)
                }
                d = c(d,dis)
              }
              breakpoint = which.min(d)-1
              logRc = c(logRc,rep(logRASPCF[probe],breakpoint),rep(logRASPCF[probe+1],totall-breakpoint))
            }
          }
          #2nd step: adapt levels!
          logRd = numeric(0)
          seg = rle(logRc)$lengths
          startprobe = 1
          endprobe = 0
          for (i in 1:length(seg)) {
            endprobe = endprobe+seg[i]
            level = mean(lr[startprobe:endprobe], na.rm=T)
            logRd = c(logRd, rep(level,seg[i]))
            startprobe = startprobe + seg[i]
          }
          logRPCFed = c(logRPCFed,logRd)
          bafPCFed = c(bafPCFed,bafASPCF)
        }
        # add a LogR segment
        else {
          level = mean(lr,na.rm=T)
          reps = length(lr)
          logRPCFed = c(logRPCFed,rep(level,reps))
        }
        # correct wrong segments in germline homozygous stretches:
        homsegs = ghs[[sample]][ghs[[sample]][,1]==chrke,]
        startchr = min(ASCATobj$chr[[chrke]])
        endchr = max(ASCATobj$chr[[chrke]])
        # to solve an annoying error when homsegs has length 1:
        if(length(homsegs)==3) {
          homsegs=t(as.matrix(homsegs))
        }
        if(!is.null(homsegs)&&!is.na(homsegs)&&dim(homsegs)[1]!=0) {
          for (i in 1:dim(homsegs)[1]) {
            # note that only the germline homozygous segment is resegmented, plus a bit extra (but this is NOT replaced)
            startpos = max(homsegs[i,2],startchr)
            endpos = min(homsegs[i,3],endchr)
            # PCF over a larger fragment
            startpos2 = max(homsegs[i,2]-100,startchr)
            endpos2 = min(homsegs[i,3]+100,endchr)
            # take into account a little extra (difference between startpos2 and startpos3 is not changed)
            startpos3 = max(homsegs[i,2]-5,startchr)
            endpos3 = min(homsegs[i,3]+5,endchr)
            # note that the parameters are arbitrary, but <100 seems to work on the ERBB2 example!
            # segmentlength is lower here, as in the full data, noise on LogR is higher!
            # do this on Winsorized data too!
            towins = ASCATobj$Tumor_LogR[startpos2:endpos2,sample]
            winsed = madWins(towins[!is.na(towins)],2.5,25)$ywin
            pcfed = vector(mode="numeric",length=length(towins))
            pcfed[!is.na(towins)] = exactPcf(winsed,6,floor(segmentlength/4))
            pcfed2 = pcfed[(startpos3-startpos2+1):(endpos3-startpos2+1)]
            dif = abs(pcfed2-logRPCFed[startpos3:endpos3])
            #0.3 is hardcoded here, in order not to have too many segments!
            #only replace if enough probes differ (in order not to get singular probes with different breakpoints)
            if(!is.na(dif)&&sum(dif>0.3)>5) {
              #replace a bit more to make sure no 'lone' probes are left (startpos3 instead of startpos)
              logRPCFed[startpos3:endpos3]=ifelse(dif>0.3,pcfed2,logRPCFed[startpos3:endpos3])
            }
          }
        }
      }
      #fill in NAs (otherwise they cause problems):
      #some NA probes are filled in with zero, replace those too:
      logRPCFed = fillNA(logRPCFed, zeroIsNA=TRUE)
      
      #adapt levels again
      seg = rle(logRPCFed)$lengths
      logRPCFed = numeric(0)
      startprobe = 1
      endprobe = 0
      prevlevel = 0
      for (i in 1:length(seg)) {
        endprobe = endprobe+seg[i]
        level = mean(ASCATobj$Tumor_LogR[startprobe:endprobe,sample], na.rm=T)
        #making sure no NA's get filled in...
        if(is.nan(level)) {
          level=prevlevel
        }
        else {
          prevlevel=level
        }
        logRPCFed = c(logRPCFed, rep(level,seg[i]))
        startprobe = startprobe + seg[i]
      }
      #put in names and write results to files
      names(logRPCFed) = rownames(ASCATobj$Tumor_LogR)
      
      # if less than 800 segments: this segmentlength is ok, otherwise, rerun with higher segmentlength
      if(length(unique(logRPCFed))<800) {
        break
      }
    }
    
    write.table(logRPCFed,logrfilename,sep="\t",col.names=F)
    write.table(bafPCFed,baffilename,sep="\t",col.names=F)
    bafPCFed = as.matrix(bafPCFed)
    Tumor_LogR_segmented[,sample] = logRPCFed
    Tumor_BAF_segmented[[sample]] = 1-bafPCFed
  }
  
  ASCATobj = list(Tumor_LogR = ASCATobj$Tumor_LogR,
                  Tumor_BAF = ASCATobj$Tumor_BAF,
                  Tumor_LogR_segmented = Tumor_LogR_segmented,
                  Tumor_BAF_segmented = Tumor_BAF_segmented,
                  Germline_LogR = ASCATobj$Germline_LogR,
                  Germline_BAF = ASCATobj$Germline_BAF,
                  SNPpos = ASCATobj$SNPpos,
                  ch = ASCATobj$ch,
                  chr = ASCATobj$chr,
                  chrs = ASCATobj$chrs,
                  probes = logRPCFed,
                  samples = colnames(ASCATobj$Tumor_LogR), gender = ASCATobj$gender,
                  sexchromosomes = ASCATobj$sexchromosomes, failedarrays = ascat.gg$failedarrays)
  return(ASCATobj)
}

ascat.runAscat = function(ASCATobj, gamma = 0.55, pdfPlot = F, y_limit = 5, circos=NA, rho_manual = NA, psi_manual = NA, # modified to create a RDS file with the CNV analysis data
                          img.dir=".", img.prefix="", writetable=T, filename="")
  {
  goodarrays=NULL
  res = vector("list",dim(ASCATobj$Tumor_LogR)[2])
  for (arraynr in 1:dim(ASCATobj$Tumor_LogR)[2]) {
    print.noquote(paste("Sample ", ASCATobj$samples[arraynr], " (",arraynr,"/",length(ASCATobj$samples),")",sep=""))
    lrr=ASCATobj$Tumor_LogR[,arraynr]
    names(lrr)=rownames(ASCATobj$Tumor_LogR)
    baf=ASCATobj$Tumor_BAF[,arraynr]
    names(baf)=rownames(ASCATobj$Tumor_BAF)
    lrrsegm = ASCATobj$Tumor_LogR_segmented[,arraynr]
    names(lrrsegm) = rownames(ASCATobj$Tumor_LogR_segmented)
    bafsegm = ASCATobj$Tumor_BAF_segmented[[arraynr]][,,drop=FALSE]
    names(bafsegm) = rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]])
    failedqualitycheck = F
    if(ASCATobj$samples[arraynr]%in%ASCATobj$failedarrays) {
      failedqualitycheck = T
    }
    ending = ifelse(pdfPlot, "pdf", "png")
    circosName=NA
    if(!is.na(circos)){
      circosName=paste(circos,"_",ASCATobj$samples[arraynr],sep="")
    }
    if(is.na(rho_manual)) {
      res[[arraynr]] = runASCAT(lrr,baf,lrrsegm,bafsegm,ASCATobj$gender[arraynr],ASCATobj$SNPpos,ASCATobj$ch,ASCATobj$chrs,ASCATobj$sexchromosomes, failedqualitycheck,
                                file.path(img.dir, paste(img.prefix, ASCATobj$samples[arraynr],".sunrise.png",sep="")),file.path(img.dir, paste(img.prefix, ASCATobj$samples[arraynr],".ASCATprofile.", ending ,sep="")),
                                file.path(img.dir, paste(img.prefix, ASCATobj$samples[arraynr],".rawprofile.", ending ,sep="")),NA,
                                gamma,NA,NA,pdfPlot, y_limit, circosName)
    } else {
      res[[arraynr]] = runASCAT(lrr,baf,lrrsegm,bafsegm,ASCATobj$gender[arraynr],ASCATobj$SNPpos,ASCATobj$ch,ASCATobj$chrs,ASCATobj$sexchromosomes, failedqualitycheck,
                                file.path(img.dir, paste(img.prefix, ASCATobj$samples[arraynr],".sunrise.png",sep="")),file.path(img.dir, paste(img.prefix, ASCATobj$samples[arraynr],".ASCATprofile.", ending,sep="")),
                                file.path(img.dir, paste(img.prefix, ASCATobj$samples[arraynr],".rawprofile.", ending,sep="")),NA,
                                gamma,rho_manual[arraynr],psi_manual[arraynr], pdfPlot, y_limit, circosName)
    }
    if(!is.na(res[[arraynr]]$rho)) {
      goodarrays[length(goodarrays)+1] = arraynr
    }
  }
  
  if(length(goodarrays)>0) {
    n1 = matrix(nrow = dim(ASCATobj$Tumor_LogR)[1], ncol = length(goodarrays))
    n2 = matrix(nrow = dim(ASCATobj$Tumor_LogR)[1], ncol = length(goodarrays))
    rownames(n1) = rownames(ASCATobj$Tumor_LogR)
    rownames(n2) = rownames(ASCATobj$Tumor_LogR)
    colnames(n1) = colnames(ASCATobj$Tumor_LogR)[goodarrays]
    colnames(n2) = colnames(ASCATobj$Tumor_LogR)[goodarrays]
    for (i in 1:length(goodarrays)) {
      n1[,i] = res[[goodarrays[i]]]$nA
      n2[,i] = res[[goodarrays[i]]]$nB
    }
    
    distance_matrix = vector("list",length(goodarrays)) 
    names(distance_matrix) <- colnames(ASCATobj$Tumor_LogR)[goodarrays]
    for (i in 1:length(goodarrays)) {
      distance_matrix[[i]] = res[[goodarrays[i]]]$distance_matrix
    }
    
    tp = vector(length=length(goodarrays))
    psi = vector(length=length(goodarrays))
    ploidy = vector(length=length(goodarrays))
    goodnessOfFit = vector(length=length(goodarrays))
    naarrays = NULL
    for (i in 1:length(goodarrays)) {
      tp[i] = res[[goodarrays[i]]]$rho
      psi[i] = res[[goodarrays[i]]]$psi
      ploidy[i] = mean(res[[goodarrays[i]]]$nA+res[[goodarrays[i]]]$nB,na.rm=T)
      goodnessOfFit[i] = res[[goodarrays[i]]]$goodnessOfFit
      if(res[[goodarrays[i]]]$nonaberrant) {
        naarrays = c(naarrays,ASCATobj$samples[goodarrays[i]])
      }
    }
    fa = colnames(ASCATobj$Tumor_LogR)[-goodarrays]
    names(tp) = colnames(n1)
    names(ploidy) = colnames(n1)
    names(psi) = colnames(n1)
    names(goodnessOfFit) = colnames(n1)
    
    Probes= NULL ## Number of probes extracted from the data
    `%do%` <- foreach::"%do%"
    probes.lr2 <- foreach::foreach(k = 1:length(ASCATobj$ch), .combine = "rbind") %do% {
      segrle <- rle(ASCATobj$Tumor_LogR_segmented[ASCATobj$ch[[k]],1])
      probes.log2 <- data.frame(Probes = segrle$lengths, Value = segrle$values, stringsAsFactors = FALSE)
      probes.log2$Chrom <- ASCATobj$chrs[[k]]
      probecs <- cumsum(probes.log2$Probes)
      probes.log2 <- probes.log2[, c(3, 2, 1)]
      
    }
    Probes = data.frame(probes.lr2)
    
    seg = NULL
    for (i in 1:length(goodarrays)) {
      segje = res[[goodarrays[i]]]$seg
      seg = rbind(seg,cbind(ASCATobj$samples[goodarrays[i]],as.vector(ASCATobj$SNPpos[segje[,1],1]),
                            ASCATobj$SNPpos[segje[,1],2],
                            ASCATobj$SNPpos[segje[,2],2],segje[,3],segje[,4]))
    }
    colnames(seg) = c("sample","chr","startpos","endpos","nMajor","nMinor")
    print(str(seg))
    
    seg = data.frame(seg,stringsAsFactors=F)
    seg[,3]=as.numeric(seg[,3])
    seg[,4]=as.numeric(seg[,4])
    seg[,5]=as.numeric(seg[,5])
    seg[,6]=as.numeric(seg[,6])
    
    
    seg_raw = NULL
    for (i in 1:length(goodarrays)) {
      segje = res[[goodarrays[i]]]$seg_raw
      seg_raw = rbind(seg_raw,cbind(ASCATobj$samples[goodarrays[i]], as.vector(ASCATobj$SNPpos[segje[,1],1]),
                                    ASCATobj$SNPpos[segje[,1],2],
                                    ASCATobj$SNPpos[segje[,2],2],segje[,3],segje[,4:ncol(segje)], 
                                    probes.lr2[,3] )) # Number of probes added to the matrix data 
      
    }
    colnames(seg_raw) = c("sample","chr","startpos","endpos","nMajor","nMinor","nAraw","nBraw", "Probes")
    print(str(seg_raw))
    
    seg_raw = data.frame(seg_raw,stringsAsFactors=F)
    seg_raw[,3]=as.numeric(seg_raw[,3])
    seg_raw[,4]=as.numeric(seg_raw[,4])
    seg_raw[,5]=as.numeric(seg_raw[,5])
    seg_raw[,6]=as.numeric(seg_raw[,6])
    seg_raw[,7]=as.numeric(seg_raw[,7])
    seg_raw[,8]=as.numeric(seg_raw[,8])
    seg_raw[,9]=as.numeric(seg_raw[,9])
    
  }
  else {
    n1 = NULL
    n2 = NULL
    tp = NULL
    ploidy = NULL
    psi = NULL
    goodnessOfFit = NULL
    fa = colnames(ASCATobj$Tumor_LogR)
    naarrays = NULL
    seg = NULL
    seg_raw = NULL
    distance_matrix = NULL
    Probes=NULL
    
    
  }
  
  finaldata<- list(nA = n1, nB = n2, aberrantcellfraction = tp, ploidy = ploidy, psi = psi, goodnessOfFit = goodnessOfFit,
                   failedarrays = fa, nonaberrantarrays = naarrays, segments = seg, segments_raw = seg_raw, distance_matrix = distance_matrix)
  if(writetable){
    saveRDS(object = finaldata, file = paste0(filename,"",".RDS"), compress = FALSE )
    
  }
  else{
    return(finaldata)
  }
} 

## Prepare data to obtain genomic scars ####
tmp.ascat.scarHRD<- function(dataRDS= "", name="", writetable=TRUE, removechr=FALSE)
  {

  #data<- readRDS(file = dataRDS)
  data<- dataRDS
  
  ####### Extract required information 
  segdata <- data$segments
  if(removechr){
    segdata$chr <- gsub("chr", "", segdata$chr)
  }
  
  ####### Create Total copy number column 
  segdata <- segdata %>%
    mutate(TCN = nMajor+nMinor)
  
  ####### Select the desired columns: for scar HRD 8 columns are needed
  segmented.data <- as.data.frame(bind_cols(segdata$sample, segdata$chr, segdata$startpos, segdata$endpos, segdata$TCN,
                                            segdata$nMajor, segdata$nMinor, data$ploidy))
  colnames(segmented.data)<-c("sample", "chr", "startpos", "endpos","TCN", "nA", "nB", "ploidy")
  segmented.data <- as.data.frame(segmented.data)
  
  if(writetable){
    write.table(x = segmented.data, file = paste0(name,"",".For.scarHRD.txt"),sep = "\t", row.names = F, col.names = T, quote = T)
  }
  else{
    return(segmented.data)
  }
  
}

tmp.TD<- function(data, Name="")
  {
  tmpl <- data 
  
  tmpl.TD <- tmpl%>%
    mutate( full.loc = paste(startpos, endpos,sep = "-") )%>%
    mutate( full.loc = paste(chr, full.loc ,sep = ":") )%>%
    mutate(Size= endpos - startpos)%>%
    mutate(Status = TCN >2) 
  tmpl.td.status<-as.data.frame(gsub(T, "Gain", tmpl.TD$Status))
  tmpl.td.status <- as.data.frame(gsub(F, "Loss", tmpl.td.status[, 1]))
  
  tmpl.TD <-as.data.frame(bind_cols(tmpl.TD[,1:10], tmpl.td.status))
  colnames(tmpl.TD)<- c("sample", "chr", "startpos", "endpos", "CN State", "nA", "nB", "ploidy", "Full Location", "Size", "Type")
  write.table(x = tmpl.TD, file = paste0(Name, ".Segments.txt" ), col.names=T,sep="\t", row.names=T,  quote=F)
  return(tmpl.TD)
  
  
}

### Obtain the GIS values ####
#To compute TAI, LOH and LST 
scar_score<-function(seg, chr.in.names=TRUE, m,seqz=FALSE, ploidy=NULL, sizelimitLOH=15e6, outputdir=NULL, chrominfo=NULL)
  {
  
  if (seqz==TRUE){
    cat('Preprocessing started...\n')
    seg<-preprocess.seqz(seg, ploidy0=ploidy, chr.in.names=chr.in.names)
    cat('Preprocessing finished \n')
  }
  else {
    # seg<-read.table(seg,header=T, check.names = F, stringsAsFactors = F, sep="\t" )
    seg<-seg
    
    seg[,9]<-rep(1,dim(seg)[1])
    
    
  }
  #prep
  cat('Determining HRD-LOH, LST, TAI \n')
  seg<-preprocess.hrd(seg)
  print(str(seg))
  #Calculating the hrd score:
  res_hrd <- calc.hrd(seg,sizelimit1=15e6) 
  #Calculating the telomeric allelic imbalance score:
  res_ai<- calc.ai_new(seg = seg, chrominfo = chrominfo) 
  #Calculating the large scale transition scores:
  res_lst <- calc.lst(seg = seg, chrominfo = chrominfo) 
  sum_HRD0<-res_lst+res_hrd+res_ai[1]
  
  if (is.null(ploidy)){
    sum_HRDc<-NA
  } else {
    sum_HRDc<-res_lst-15.5*ploidy+res_hrd+res_ai[1]
  }
  
  HRDresulst<-c(res_hrd,res_ai[1],res_lst,sum_HRD0)
  names(HRDresulst)<-c("LOH",colnames(res_ai)[1],"LST", "HRDsum")
  run_name<-names(sum_HRD0)
  #write.table(t(HRDresulst),paste0(outputdir,"/",run_name,"_HRDresults.txt"),col.names=NA,sep="\t",row.names=unique(seg[,1]))
  return(t(HRDresulst))
}

#Return final score values 

final.scores.openHRD<- function(hrd=NULL, qc=qc, Name="")
  {
  #Formula y=xm-b -> b=16.4643 & m=0.8818 -> x=y-b/m
  b=13.4588
  m=0.8968
  
  tmphrd<-as.data.frame(hrd)%>%
    mutate(teorval= round((as.numeric(HRDsum)-b)/m, digits = 0))%>%
    mutate(dif=HRDsum-teorval)%>%
    mutate(HRDcorrected=HRDsum-dif)
  tmphrd$HRDcorrected[tmphrd$HRDcorrected<0]<- 0
  HRD<- cbind(hrd, tmphrd$HRDcorrected)
  
  f.s<- cbind(qc, HRD)
  colnames(f.s)<- c("cel_files", "Y.Gender", "MAPD: <0.30", "ndSNPQC: > 26", "LOH", "TAI", "LST", "HRDsum", "HRDcorrected")
  write.table(f.s,paste0(Name,"_GenomicScars.results.txt"),col.names=T,sep="\t", row.names = F,quote = F )
  return(f.s)
}

################################################################################
#### Sumplementary Functions ####

### Normalization (EaCoN) ####
compressed_handler <- function(CELz = NULL)
  {
  `%do%` <- foreach::"%do%"
  CELz2 <- foreach(CEL = CELz, .combine = "c") %do% {
    tmsg(paste0("Decompressing ", CEL, " ..."))
    if (tolower(tools::file_ext(CEL)) == "bz2") {
      uncomp_file <- tempfile(fileext = ".CEL")
      R.utils::bunzip2(filename = CEL, destname = uncomp_file, FUN = bzfile, remove = FALSE)
      CEL <- uncomp_file
    } else if (tolower(tools::file_ext(CEL)) == "gz") {
      uncomp_file <- tempfile(fileext = ".CEL")
      R.utils::gunzip(filename = CEL, destname = uncomp_file, FUN = gzfile, remove = FALSE)
      CEL <- uncomp_file
    } else if (tolower(tools::file_ext(CEL)) == "zip") {
      zlist <- utils::unzip(CEL, list = TRUE)
      if (length(grep(zlist$Name, pattern = "\\.CEL", ignore.case = TRUE)) != 1) stop(tmsg(paste0(CEL, "archive file does not contain a single and unique CEL file !")), call. = FALSE)
      zname <- zlist$Name[1]
      utils::unzip(zipfile = CEL, files = zname, exdir = tempdir(), overwrite = TRUE)
      CEL <- file.path(tempdir(), zname)
    } else if (tolower(tools::file_ext(CEL)) != "cel") stop(tmsg(paste0("File ", CEL, " is not recognized as raw nor compressed (gz, bz2, zip) CEL file !")), call. = FALSE)
    return(CEL)
  }
  return(CELz2)
}

tmsg <- function(text = NULL) { message(paste0(" [", Sys.info()[['nodename']], ":", Sys.getpid(), "] ", text)) }

oschp.load <- function(file = NULL)
  {
  if (is.null(file)) stop(tmsg("Please provide an OSCHP file !"), call. = FALSE)
  if (!file.exists(file)) stop(tmsg("Provided OSCHP file does not exist !"), call. = FALSE)
  h5.data <- rhdf5::h5read(file = file, name = "/")
  h5.mlist <- h5.data$Dset_IO_HDF5_Gdh
  if(length(h5.mlist) > 1) {
    `%do%` <- foreach::"%do%"
    h5.meta <- foreach::foreach (a = 1:(length(h5.mlist)-1)) %do% {
      h5.meta.c <- foreach (l = 1:(list.depth(h5.mlist[[a]])-1)) %do% {
        tmp.list <- h5.mlist[[a]][["_&keyvals"]]
        h5.mlist[[a]] <- h5.mlist[[a]][[1]]
        return(meta.df2list(tmp.list))
      }
      names(h5.meta.c) <- foreach (l = h5.meta.c, .combine = "c") %do% {
        return(rev(unlist(strsplit(x = l[["data_source"]], split = "-")))[1])
      }
      return(h5.meta.c)
    }
    names(h5.meta) <- paste0("CEL", 1:(length(h5.mlist)-1))
  } else h5.meta <- list()
  h5.meta$analysis = meta.df2list(h5.mlist[["_&keyvals"]])
  h5.data$Meta <- h5.meta
  h5.data$Dset_IO_HDF5_Gdh <- NULL
  return(h5.data)
}

list.depth <- function(this) ifelse(is.list(this), 1L + max(sapply(this, list.depth)), 0L)

meta.df2list <- function(meta.df = NULL)
  {
  return(sapply(seq_len(nrow(meta.df)), function(x) {
    l <- meta.df[x,]
    return(setNames(l[2], l[1]))
  }))
}

#Create a chromosomes-like object from a BSgenome object
chromobjector <- function(BSg = NULL)
  {
  if (is.null(BSg)) stop("NULL object !", call. = FALSE)
  # chromobj <- list(species = GenomeInfoDb::organism(BSg), genomebuild = BSgenome::providerVersion(BSg))
  chromobj <- list(species = GenomeInfoDb::organism(BSg), genomebuild = metadata(BSg)$genome)
  chromdf <- data.frame(chrom = BSgenome::seqnames(BSg), chrN = seq_along(BSgenome::seqnames(BSg)), chr.length = GenomeInfoDb::seqlengths(BSg), stringsAsFactors = FALSE)
  chromdf$chr.length.sum <- cumsum(as.numeric(chromdf$chr.length))
  chromdf$chr.length.toadd <- c(0, chromdf$chr.length.sum[-nrow(chromdf)])
  chromdf$mid.chr <- round(diff(c(0, chromdf$chr.length.sum)) /2)
  chromdf$mid.chr.geno <- chromdf$mid.chr + chromdf$chr.length.toadd
  chromobj$chromosomes <- chromdf
  rm(chromdf)
  chromobj$chrom2chr <- sapply(chromobj$chromosomes$chrom, function(k) { chromobj$chromosomes$chrN[chromobj$chromosomes$chrom == k]}, simplify = FALSE)
  chromobj$chr2chrom <- sapply(chromobj$chromosomes$chrN, function(k) { chromobj$chromosomes$chrom[chromobj$chromosomes$chrN == k]}, simplify = FALSE)
  names(chromobj$chr2chrom) <- chromobj$chromosomes$chrN
  chromobj$genome.length <- sum(as.numeric(chromobj$chromosomes$chr.length), na.rm = TRUE)
  return(chromobj)
}

getmeta <- function(key = NULL, meta = NULL)
  {
  if (key %in% names(meta)) {
    val <- meta[[key]]
    if(is.character(val)) val <- sub(replacement = "   ", x = val, pattern = " // ")
    return(val)
  } else return(NA)
}

setmeta <- function(key = NULL, val = NULL, meta = NULL)
  {
  for (x in 1:length(key)) meta[[key[x]]] <- val[x]
  return(meta)
}

#Main renormalization function
renorm.go <- function(input.data = NULL, renorm.rda = NULL, track.type = "GC", smo = 399, arraytype = NULL, genome = NULL)
  {
  if (!is.null(renorm.rda)) {
    load(renorm.rda, envir = environment())
    rn.arraytype <- renorm.data$info$value[renorm.data$info$key == "array_type"]
    rn.genome <- renorm.data$info$value[renorm.data$info$key == "genome-version"]
    rn.track.type <- renorm.data$info$value[renorm.data$info$key == "track_type"]
    if ((rn.track.type != track.type) | (rn.arraytype != arraytype) | (rn.genome != genome)) stop(tmsg(paste0("Provided renormalization pack is not as intended ! Expected [", track.type, ", ", arraytype, ", ", genome, "], got [", rn.track.type, ", ", rn.arraytype, ", ", rn.genome, "] !")), call. = FALSE)
  } else {
    RN.pkg.name <- "affy.CN.norm.data"
    if (!(RN.pkg.name %in% installed.packages())) stop(tmsg(paste0("Package ", RN.pkg.name, " not found !")), call. = FALSE)
    RN.file <- system.file(paste0("data/", arraytype, ".", genome, ".", track.type, ".rda"), package = RN.pkg.name)
    if (RN.file == "") stop(tmsg(paste0("Could not find a ", track.type, " data package for [", arraytype, ", ", genome, "] in package '", RN.pkg.name, "' ! Please build your own Wave data pack with ", RN.pkg.name, "::affy.wave.compute() and submit it using the 'renorm.rda' option.")), call. = FALSE)
    data(list = paste0(arraytype, ".", genome, ".", track.type), package = RN.pkg.name, envir = environment())
  }
  # print(str(RNdata))
  # print(str(rownames(input.data)))
  RNdata <- renorm.data$tracks[renorm.data$tracks$ProbeSetName %in% input.data$ProbeSetName,]
  input.data <- input.data[input.data$ProbeSetName %in% RNdata$ProbeSetName,]
  # print(str(input.data))
  if (!all(unique(RNdata$ProbeSetName == input.data$ProbeSetName))) stop(tmsg(paste0(track.type, " data and L2R data are not synched, or ordered differently !")), call. = FALSE)
  # ndata <- data.frame(chr = paste0("chr", input.data$chrs), start = input.data$pos, end = input.data$pos, name = rownames(input.data), RNdata[,-c(1:4), drop = FALSE], stringsAsFactors = FALSE)
  ndata <- data.frame(chr = input.data$chr, start = input.data$pos, end = input.data$pos, name = input.data$ProbeSetName, RNdata[,-c(1:4), drop = FALSE], stringsAsFactors = FALSE)
  rm(RNdata, renorm.data)
  # print(str(ndata))
  rm.diff <- diff(as.numeric(runmed(input.data$L2R[!is.na(input.data$L2R)], smo)))
  my.rm.mad <- sum(abs(rm.diff[rm.diff != 0]))
  # print(paste0("RMMAD ", my.rm.mad))
  # print(paste0(summary(my.rm.mad)))
  normloop.res <- list(data = input.data, renorm = l2r.fitloop(l2rObj = list(l2r=input.data$L2R, rm.mad = my.rm.mad), tfd = ndata, smo = smo))
  
  return(normloop.res)
  
  # input.data[[paste0("L2R.", pack.type)]] <- normloop.res$l2r$l2r + median(input.data$L2R, na.rm = TRUE)
  # input.data$L2R <- input.data[[paste0("L2R.", pack.type)]]
  # return(input.data)
  
}

#L2R f fit loop function
l2r.fitloop <- function(l2rObj, tfd, smo = 399, method = "loess")
  {
  
  ### FITLOOP
  minitf <- tfd[,-c(1:4), drop = FALSE]
  tfheads <- colnames(minitf)
  b <- ncol(minitf)+1
  posfit <- c()
  
  tmsg(paste0("Init (", l2rObj$rm.mad, ")"))
  
  while ( (b != 1) & (ncol(minitf) != 0) ) {
    
    biggy <- list()
    biggy <- append(biggy, list(l2rObj))
    rmtest <- l2rObj$rm.mad
    for (z in 1:length(minitf)) {
      if (method == "loess") tempfit <- l2r.fit(biggy[[1]]$l2r, minitf[,z], smo)
      if (method == "pcs") tempfit <- l2r.pcs(biggy[[1]]$l2r, minitf[,z], smo)
      biggy <- append(biggy, list(tempfit))
      rmtest <- c(rmtest, tempfit$rm.mad)
      # message(paste0(z, " / ", tempfit$rm.mad))
    }
    b <- which.min(rmtest)
    if (b > 1) {
      tmsg(paste0(" Positive fit with ", tfheads[b-1], " (", min(rmtest), ")"))
      posfit <- c(posfit, tfheads[b-1])
      l2rObj <- biggy[[b]]
      minitf <- as.data.frame(minitf[,-c(b-1)])
      tfheads <- tfheads[-c(b-1)]
    }
  }
  return(list(l2r = l2rObj, pos = posfit))
}

l2r.fit <- function(l2r, tf, smo)
  {
  l2fN <- limma::loessFit(l2r, tf)
  l2N <- l2r-l2fN$fitted
  rm.diff <- diff(as.numeric(runmed(l2N[!is.na(l2N)], smo)))
  Nrmspread <- sum(abs(rm.diff[rm.diff != 0]))
  return(list(l2r=l2N, rm.mad=Nrmspread))
}

### ASCAT ####
split_genome = function(SNPpos)
  {
  
  # look for gaps of more than 5Mb (arbitrary treshold to account for big centremeres or other gaps) and chromosome borders
  bigHoles = which(diff(SNPpos[,2])>=5000000)+1
  chrBorders = which(SNPpos[1:(dim(SNPpos)[1]-1),1]!=SNPpos[2:(dim(SNPpos)[1]),1])+1
  
  holes = unique(sort(c(bigHoles,chrBorders)))
  
  startseg = c(1,holes)
  endseg = c(holes-1,dim(SNPpos)[1])
  
  chr=list()
  for (i in 1:length(startseg)) {
    chr[[i]]=startseg[i]:endseg[i]
  }
  
  return(chr)
}

ascat.plotRawData = function(ASCATobj, img.dir=".", img.prefix="")
  {
  print.noquote("Plotting tumor data")
  for (i in 1:dim(ASCATobj$Tumor_LogR)[2]) {
    png(filename = file.path(img.dir, paste(img.prefix, ASCATobj$samples[i],".tumour.png",sep="")), width = 2000, height = 1000, res = 200)
    par(mar = c(0.5,5,5,0.5), mfrow = c(2,1), cex = 0.4, cex.main=3, cex.axis = 2, pch = ifelse(dim(ASCATobj$Tumor_LogR)[1]>100000,".",20))
    plot(c(1,dim(ASCATobj$Tumor_LogR)[1]), c(-1,1), type = "n", xaxt = "n", main = paste(ASCATobj$samples[i], ", tumor data, LogR", sep = ""), xlab = "", ylab = "")
    points(ASCATobj$Tumor_LogR[,i],col="red")
    #points(ASCATobj$Tumor_LogR[,i],col=rainbow(24)[ASCATobj$SNPpos$Chr])
    abline(v=0.5,lty=1,col="lightgrey")
    chrk_tot_len = 0
    for (j in 1:length(ASCATobj$ch)) {
      chrk = ASCATobj$ch[[j]];
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk)
      vpos = chrk_tot_len;
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
      text(tpos,1,ASCATobj$chrs[j], pos = 1, cex = 2)
      abline(v=vpos+0.5,lty=1,col="lightgrey")
    }
    plot(c(1,dim(ASCATobj$Tumor_BAF)[1]), c(0,1), type = "n", xaxt = "n", main = paste(ASCATobj$samples[i], ", tumor data, BAF", sep = ""), xlab = "", ylab = "")
    points(ASCATobj$Tumor_BAF[,i],col="red")
    abline(v=0.5,lty=1,col="lightgrey")
    chrk_tot_len = 0
    for (j in 1:length(ASCATobj$ch)) {
      chrk = ASCATobj$ch[[j]];
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk)
      vpos = chrk_tot_len;
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
      text(tpos,1,ASCATobj$chrs[j], pos = 1, cex = 2)
      abline(v=vpos+0.5,lty=1,col="lightgrey")
    }
    dev.off()
  }
  
  if(!is.null(ASCATobj$Germline_LogR)) {
    print.noquote("Plotting germline data")
    for (i in 1:dim(ASCATobj$Germline_LogR)[2]) {
      png(filename = file.path(img.dir, paste(img.prefix, ASCATobj$samples[i],".germline.png",sep="")), width = 2000, height = 1000, res = 200)
      par(mar = c(0.5,5,5,0.5), mfrow = c(2,1), cex = 0.4, cex.main=3, cex.axis = 2, pch = ifelse(dim(ASCATobj$Tumor_LogR)[1]>100000,".",20))
      plot(c(1,dim(ASCATobj$Germline_LogR)[1]), c(-1,1), type = "n", xaxt = "n", main = paste(ASCATobj$samples[i], ", germline data, LogR", sep = ""), xlab = "", ylab = "")
      points(ASCATobj$Germline_LogR[,i],col="red")
      abline(v=0.5,lty=1,col="lightgrey")
      chrk_tot_len = 0
      for (j in 1:length(ASCATobj$ch)) {
        chrk = ASCATobj$ch[[j]];
        chrk_tot_len_prev = chrk_tot_len
        chrk_tot_len = chrk_tot_len + length(chrk)
        vpos = chrk_tot_len;
        tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
        text(tpos,1,ASCATobj$chrs[j], pos = 1, cex = 2)
        abline(v=vpos+0.5,lty=1,col="lightgrey")
      }
      plot(c(1,dim(ASCATobj$Germline_BAF)[1]), c(0,1), type = "n", xaxt = "n", main = paste(ASCATobj$samples[i], ", germline data, BAF", sep = ""), xlab = "", ylab = "")
      points(ASCATobj$Germline_BAF[,i],col="red")
      abline(v=0.5,lty=1,col="lightgrey")
      chrk_tot_len = 0
      for (j in 1:length(ASCATobj$ch)) {
        chrk = ASCATobj$ch[[j]];
        chrk_tot_len_prev = chrk_tot_len
        chrk_tot_len = chrk_tot_len + length(chrk)
        vpos = chrk_tot_len;
        tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
        text(tpos,1,ASCATobj$chrs[j], pos = 1, cex = 2)
        abline(v=vpos+0.5,lty=1,col="lightgrey")
      }
      dev.off()
    }
  }
}

predictGermlineHomozygousStretches = function(chr, hom)
  {
  
  # contains the result: a list of vectors of probe numbers in homozygous stretches for each sample
  HomoStretches = list()
  
  for (sam in 1:dim(hom)[2]) {
    homsam = hom[,sam]
    
    perchom = sum(homsam,na.rm=T)/sum(!is.na(homsam))
    
    # NOTE THAT A P-VALUE THRESHOLD OF 0.001 IS HARDCODED HERE
    homthres = ceiling(log(0.001,perchom))
    
    allhprobes = NULL
    for (chrke in 1:length(chr)) {
      hschr = homsam[chr[[chrke]]]
      
      hprobes = vector(length=0)
      for(probe in 1:length(hschr)) {
        if(!is.na(hschr[probe])) {
          if(hschr[probe]) {
            hprobes = c(hprobes,probe)
          }
          else {
            if(length(hprobes)>=homthres) {
              allhprobes = rbind(allhprobes,c(chrke,chr[[chrke]][min(hprobes)],chr[[chrke]][max(hprobes)]))
            }
            hprobes = vector(length=0)
          }
        }
      }
      # if the last probe is homozygous, this is not yet accounted for
      if(!is.na(hschr[probe]) & hschr[probe]) {
        if(length(hprobes)>=homthres) {
          allhprobes = rbind(allhprobes,c(chrke,chr[[chrke]][min(hprobes)],chr[[chrke]][max(hprobes)]))
        }
      }
      
    }
    
    HomoStretches[[sam]]=allhprobes
    
  }
  
  return(HomoStretches)
}
#FOR RUN ASCAT FUNCTION
runASCAT = function(lrr, baf, lrrsegmented, bafsegmented, gender, SNPpos, chromosomes, chrnames, sexchromosomes, failedqualitycheck = F,
                    distancepng = NA, copynumberprofilespng = NA, nonroundedprofilepng = NA, aberrationreliabilitypng = NA, gamma = 0.55,
                    rho_manual = NA, psi_manual = NA, pdfPlot = F, y_limit = 5, circos=NA)
  {
  ch = chromosomes
  chrs = chrnames
  b = bafsegmented
  r = lrrsegmented[names(bafsegmented)]
  
  SNPposhet = SNPpos[names(bafsegmented),]
  autoprobes = !(SNPposhet[,1]%in%sexchromosomes)
  
  b2 = b[autoprobes]
  r2 = r[autoprobes]
  
  s = make_segments(r2,b2)
  d = create_distance_matrix(s, gamma)
  plot_d=d
  
  TheoretMaxdist = sum(rep(0.25,dim(s)[1]) * s[,"length"] * ifelse(s[,"b"]==0.5,0.05,1),na.rm=T)
  
  # flag the sample as non-aberrant if necessary
  nonaberrant = F
  MINABB = 0.03
  MINABBREGION = 0.005
  
  percentAbb = sum(ifelse(s[,"b"]==0.5,0,1)*s[,"length"])/sum(s[,"length"])
  maxsegAbb = max(ifelse(s[,"b"]==0.5,0,s[,"length"]))/sum(s[,"length"])
  if(percentAbb <= MINABB & maxsegAbb <= MINABBREGION) {
    nonaberrant = T
  }
  
  
  MINPLOIDY = 1.5
  MAXPLOIDY = 5.5
  MINRHO = 0.2
  MINGOODNESSOFFIT = 80
  MINPERCZERO = 0.02
  MINPERCZEROABB = 0.1
  MINPERCODDEVEN = 0.05
  MINPLOIDYSTRICT = 1.7
  MAXPLOIDYSTRICT = 2.3
  
  nropt = 0
  localmin = NULL
  optima = list()
  
  if(!failedqualitycheck && is.na(rho_manual)) {
    
    # first, try with all filters
    for (i in 4:(dim(d)[1]-3)) {
      for (j in 4:(dim(d)[2]-3)) {
        m = d[i,j]
        seld = d[(i-3):(i+3),(j-3):(j+3)]
        seld[4,4] = max(seld)
        if(min(seld) > m) {
          psi = as.numeric(rownames(d)[i])
          rho = as.numeric(colnames(d)[j])
          nA = (rho-1-(s[,"b"]-1)*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
          nB = (rho-1+s[,"b"]*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
          
          # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
          ploidy = sum((nA+nB) * s[,"length"]) / sum(s[,"length"]);
          
          percentzero = (sum((round(nA)==0)*s[,"length"])+sum((round(nB)==0)*s[,"length"]))/sum(s[,"length"])
          
          goodnessOfFit = (1-m/TheoretMaxdist) * 100
          
          if (!nonaberrant & ploidy > MINPLOIDY & ploidy < MAXPLOIDY & rho >= MINRHO & goodnessOfFit > MINGOODNESSOFFIT & percentzero > MINPERCZERO) {
            nropt = nropt + 1
            optima[[nropt]] = c(m,i,j,ploidy,goodnessOfFit)
            localmin[nropt] = m
          }
        }
      }
    }
    
    # if no solution, drop the percentzero > MINPERCZERO filter (allow non-aberrant solutions - but limit the ploidy options)
    if (nropt == 0) {
      for (i in 4:(dim(d)[1]-3)) {
        for (j in 4:(dim(d)[2]-3)) {
          m = d[i,j]
          seld = d[(i-3):(i+3),(j-3):(j+3)]
          seld[4,4] = max(seld)
          if(min(seld) > m) {
            psi = as.numeric(rownames(d)[i])
            rho = as.numeric(colnames(d)[j])
            nA = (rho-1-(s[,"b"]-1)*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
            nB = (rho-1+s[,"b"]*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
            
            # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
            ploidy = sum((nA+nB) * s[,"length"]) / sum(s[,"length"]);
            
            perczeroAbb = (sum((round(nA)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1))+sum((round(nB)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1)))/sum(s[,"length"]*ifelse(s[,"b"]==0.5,0,1))
            # the next can happen if BAF is a flat line at 0.5
            if (is.na(perczeroAbb)) {
              perczeroAbb = 0
            }
            
            goodnessOfFit = (1-m/TheoretMaxdist) * 100
            
            if (ploidy > MINPLOIDYSTRICT & ploidy < MAXPLOIDYSTRICT & rho >= MINRHO & goodnessOfFit > MINGOODNESSOFFIT & perczeroAbb > MINPERCZEROABB) {
              nropt = nropt + 1
              optima[[nropt]] = c(m,i,j,ploidy,goodnessOfFit)
              localmin[nropt] = m
            }
          }
        }
      }
    }
    
    # if still no solution, allow solutions with 100% aberrant cells (include the borders with rho = 1), but in first instance, keep the percentzero > 0.01 filter
    if (nropt == 0) {
      #first, include borders
      cold = which(as.numeric(colnames(d))>1)
      d[,cold]=1E20
      for (i in 4:(dim(d)[1]-3)) {
        for (j in 4:(dim(d)[2]-3)) {
          m = d[i,j]
          seld = d[(i-3):(i+3),(j-3):(j+3)]
          seld[4,4] = max(seld)
          if(min(seld) > m) {
            psi = as.numeric(rownames(d)[i])
            rho = as.numeric(colnames(d)[j])
            nA = (rho-1-(s[,"b"]-1)*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
            nB = (rho-1+s[,"b"]*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
            
            # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
            ploidy = sum((nA+nB) * s[,"length"]) / sum(s[,"length"]);
            
            percentzero = (sum((round(nA)==0)*s[,"length"])+sum((round(nB)==0)*s[,"length"]))/sum(s[,"length"])
            percOddEven = sum((round(nA)%%2==0&round(nB)%%2==1|round(nA)%%2==1&round(nB)%%2==0)*s[,"length"])/sum(s[,"length"])
            perczeroAbb = (sum((round(nA)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1))+sum((round(nB)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1)))/sum(s[,"length"]*ifelse(s[,"b"]==0.5,0,1))
            if (is.na(perczeroAbb)) {
              perczeroAbb = 0
            }
            
            goodnessOfFit = (1-m/TheoretMaxdist) * 100
            
            if (!nonaberrant & ploidy > MINPLOIDY & ploidy < MAXPLOIDY & rho >= MINRHO & goodnessOfFit > MINGOODNESSOFFIT &
                (perczeroAbb > MINPERCZEROABB | percentzero > MINPERCZERO | percOddEven > MINPERCODDEVEN)) {
              nropt = nropt + 1
              optima[[nropt]] = c(m,i,j,ploidy,goodnessOfFit)
              localmin[nropt] = m
            }
          }
        }
      }
    }
    
    # if still no solution, drop the percentzero > MINPERCENTZERO filter, but strict ploidy borders
    if (nropt == 0) {
      for (i in 4:(dim(d)[1]-3)) {
        for (j in 4:(dim(d)[2]-3)) {
          m = d[i,j]
          seld = d[(i-3):(i+3),(j-3):(j+3)]
          seld[4,4] = max(seld)
          if(min(seld) > m) {
            psi = as.numeric(rownames(d)[i])
            rho = as.numeric(colnames(d)[j])
            nA = (rho-1-(s[,"b"]-1)*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
            nB = (rho-1+s[,"b"]*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
            
            # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
            ploidy = sum((nA+nB) * s[,"length"]) / sum(s[,"length"]);
            
            perczeroAbb = (sum((round(nA)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1))+sum((round(nB)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1)))/sum(s[,"length"]*ifelse(s[,"b"]==0.5,0,1))
            # the next can happen if BAF is a flat line at 0.5
            if (is.na(perczeroAbb)) {
              perczeroAbb = 0
            }
            
            goodnessOfFit = (1-m/TheoretMaxdist) * 100
            
            if (ploidy > MINPLOIDYSTRICT & ploidy < MAXPLOIDYSTRICT & rho >= MINRHO & goodnessOfFit > MINGOODNESSOFFIT) {
              nropt = nropt + 1
              optima[[nropt]] = c(m,i,j,ploidy,goodnessOfFit)
              localmin[nropt] = m
            }
          }
        }
      }
    }
  }
  
  if (!is.na(rho_manual)) {
    
    rho = rho_manual
    psi = psi_manual
    
    nA = (rho-1-(s[,"b"]-1)*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
    nB = (rho-1+s[,"b"]*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
    
    # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
    ploidy = sum((nA+nB) * s[,"length"]) / sum(s[,"length"]);
    
    nMinor = NULL
    if (sum(nA,na.rm=T) < sum(nB,na.rm=T)) {
      nMinor = nA
    }
    else {
      nMinor = nB
    }
    m = sum(abs(nMinor - pmax(round(nMinor),0))^2 * s[,"length"] * ifelse(s[,"b"]==0.5,0.05,1), na.rm=T)
    
    goodnessOfFit = (1-m/TheoretMaxdist) * 100
    
    nropt = 1
    optima[[1]] = c(m,rho,psi,ploidy,goodnessOfFit)
    localmin[1] = m
    
  }
  
  
  if (nropt>0) {
    if (is.na(rho_manual)) {
      optlim = sort(localmin)[1]
      for (i in 1:length(optima)) {
        if(optima[[i]][1] == optlim) {
          psi_opt1 = as.numeric(rownames(d)[optima[[i]][2]])
          rho_opt1 = as.numeric(colnames(d)[optima[[i]][3]])
          if(rho_opt1 > 1) {
            rho_opt1 = 1
          }
          ploidy_opt1 = optima[[i]][4]
          goodnessOfFit_opt1 = optima[[i]][5]
        }
      }
    } else {
      rho_opt1 = optima[[1]][2]
      psi_opt1 = optima[[1]][3]
      ploidy_opt1 = optima[[1]][4]
      goodnessOfFit_opt1 = optima[[1]][5]
    }
  }
  
  if(nropt>0) {
    #plot Sunrise
    if (!is.na(distancepng)) {
      png(filename = distancepng, width = 1000, height = 1000, res = 1000/7)
    }
    ascat.plotSunrise(plot_d,psi_opt1,rho_opt1)
    if (!is.na(distancepng)) {
      dev.off()
    }
    
    rho = rho_opt1
    psi = psi_opt1
    SNPposhet = SNPpos[names(bafsegmented),]
    haploidchrs = unique(c(substring(gender,1,1),substring(gender,2,2)))
    if(substring(gender,1,1)==substring(gender,2,2)) {
      haploidchrs = setdiff(haploidchrs,substring(gender,1,1))
    }
    diploidprobes = !(SNPposhet[,1]%in%haploidchrs)
    nullchrs = setdiff(sexchromosomes,unique(c(substring(gender,1,1),substring(gender,2,2))))
    nullprobes = SNPposhet[,1]%in%nullchrs
    
    nAfull = ifelse(diploidprobes,
                    (rho-1-(b-1)*2^(r/gamma)*((1-rho)*2+rho*psi))/rho,
                    ifelse(nullprobes,0,
                           ifelse(b<0.5,(rho-1+((1-rho)*2+rho*psi)*2^(r/gamma))/rho,0)))
    nBfull = ifelse(diploidprobes,
                    (rho-1+b*2^(r/gamma)*((1-rho)*2+rho*psi))/rho,
                    ifelse(nullprobes,0,
                           ifelse(b<0.5,0,(rho-1+((1-rho)*2+rho*psi)*2^(r/gamma))/rho)))
    nA = pmax(round(nAfull),0)
    nB = pmax(round(nBfull),0)
    
    if(!is.na(circos)){
      frame<-cbind(SNPposhet,nAfull,nBfull)
      chrSegmA<-rle(frame$nAfull)
      chrSegmB<-rle(frame$nBfull)
      if(all(chrSegmA$lengths==chrSegmB$lengths)){
        start=1
        for(i in 1:length(chrSegmA$values)){
          valA<-chrSegmA$values[i]
          valB<-chrSegmB$values[i]
          size<-chrSegmA$lengths[i]
          write(c(paste("hs",frame[start,1],sep=""),frame[start,2],frame[(start+size-1),2],valA), file = paste(circos,"_major",sep=""), ncolumns = 4, append = TRUE, sep = "\t")
          write(c(paste("hs",frame[start,1],sep=""),frame[start,2],frame[(start+size-1),2],valB), file = paste(circos,"_minor",sep=""), ncolumns = 4, append = TRUE, sep = "\t")
          start=start+size
        }
      }
      else{
        print("Major and minor allele copy numbers are segmented differently.")
      }
    }
    
    if (is.na(nonroundedprofilepng)) {
      dev.new(10,5)
    }
    else {
      if(pdfPlot){
        pdf(file = nonroundedprofilepng, width = 20, height = y_limit, pointsize=20)
      }
      else{
        png(filename = nonroundedprofilepng, width = 2000, height = (y_limit*100), res = 200)
      }
    }
    
    ascat.plotNonRounded(ploidy_opt1, rho_opt1, goodnessOfFit_opt1, nonaberrant, nAfull, nBfull, y_limit, bafsegmented, ch,lrr, chrnames)
    
    if (!is.na(nonroundedprofilepng)) {
      dev.off()
    }
    
    rho = rho_opt1
    psi = psi_opt1
    
    diploidprobes = !(SNPpos[,1]%in%haploidchrs)
    nullprobes = SNPpos[,1]%in%nullchrs
    
    #this replaces an occurrence of unique that caused problems
    #introduces segment spanning over chr ends, when two consecutive probes from diff chr have same logR!
    # build helping vector
    chrhelp = vector(length=length(lrrsegmented))
    for (chrnr in 1:length(ch)) {
      chrke = ch[[chrnr]]
      chrhelp[chrke] = chrnr
    }
    
    tlr2 = rle(lrrsegmented)
    tlr.chr= rle(chrhelp)
    
    tlrstart = c(1,cumsum(tlr2$lengths)+1)
    tlrstart = tlrstart[1:(length(tlrstart)-1)]
    tlrend = cumsum(tlr2$lengths)
    
    tlrstart.chr= c(1,cumsum(tlr.chr$lengths)+1)
    tlrstart.chr = tlrstart.chr[1:(length(tlrstart.chr)-1)]
    tlrend.chr = cumsum(tlr.chr$lengths)
    
    tlrend<-sort(union(tlrend, tlrend.chr))
    tlrstart<-sort(union(tlrstart, tlrstart.chr))
    
    tlr=NULL
    for(ind in tlrstart){
      val<-lrrsegmented[ind]
      tlr<-c(tlr, val)
    }
    
    # For each LRR probe, find the matching BAF probe
    # and its position in bafsegmented
    probeLookup = data.frame(
      lrrprobe = names(lrrsegmented),
      bafpos = match(names(lrrsegmented), names(bafsegmented)),
      stringsAsFactors=F
    )
    
    seg = NULL
    for (i in 1:length(tlr)) {
      logR = tlr[i]
      #pr = which(lrrsegmented==logR) # this was a problem
      pr = tlrstart[i]:tlrend[i]
      start = min(pr)
      end = max(pr)
      
      bafpos = probeLookup$bafpos[pr]
      bafpos = bafpos[!is.na(bafpos)]
      bafke  = bafsegmented[bafpos][1]
      
      #if bafke is NA, this means that we are dealing with a germline homozygous stretch with a copy number change within it.
      #in this case, nA and nB are irrelevant, just their sum matters
      if(is.na(bafke)) {
        bafke=0
      }
      
      nAraw = ifelse(diploidprobes[start],
                     (rho-1-(bafke-1)*2^(logR/gamma)*((1-rho)*2+rho*psi))/rho,
                     ifelse(nullprobes[start],0,
                            (rho-1+((1-rho)*2+rho*psi)*2^(logR/gamma))/rho))
      nBraw = ifelse(diploidprobes[start],(rho-1+bafke*2^(logR/gamma)*((1-rho)*2+rho*psi))/rho,0)
      # correct for negative values:
      if (nAraw+nBraw<0) {
        nAraw = 0
        nBraw = 0
      }
      else if (nAraw<0) {
        nBraw = nAraw+nBraw
        nAraw = 0
      }
      else if (nBraw<0) {
        nAraw = nAraw+nBraw
        nBraw = 0
      }
      # when evidence for odd copy number in segments of BAF = 0.5, assume a deviation..
      limitround = 0.5
      nA = ifelse(bafke==0.5,
                  ifelse(nAraw+nBraw>round(nAraw)+round(nBraw)+limitround,
                         round(nAraw)+1,
                         ifelse(nAraw+nBraw<round(nAraw)+round(nBraw)-limitround,
                                round(nAraw),
                                round(nAraw))),
                  round(nAraw))
      nB = ifelse(bafke==0.5,
                  ifelse(nAraw+nBraw>round(nAraw)+round(nBraw)+limitround,
                         round(nBraw),
                         ifelse(nAraw+nBraw<round(nAraw)+round(nBraw)-limitround,
                                round(nBraw)-1,
                                round(nBraw))),
                  round(nBraw))
      if (is.null(seg)) {
        seg = t(as.matrix(c(start,end,nA,nB)))
        seg_raw = t(as.matrix(c(start,end,nA,nB,nAraw,nBraw)))
      }
      else {
        seg = rbind(seg,c(start,end,nA,nB))
        seg_raw = rbind(seg_raw,c(start,end,nA,nB,nAraw,nBraw))
      }
    }
    colnames(seg)=c("start","end","nA","nB")
    colnames(seg_raw)=c("start","end","nA","nB","nAraw","nBraw")
    
    # every repeat joins 2 ends. 20 repeats will join about 1 million ends..
    for (rep in 1:20) {
      seg2=seg
      seg = NULL
      skipnext = F
      for(i in 1:dim(seg2)[1]) {
        if(!skipnext) {
          if(i != dim(seg2)[1] && seg2[i,"nA"]==seg2[i+1,"nA"] && seg2[i,"nB"]==seg2[i+1,"nB"] &&
             chrhelp[seg2[i,"end"]]==chrhelp[seg2[i+1,"start"]]) {
            segline = c(seg2[i,"start"],seg2[i+1,"end"],seg2[i,3:4])
            skipnext = T
          }
          else {
            segline = seg2[i,]
          }
          
          if (is.null(seg)) {
            seg = t(as.matrix(segline))
          }
          else {
            seg = rbind(seg,segline)
          }
        }
        else {
          skipnext = F
        }
      }
      colnames(seg)=colnames(seg2)
    }
    rownames(seg)=NULL
    
    nMajor = vector(length = length(lrrsegmented))
    names(nMajor) = names(lrrsegmented)
    nMinor = vector(length = length(lrrsegmented))
    names(nMinor) = names(lrrsegmented)
    
    for (i in 1:dim(seg)[1]) {
      nMajor[seg[i,"start"]:seg[i,"end"]] = seg[i,"nA"]
      nMinor[seg[i,"start"]:seg[i,"end"]] = seg[i,"nB"]
    }
    
    n1all = vector(length = length(lrrsegmented))
    names(n1all) = names(lrrsegmented)
    n2all = vector(length = length(lrrsegmented))
    names(n2all) = names(lrrsegmented)
    
    # note: any of these can have length 0
    NAprobes = which(is.na(lrr))
    heteroprobes = setdiff(which(names(lrrsegmented)%in%names(bafsegmented)),NAprobes)
    homoprobes = setdiff(setdiff(which(!is.na(baf)),heteroprobes),NAprobes)
    CNprobes = setdiff(which(is.na(baf)),NAprobes)
    
    n1all[NAprobes] = NA
    n2all[NAprobes] = NA
    n1all[CNprobes] = nMajor[CNprobes]+nMinor[CNprobes]
    n2all[CNprobes] = 0
    heteroprobes2 = names(lrrsegmented)[heteroprobes]
    n1all[heteroprobes] = ifelse(baf[heteroprobes2]<=0.5,nMajor[heteroprobes], nMinor[heteroprobes])
    n2all[heteroprobes] = ifelse(baf[heteroprobes2]>0.5,nMajor[heteroprobes], nMinor[heteroprobes])
    n1all[homoprobes] = ifelse(baf[homoprobes]<=0.5,nMajor[homoprobes]+nMinor[homoprobes],0)
    n2all[homoprobes] = ifelse(baf[homoprobes]>0.5,nMajor[homoprobes]+nMinor[homoprobes],0)
    
    
    # plot ASCAT profile
    if (is.na(copynumberprofilespng)) {
      dev.new(10,2.5)
    }
    else {
      if(pdfPlot){
        pdf(file = copynumberprofilespng, width = 20, height = y_limit, pointsize=20)
      }
      else{
        png(filename = copynumberprofilespng, width = 2000, height = (y_limit*100), res = 200)
      }
    }
    #plot ascat profile
    ascat.plotAscatProfile(n1all, n2all, heteroprobes, ploidy_opt1, rho_opt1, goodnessOfFit_opt1, nonaberrant,y_limit, ch, lrr, bafsegmented, chrnames)
    
    
    if (!is.na(copynumberprofilespng)) {
      dev.off()
    }
    
    
    if (!is.na(aberrationreliabilitypng)) {
      png(filename = aberrationreliabilitypng, width = 2000, height = 500, res = 200)
      par(mar = c(0.5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2.5)
      
      diploidprobes = !(SNPposhet[,1]%in%haploidchrs)
      nullprobes = SNPposhet[,1]%in%nullchrs
      
      rBacktransform = ifelse(diploidprobes,
                              gamma*log((rho*(nA+nB)+(1-rho)*2)/((1-rho)*2+rho*psi),2),
                              # the value for nullprobes is arbitrary (but doesn't matter, as these are not plotted anyway because BAF=0.5)
                              ifelse(nullprobes,-10,gamma*log((rho*(nA+nB)+(1-rho))/((1-rho)*2+rho*psi),2)))
      
      bBacktransform = ifelse(diploidprobes,
                              (1-rho+rho*nB)/(2-2*rho+rho*(nA+nB)),
                              ifelse(nullprobes,0.5,0))
      
      rConf = ifelse(abs(rBacktransform)>0.15,pmin(100,pmax(0,100*(1-abs(rBacktransform-r)/abs(r)))),NA)
      bConf = ifelse(diploidprobes & bBacktransform!=0.5, pmin(100,pmax(0,ifelse(b==0.5,100,100*(1-abs(bBacktransform-b)/abs(b-0.5))))), NA)
      confidence = ifelse(is.na(rConf),bConf,ifelse(is.na(bConf),rConf,(rConf+bConf)/2))
      maintitle = paste("Aberration reliability score (%), average: ", sprintf("%2.0f",mean(confidence,na.rm=T)),"%",sep="")
      plot(c(1,length(nAfull)), c(0,100), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
      points(confidence,col="blue",pch = "|")
      abline(v=0,lty=1,col="lightgrey")
      chrk_tot_len = 0
      for (i in 1:length(ch)) {
        chrk = ch[[i]];
        chrk_hetero = intersect(names(lrr)[chrk],names(bafsegmented))
        chrk_tot_len_prev = chrk_tot_len
        chrk_tot_len = chrk_tot_len + length(chrk_hetero)
        vpos = chrk_tot_len;
        tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
        text(tpos,5,chrs[i], pos = 1, cex = 2)
        abline(v=vpos,lty=1,col="lightgrey")
      }
      dev.off()
    }
    
    return(list(rho = rho_opt1, psi = psi_opt1, goodnessOfFit = goodnessOfFit_opt1, nonaberrant = nonaberrant,
                nA = n1all, nB = n2all, seg = seg, seg_raw = seg_raw, distance_matrix = d))
    
  }
  
  else {
    
    name=gsub(".sunrise.png","",basename(distancepng))
    
    png(filename = distancepng, width = 1000, height = 1000, res = 1000/7)
    ascat.plotSunrise(plot_d,0,0)
    dev.off()
    
    warning(paste("ASCAT could not find an optimal ploidy and cellularity value for sample ", name, ".\n", sep=""))
    return(list(rho = NA, psi = NA, goodnessOfFit = NA, nonaberrant = F, nA = NA, nB = NA, seg = NA, seg_raw = NA, distance_matrix = NA))
  }
  
}

make_segments = function(r,b)
  {
  m = matrix(ncol = 2, nrow = length(b))
  m[,1] = r
  m[,2] = b
  m = as.matrix(na.omit(m))
  pcf_segments = matrix(ncol = 3, nrow = dim(m)[1])
  colnames(pcf_segments) = c("r","b","length");
  index = 0;
  previousb = -1;
  previousr = 1E10;
  for (i in 1:dim(m)[1]) {
    if (m[i,2] != previousb || m[i,1] != previousr) {
      index=index+1;
      count=1;
      pcf_segments[index, "r"] = m[i,1];
      pcf_segments[index, "b"] = m[i,2];
    }
    else {
      count = count + 1;
    }
    pcf_segments[index, "length"] = count;
    previousb = m[i,2];
    previousr = m[i,1];
  }
  pcf_segments = as.matrix(na.omit(pcf_segments))[,,drop=FALSE]
  return(pcf_segments);
}

create_distance_matrix = function(segments, gamma)
  {
  s = segments
  psi_pos = seq(1,6,0.05)
  rho_pos = seq(0.1,1.05,0.01)
  d = matrix(nrow = length(psi_pos), ncol = length(rho_pos))
  rownames(d) = psi_pos
  colnames(d) = rho_pos
  dmin = 1E20;
  for(i in 1:length(psi_pos)) {
    psi = psi_pos[i]
    for(j in 1:length(rho_pos)) {
      rho = rho_pos[j]
      nA = (rho-1-(s[,"b"]-1)*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
      nB = (rho-1+s[,"b"]*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
      # choose the minor allele
      nMinor = NULL
      if (sum(nA,na.rm=T) < sum(nB,na.rm=T)) {
        nMinor = nA
      }
      else {
        nMinor = nB
      }
      d[i,j] = sum(abs(nMinor - pmax(round(nMinor),0))^2 * s[,"length"] * ifelse(s[,"b"]==0.5,0.05,1), na.rm=T)
    }
  }
  return(d)
}

ascat.plotSunrise<-function(d, psi_opt1, rho_opt1, minim=T)
  {
  
  par(mar = c(5,5,0.5,0.5), cex=0.75, cex.lab=2, cex.axis=2)
  
  if(minim){
    hmcol = rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256))
  } else {
    hmcol = colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256)
  } 
  image(log(d), col = hmcol, axes = F, xlab = "Ploidy", ylab = "Aberrant cell fraction")
  
  ploidy_min<-as.numeric(rownames(d)[1])
  ploidy_max<-as.numeric(rownames(d)[nrow(d)])
  purity_min<-as.numeric(colnames(d)[1])
  purity_max<-as.numeric(colnames(d)[ncol(d)])
  
  axis(1, at = seq(0, 1, by = 1/(ploidy_max-1)), labels = seq(ploidy_min, ploidy_max, by = 1))
  axis(2, at = seq(0, 1/purity_max, by = 1/3/purity_max), labels = seq(purity_min, purity_max, by = 3/10))
  
  if(psi_opt1>0 && rho_opt1>0){
    points((psi_opt1-ploidy_min)/(ploidy_max-1),(rho_opt1-purity_min)/(1/purity_max),col="green",pch=4, cex = 2)
  }
}

ascat.plotNonRounded <- function(ploidy, rho, goodnessOfFit, nonaberrant, nAfull, nBfull,y_limit=5,bafsegmented,ch,lrr, chrs)
  {
  maintitle = paste("Ploidy: ",sprintf("%1.2f",ploidy),", aberrant cell fraction: ",sprintf("%2.0f",rho*100),"%, goodness of fit: ",sprintf("%2.1f",goodnessOfFit),"%", ifelse(nonaberrant,", non-aberrant",""),sep="")
  nBfullPlot<-ifelse(nBfull<y_limit, nBfull, y_limit+0.1)
  nAfullPlot<-ifelse((nAfull+nBfull)<y_limit, nAfull+nBfull, y_limit+0.1)
  colourTotal = "purple"
  colourMinor = "blue"
  base.gw.plot(bafsegmented,nAfullPlot,nBfullPlot,colourTotal,colourMinor,maintitle,ch,lrr,chrs,y_limit,twoColours=TRUE)
}

base.gw.plot = function(bafsegmented,nAfullPlot,nBfullPlot,colourTotal,colourMinor,maintitle,chr.segs,lrr,chr.names,y_limit,twoColours=FALSE)
  {
  par(mar = c(0.5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2.5)
  ticks=seq(0, y_limit, 1)
  plot(c(1,length(nAfullPlot)), c(0,y_limit), type = "n", xaxt = "n", yaxt="n", main = maintitle, xlab = "", ylab = "")
  axis(side = 2, at = ticks)
  abline(h=ticks, col="lightgrey", lty=1)
  
  A_rle<-rle(nAfullPlot)
  start=0
  #plot total copy number
  for(i in 1:length(A_rle$values)){
    val<-A_rle$values[i]
    size<-A_rle$lengths[i]
    rect(start, (val-0.07), (start+size-1), (val+0.07), col=ifelse((twoColours & val>=y_limit), adjustcolor(colourTotal,red.f=0.75,green.f=0.75,blue.f=0.75), colourTotal), border=ifelse((twoColours & val>=y_limit), adjustcolor(colourTotal,red.f=0.75,green.f=0.75,blue.f=0.75), colourTotal))
    start=start+size
  }
  
  B_rle<-rle(nBfullPlot)
  start=0
  #plot minor allele copy number
  for(i in 1:length(B_rle$values)){
    val<-B_rle$values[i]
    size<-B_rle$lengths[i]
    rect(start, (val-0.07), (start+size-1), (val+0.07), col=ifelse((twoColours & val>=y_limit), adjustcolor(colourMinor, red.f=0.75,green.f=0.75,blue.f=0.75), colourMinor), border=ifelse((twoColours & val>=y_limit), adjustcolor(colourMinor, red.f=0.75,green.f=0.75,blue.f=0.75), colourMinor))
    start=start+size
  }
  
  chrk_tot_len = 0
  abline(v=0,lty=1,col="lightgrey")
  for (i in 1:length(chr.segs)) {
    chrk = chr.segs[[i]];
    chrk_hetero = intersect(names(lrr)[chrk],names(bafsegmented))
    chrk_tot_len_prev = chrk_tot_len
    chrk_tot_len = chrk_tot_len + length(chrk_hetero)
    vpos = chrk_tot_len;
    tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
    text(tpos,y_limit,chr.names[i], pos = 1, cex = 2)
    abline(v=vpos,lty=1,col="lightgrey")
  }
  
  #   #add text to too high fragments
  #   if(textFlag){
  #     #      rleB<-rle(nBfullPlot>y_limit)
  #     #      pos<-0
  #     #      for(i in 1:length(rleB$values)){
  #     #        if(rleB$values[i]){
  #     #          xpos=pos+(rleB$lengths[i]/2)
  #     #          text(xpos,y_limit+0.1,sprintf("%1.2f",nBfull[pos+1]), pos = 1, cex = 0.7)
  #     #        }
  #     #        pos=pos+rleB$lengths[i]
  #     #      }
  #
  #     rleA<-rle(nAfullPlot>y_limit)
  #     pos<-0
  #     for(i in 1:length(rleA$values)){
  #       if(rleA$values[i]){
  #         xpos=pos+(rleA$lengths[i]/2)
  #         text(xpos,y_limit+0.1,sprintf("%1.2f",(nAfull[pos+1]+nBfull[pos+1])), pos = 1, cex = 0.7)
  #       }
  #       pos=pos+rleA$lengths[i]
  #     }
  #   }
}

ascat.plotAscatProfile<-function(n1all, n2all, heteroprobes, ploidy, rho, goodnessOfFit, nonaberrant, y_limit=5, ch, lrr, bafsegmented, chrs)
  {
  nA2 = n1all[heteroprobes]
  nB2 = n2all[heteroprobes]
  nA = ifelse(nA2>nB2,nA2,nB2)
  nB = ifelse(nA2>nB2,nB2,nA2)
  
  nBPlot<-ifelse(nB<=y_limit, nB+0.1, y_limit+0.1)
  nAPlot<-ifelse(nA<=y_limit, nA-0.1, y_limit+0.1)
  
  colourTotal="red"
  colourMinor="green"
  
  maintitle = paste("Ploidy: ",sprintf("%1.2f",ploidy),", aberrant cell fraction: ",sprintf("%2.0f",rho*100),"%, goodness of fit: ",sprintf("%2.1f",goodnessOfFit),"%", ifelse(nonaberrant,", non-aberrant",""),sep="")
  base.gw.plot(bafsegmented,nAPlot,nBPlot,colourTotal,colourMinor,maintitle,ch,lrr,chrs,y_limit,twoColours=TRUE)
}

madWins <- function(x,tau,k)
  {
  xhat <- medianFilter(x,k)
  d <- x-xhat
  SD <- mad(d)
  z <- tau*SD
  xwin <- xhat + psi(d, z)
  outliers <- rep(0, length(x))
  outliers[x > xwin] <- 1
  outliers[x < xwin] <- -1
  return(list(ywin=xwin,sdev=SD,outliers=outliers))
}

medianFilter <- function(x,k)
  {
  n <- length(x)
  filtWidth <- 2*k + 1
  
  #Make sure filtWidth does not exceed n
  if(filtWidth > n){
    if(n==0){
      filtWidth <- 1
    }else if(n%%2 == 0){
      #runmed requires filtWidth to be odd, ensure this:
      filtWidth <- n - 1
    }else{
      filtWidth <- n
    }
  }
  
  runMedian <- runmed(x,k=filtWidth,endrule="median")
  
  return(runMedian)
  
}

psi <- function(x,z)
  {
  xwin <- x
  xwin[x < -z] <- -z
  xwin[x > z] <- z
  return(xwin)
}

fastAspcf <- function(logR, allB, kmin, gamma)
  {
  
  N <- length(logR)
  w <- 1000 #w: windowsize
  d <- 100
  
  startw = -d
  stopw = w-d
  
  nseg = 0
  var2 = 0
  var3 = 0
  breakpts = 0
  larger = TRUE
  repeat{
    from <- max(c(1,startw))
    to  <-  min(c(stopw,N))
    logRpart <- logR[from:to]
    allBpart <- allB[from:to]
    allBflip <- allBpart
    allBflip[allBpart > 0.5] <- 1 - allBpart[allBpart > 0.5]
    
    sd1 <- getMad(logRpart)
    sd2 <- getMad(allBflip)
    sd3 <- getMad(allBpart) 
    
    #Must check that sd1 and sd2 are defined and != 0:
    sd.valid <- c(!is.na(sd1),!is.na(sd2),sd1!=0,sd2!=0)
    if(all(sd.valid)){
      #run aspcfpart:
      #part.res <- aspcfpart(logRpart=logRpart, allBflip=allBflip, a=startw, b=stopw, d=d, sd1=sd1, sd2=sd2, N=N, kmin=kmin, gamma=gamma)
      part.res <- aspcfpart(logRpart=logRpart, allBflip=allBflip, a=startw, b=stopw, d=d, sd1=sd1, sd2=sd2, N=N, kmin=kmin, gamma=gamma)
      breakptspart <- part.res$breakpts
      # the 'larger' is (occasionally) necessary in the last window of the segmentation!
      larger = breakptspart>breakpts[length(breakpts)]
      breakpts <- c(breakpts, breakptspart[larger])
      var2 <- var2 + sd2^2
      var3 <- var3 + sd3^2
      nseg = nseg+1
    }
    
    if(stopw < N+d){
      startw <- min(stopw-2*d + 1,N-2*d)
      stopw <- startw + w
    }else{
      break
    }
    
  }#end repeat
  breakpts <- unique(c(breakpts, N))
  if(nseg==0){nseg=1}  #just in case the sd-test never passes.
  sd2 <- sqrt(var2/nseg)
  sd3 <- sqrt(var3/nseg)
  
  # On each segment calculate mean of unflipped B allele data
  frst <- breakpts[1:length(breakpts)-1] + 1
  last <- breakpts[2:length(breakpts)]
  nseg <- length(frst)
  
  yhat1 <- rep(NA,N)
  yhat2 <- rep(NA,N)
  
  for(i in 1:nseg){
    yhat1[frst[i]:last[i]] <- rep(mean(logR[frst[i]:last[i]]), last[i]-frst[i]+1)
    yi2 <- allB[frst[i]:last[i]]
    # Center data around zero (by subtracting 0.5) and estimate mean
    if(length(yi2)== 0){
      mu <- 0
    }else{
      mu <- mean(abs(yi2-0.5))
    }
    
    # Make a (slightly arbitrary) decision concerning branches
    # This may be improved by a test of equal variances
    if(sqrt(sd2^2+mu^2) < 2*sd2){
      # if(sd3 < 1.8*sd2){
      mu <- 0
    }
    yhat2[frst[i]:last[i]] <- rep(mu+0.5,last[i]-frst[i]+1)
  }
  
  return(list(yhat1=yhat1,yhat2=yhat2))
  
}

getMad <- function(x,k=25)
  {
  
  #Remove observations that are equal to zero; are likely to be imputed, should not contribute to sd:
  x <- x[x!=0]
  
  #Calculate runMedian
  runMedian <- medianFilter(x,k)
  
  dif <- x-runMedian
  SD <- mad(dif)
  
  return(SD)
}

aspcfpart <- function(logRpart, allBflip, a, b, d, sd1, sd2, N, kmin, gamma)
  {
  
  from <- max(c(1,a))
  usefrom <- max(c(1,a+d))
  useto <- min(c(N,b-d))
  
  N <- length(logRpart)
  y1 <- logRpart
  y2 <- allBflip
  
  #Check that vectors are long enough to run algorithm:
  if(N < 2*kmin){
    breakpts <- 0
    return(list(breakpts=breakpts))
  }
  
  # Find initSum, initKvad, initAve for segment y[1..kmin]
  initSum1 <- sum(y1[1:kmin])
  initKvad1 <- sum(y1[1:kmin]^2)
  initAve1 <- initSum1/kmin
  initSum2 <- sum(y2[1:kmin])
  initKvad2 <- sum(y2[1:kmin]^2)
  initAve2 <- initSum2/kmin
  
  # Define vector of best costs
  bestCost <- rep(0,N)
  cost1 <- (initKvad1 - initSum1*initAve1)/sd1^2
  cost2 <- (initKvad2 - initSum2*initAve2)/sd2^2
  bestCost[kmin] <- cost1 + cost2
  
  # Define vector of best splits
  bestSplit <- rep(0,N)
  
  # Define vector of best averages
  bestAver1 <- rep(0,N)
  bestAver2 <- rep(0,N)
  bestAver1[kmin] <- initAve1
  bestAver2[kmin] <- initAve2
  
  
  #Initialize
  Sum1 <- rep(0,N)
  Sum2 <- rep(0,N)
  Kvad1 <- rep(0,N)
  Kvad2 <- rep(0,N)
  Aver1 <- rep(0,N)
  Aver2 <- rep(0,N)
  Cost <- rep(0,N)
  
  # We have to treat the region y(1..2*kmin-1) separately, as it
  # cannot be split into two full segments
  kminP1 <- kmin+1
  for (k in (kminP1):(2*kmin-1)) {
    Sum1[kminP1:k] <- Sum1[kminP1:k]+y1[k]
    Aver1[kminP1:k] <- Sum1[kminP1:k]/((k-kmin):1)
    Kvad1[kminP1:k] <- Kvad1[kminP1:k]+y1[k]^2
    Sum2[kminP1:k] <- Sum2[kminP1:k]+y2[k]
    Aver2[kminP1:k] <- Sum2[kminP1:k]/((k-kmin):1)
    Kvad2[kminP1:k] <- Kvad2[kminP1:k]+y2[k]^2
    
    
    bestAver1[k] <- (initSum1+Sum1[kminP1])/k
    bestAver2[k] <- (initSum2+Sum2[kminP1])/k
    cost1 <- ((initKvad1+Kvad1[kminP1])-k*bestAver1[k]^2)/sd1^2
    cost2 <- ((initKvad2+Kvad2[kminP1])-k*bestAver2[k]^2)/sd2^2
    
    bestCost[k] <- cost1 + cost2
    
  }
  
  
  for (n in (2*kmin):N) {
    
    nMkminP1=n-kmin+1
    
    Sum1[kminP1:n] <- Sum1[kminP1:n]+ y1[n]
    Aver1[kminP1:n] <- Sum1[kminP1:n]/((n-kmin):1)
    Kvad1[kminP1:n] <- Kvad1[kminP1:n]+ (y1[n])^2
    
    cost1 <- (Kvad1[kminP1:nMkminP1]-Sum1[kminP1:nMkminP1]*Aver1[kminP1:nMkminP1])/sd1^2
    
    Sum2[kminP1:n] <- Sum2[kminP1:n]+ y2[n]
    Aver2[kminP1:n] <- Sum2[kminP1:n]/((n-kmin):1)
    Kvad2[kminP1:n] <- Kvad2[kminP1:n]+ (y2[n])^2
    cost2 <- (Kvad2[kminP1:nMkminP1]-Sum2[kminP1:nMkminP1]*Aver2[kminP1:nMkminP1])/sd2^2
    
    Cost[kminP1:nMkminP1] <- bestCost[kmin:(n-kmin)] + cost1 + cost2
    
    Pos <- which.min(Cost[kminP1:nMkminP1])+kmin
    cost <- Cost[Pos] + gamma
    
    aver1 <- Aver1[Pos]
    aver2 <- Aver2[Pos]
    totAver1 <- (Sum1[kminP1]+initSum1)/n
    totCost1 <- ((Kvad1[kminP1]+initKvad1) - n*totAver1*totAver1)/sd1^2
    totAver2 <- (Sum2[kminP1]+initSum2)/n
    totCost2 <- ((Kvad2[kminP1]+initKvad2) - n*totAver2*totAver2)/sd2^2
    totCost <- totCost1 + totCost2
    
    if (totCost < cost) {
      Pos <- 1
      cost <- totCost
      aver1 <- totAver1
      aver2 <- totAver2
    }
    bestCost[n] <- cost
    bestAver1[n] <- aver1
    bestAver2[n] <- aver2
    bestSplit[n] <- Pos-1
    
    
  }#endfor
  
  
  # Trace back
  n <- N
  breakpts <- n
  while(n > 0){
    breakpts <- c(bestSplit[n], breakpts)
    n <- bestSplit[n]
  }#endwhile
  
  breakpts <- breakpts + from -1
  breakpts <- breakpts[breakpts>=usefrom & breakpts<=useto]
  
  return(list(breakpts=breakpts))
  
}

exactPcf <- function(y, kmin, gamma)
  {
  ## Implementaion of exact PCF by Potts-filtering
  ## x: input array of (log2) copy numbers
  ## kmin: Mininal length of plateaus
  ## gamma: penalty for each discontinuity
  N <- length(y)
  yhat <- rep(0,N);
  if (N < 2*kmin) {
    yhat <- rep(mean(y),N)
    return(yhat)
  }
  initSum <- sum(y[1:kmin])
  initKvad <- sum(y[1:kmin]^2)
  initAve <- initSum/kmin;
  bestCost <- rep(0,N)
  bestCost[kmin] <- initKvad - initSum*initAve
  bestSplit <- rep(0,N)
  bestAver <- rep(0,N)
  bestAver[kmin] <- initAve
  Sum <- rep(0,N)
  Kvad <- rep(0,N)
  Aver <- rep(0,N)
  Cost <- rep(0,N)
  kminP1=kmin+1
  for (k in (kminP1):(2*kmin-1)) {
    Sum[kminP1:k]<-Sum[kminP1:k]+y[k]
    Aver[kminP1:k] <- Sum[kminP1:k]/((k-kmin):1)
    Kvad[kminP1:k] <- Kvad[kminP1:k]+y[k]^2
    bestAver[k] <- (initSum+Sum[kminP1])/k
    bestCost[k] <- (initKvad+Kvad[kminP1])-k*bestAver[k]^2
  }
  for (n in (2*kmin):N) {
    yn <- y[n]
    yn2 <- yn^2
    Sum[kminP1:n] <- Sum[kminP1:n]+yn
    Aver[kminP1:n] <- Sum[kminP1:n]/((n-kmin):1)
    Kvad[kminP1:n] <- Kvad[kminP1:n]+yn2
    nMkminP1=n-kmin+1
    Cost[kminP1:nMkminP1] <- bestCost[kmin:(n-kmin)]+Kvad[kminP1:nMkminP1]-Sum[kminP1:nMkminP1]*Aver[kminP1:nMkminP1]+gamma
    Pos <- which.min(Cost[kminP1:nMkminP1])+kmin
    cost <- Cost[Pos]
    aver <- Aver[Pos]
    totAver <- (Sum[kminP1]+initSum)/n
    totCost <- (Kvad[kminP1]+initKvad) - n*totAver*totAver
    if (totCost < cost) {
      Pos <- 1
      cost <- totCost
      aver <- totAver
    }
    bestCost[n] <- cost
    bestAver[n] <- aver
    bestSplit[n] <- Pos-1
  }
  n <- N
  while (n > 0) {
    yhat[(bestSplit[n]+1):n] <- bestAver[n]
    n <- bestSplit[n]
  }
  return(yhat)
}

fillNA = function(vec, zeroIsNA=TRUE)
  {
  if (zeroIsNA) {vec[vec==0] <- NA}
  nas = which(is.na(vec))
  
  if(length(nas) == 0) {
    return(vec)
  }
  
  # Find stretches of contiguous NAs
  starts = c(1, which(diff(nas)>1)+1)
  ends = c(starts[-1] - 1, length(nas))
  
  starts = nas[starts]
  ends = nas[ends]
  
  # Special-case: vec[1] is NA
  startAt = 1
  if(starts[1]==1) {
    vec[1:ends[1]] = vec[ends[1]+1]
    startAt = 2
  }
  
  if (startAt > length(starts)) {
    return(vec)
  }  
  
  # Special-case: last element in vec is NA
  endAt = length(starts)
  if(is.na(vec[length(vec)])) {
    vec[starts[endAt]:ends[endAt]] = vec[starts[endAt]-1]
    endAt = endAt-1
  }
  
  if (endAt < startAt) {
    return(vec)
  }  
  
  # For each stretch of NAs, set start:midpoint to the value before,
  # and midpoint+1:end to the value after.
  for(i in startAt:endAt) {
    start = starts[i]
    end = ends[i]
    N = 1 + end-start
    if (N==1) {
      vec[start] = vec[start-1]
    } else {
      midpoint = start+ceiling(N/2)
      vec[start:midpoint] = vec[start-1]
      vec[(midpoint+1):end] = vec[end+1]
    }
  }
  
  return(vec)
}

ascat.plotGenotypes<-function(ASCATobj, title, Tumor_BAF_noNA, Hom, ch_noNA)
  {
  par(mar = c(0.5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2, pch = ifelse(dim(ASCATobj$Tumor_LogR)[1]>100000,".",20))
  plot(c(1,length(Tumor_BAF_noNA)), c(0,1), type = "n", xaxt = "n", main = title, xlab = "", ylab = "")
  points(Tumor_BAF_noNA,col=ifelse(is.na(Hom),"green",ifelse(Hom,"blue","red")))
  
  abline(v=0.5,lty=1,col="lightgrey")
  chrk_tot_len = 0
  for (j in 1:length(ch_noNA)) {
    chrk = ch_noNA[[j]];
    chrk_tot_len_prev = chrk_tot_len
    chrk_tot_len = chrk_tot_len + length(chrk)
    vpos = chrk_tot_len;
    tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
    text(tpos,1,ASCATobj$chrs[j], pos = 1, cex = 2)
    abline(v=vpos+0.5,lty=1,col="lightgrey")
  }
}

### Scar scores ####
# LOH, TAI, LST functions 
chrominfoHRD<- function (reference="hg19", chr.in.names=TRUE)
  {
  if (reference == "grch38"){
    chrominfo = chrominfo_grch38
    if(!chr.in.names){
      chrominfo$chrom = gsub("chr","",chrominfo$chr)
      rownames(chrominfo) = chrominfo$chr
    }
  }
   else if(reference == "hg19"){ # Sames as hg19
    chrominfo = chrominfo_hg19()
    if(!chr.in.names){
      chrominfo$chrom = gsub("chr","",chrominfo$chr)
      rownames(chrominfo) = chrominfo$chr
    }
  } 
  else {
    stop()
  }
  return(chrominfo)
}

calc.ai_new<-function(seg, chrominfo, min.size=1e6, cont = 0,ploidyByChromosome=TRUE, shrink=TRUE)
  {
  seg <- seg[seg[,4]- seg[,3] >= min.size,]
  seg <- seg[seg[,9] >= cont,]
  if(shrink){
    seg <- shrink.seg.ai.wrapper(seg)
  }
  AI <- rep(NA, nrow(seg))
  seg <- cbind(seg, AI)
  samples <- as.character(unique(seg[,1]))
  
  ascat.ploidy <- setNames(seg[!duplicated(seg[,1]),8], seg[!duplicated(seg[,1]),1])
  for(j in samples){
    sample.seg <- seg[seg[,1] %in% j,]
    if(!ploidyByChromosome){
      ploidy <- vector()
      for(k in unique(sample.seg[,6])){
        tmp <- sample.seg[sample.seg[,6] %in% k,]
        ploidy <- c(ploidy, setNames(sum(tmp[,4]-tmp[,3]), k))
      }
      ploidy <- as.numeric(names(ploidy[order(ploidy,decreasing=T)]))[1]
      sample.seg[,8] <- ploidy
      if(ploidy %in% c(1,seq(2, 200,by=2))){
        sample.seg[,'AI'] <- c(0,2)[match(sample.seg[,6] == sample.seg[,7], c('TRUE', 'FALSE'))]
      }
      if(!ploidy %in%  c(1,seq(2, 200,by=2))){
        sample.seg[,'AI'] <- c(0,2)[match(sample.seg[,6] + sample.seg[,7] == ploidy & sample.seg[,7] != ploidy, c('TRUE', 'FALSE'))]
      }
    }
    new.sample.seg<- sample.seg[1,]
    for(i in unique(sample.seg[,2])){
      sample.chrom.seg <- sample.seg[sample.seg[,2] %in% i,,drop=F]
      if(nrow(sample.chrom.seg) == 0){ next}
      if(ploidyByChromosome){
        ploidy <- vector()
        for(k in unique(sample.seg[,6])){
          tmp <- sample.chrom.seg[sample.chrom.seg[,6] %in% k,]
          ploidy <- c(ploidy, setNames(sum(tmp[,4]-tmp[,3]), k))
          ploidy <- ploidy[!names(ploidy) %in% 0]
        }
        ploidy <- as.numeric(names(ploidy[order(ploidy,decreasing=T)]))[1]
        sample.chrom.seg[,8] <- ploidy # update "ploidy" column, so the new calculated value can be returned
        if(ploidy %in% c(1,seq(2, 200,by=2))){
          sample.chrom.seg[,'AI'] <- c(0,2)[match(sample.chrom.seg[,6] == sample.chrom.seg[,7], c('TRUE', 'FALSE'))]
        }
        if(!ploidy %in%  c(1,seq(2, 200,by=2))){
          sample.chrom.seg[,'AI'] <- c(0,2)[match(sample.chrom.seg[,6] + sample.chrom.seg[,7] == ploidy & sample.chrom.seg[,7] != 0, c('TRUE', 'FALSE'))]
        }
        sample.seg[sample.seg[,2] %in% i,8] <-ploidy
        sample.seg[sample.seg[,2] %in% i,'AI'] <-sample.chrom.seg[,'AI']
      }
      
      if(class(chrominfo) != 'logical'){# Here we consider the centromere
        if(sample.chrom.seg[1,'AI'] == 2 & nrow(sample.chrom.seg) != 1 & sample.chrom.seg[1,4] < (chrominfo[i,2])){
          sample.seg[sample.seg[,2]==i,'AI'][1] <- 1
        }
        if(sample.chrom.seg[nrow(sample.chrom.seg),'AI'] == 2 & nrow(sample.chrom.seg) != 1 & sample.chrom.seg[nrow(sample.chrom.seg),3] > (chrominfo[i,3])){
          sample.seg[sample.seg[,2]==i,'AI'][nrow(sample.seg[sample.seg[,2]==i,])] <- 1
        }
      }
      if(nrow(sample.seg[sample.seg[,2]==i,]) == 1 & sample.seg[sample.seg[,2]==i,'AI'][1] != 0){
        sample.seg[sample.seg[,2]==i,'AI'][1] <- 3
      }
    }
    seg[seg[,1] %in% j,] <- sample.seg
  }
  samples <- as.character(unique(seg[,1]))
  #0 = no AI, 1=telomeric AI, 2=interstitial AI, 3= whole chromosome AI
  no.events <- matrix(0, nrow=length(samples), ncol=12)
  rownames(no.events) <- samples
  colnames(no.events) <- c("Telomeric AI", "Mean size", "Interstitial AI", "Mean Size", "Whole chr AI", "Telomeric LOH",  "Mean size", "Interstitial LOH", "Mean Size", "Whole chr LOH", "Ploidy", "Aberrant cell fraction")
  for(j in samples){
    sample.seg <- seg[seg[,1] %in% j,]
    no.events[j,1] <- nrow(sample.seg[sample.seg[,'AI'] == 1,])
    no.events[j,2] <- mean(sample.seg[sample.seg[,'AI'] == 1,4] - sample.seg[sample.seg[,'AI'] == 1,3])
    no.events[j,3] <- nrow(sample.seg[sample.seg[,'AI'] == 2,])
    no.events[j,4] <- mean(sample.seg[sample.seg[,'AI'] == 2,4] - sample.seg[sample.seg[,'AI'] == 2,3])
    no.events[j,5] <- nrow(sample.seg[sample.seg[,'AI'] == 3,])
    no.events[j,11] <- ascat.ploidy[j]
    no.events[j,12] <- unique(sample.seg[,9]) # aberrant cell fraction
    #Here we restrict ourselves to real LOH
    sample.seg <- sample.seg[sample.seg[,7] == 0,]
    no.events[j,6] <- nrow(sample.seg[sample.seg[,'AI'] == 1,])
    no.events[j,7] <- mean(sample.seg[sample.seg[,'AI'] == 1,4] - sample.seg[sample.seg[,'AI'] == 1,3])
    no.events[j,8] <- nrow(sample.seg[sample.seg[,'AI'] == 2,])
    no.events[j,9] <- mean(sample.seg[sample.seg[,'AI'] == 2,4] - sample.seg[sample.seg[,'AI'] == 2,3])
    no.events[j,10] <- nrow(sample.seg[sample.seg[,'AI'] == 3,])
  }
  return(no.events)
}

calc.hrd<-function(seg, nA=6, return.loc=F,sizelimit1)
  {
  nB <- nA+1
  output <- rep(0, length(unique(seg[,1])))
  names(output) <- unique(seg[,1])
  if(return.loc) {
    out.seg <- matrix(0,0,9)
    colnames(out.seg) <- c(colnames(seg)[1:8],'HRD breakpoint')
  }
  #For multiple patients
  for(i in unique(seg[,1])){
    segSamp <- seg[seg[,1] %in% i,]
    chrDel <-vector()
    for(j in unique(segSamp[,2])){
      if(all(segSamp[segSamp[,2] == j,nB] == 0)) {
        chrDel <- c(chrDel, j)
      }
    }
    segSamp[segSamp[,nA] > 1,nA] <- 1
    segSamp <- shrink.seg.ai.wrapper(segSamp)
    segLOH <- segSamp[segSamp[,nB] == 0 & segSamp[,nA] != 0,,drop=F]
    segLOH <- segLOH[segLOH[,4]-segLOH[,3] > sizelimit1,,drop=F]
    segLOH <- segLOH[!segLOH[,2] %in% chrDel,,drop=F]
    output[i] <- nrow(segLOH)
    if(return.loc){
      if(nrow(segLOH) < 1){next}
      segLOH <- cbind(segLOH[,1:7], 1)
      colnames(segLOH)[8] <- 'HRD breakpoint'
      out.seg <- rbind(out.seg, segLOH)
    }
  }
  
  if(return.loc){
    return(out.seg)
  } else {
    return(output)
  }
  print(paste("return.loc", return.loc))
  print(paste("output", output))
}

calc.lst<-function(seg, chrominfo=chrominfo,nA=6,chr.arm='no')
  {
  nB <- nA+1
  samples <- unique(seg[,1])
  output <- setNames(rep(0,length(samples)), samples)
  for(j in samples){
    sample.seg <- seg[seg[,1] %in% j,]
    sample.lst <- c()
    chroms <- unique(sample.seg[,2])
    chroms <- chroms[!chroms %in% c(23,24,'chr23','chr24','chrX','chrx','chrY','chry')]
    for(i in chroms){
      sample.chrom.seg <- sample.seg[sample.seg[,2] %in% i,,drop=F]
      if(chr.arm !='no'){
        p.max <- if(any(sample.chrom.seg[,chr.arm] == 'p')){max(sample.chrom.seg[sample.chrom.seg[,chr.arm] == 'p',4])}
        q.min <- min(sample.chrom.seg[sample.chrom.seg[,chr.arm] == 'q',3])
      }
      if(nrow(sample.chrom.seg) < 2) {next}
      sample.chrom.seg.new <- sample.chrom.seg
      if(chr.arm == 'no'){
        p.arm <- sample.chrom.seg.new[sample.chrom.seg.new[,3] <= chrominfo[i,2],,drop=F] # split into chromosome arms
        q.arm <- sample.chrom.seg.new[sample.chrom.seg.new[,4] >= chrominfo[i,3],,drop=F]
        q.arm<- shrink.seg.ai(q.arm)
        q.arm[1,3] <- chrominfo[i,3]
        if(nrow(p.arm) > 0){
          p.arm<- shrink.seg.ai(p.arm)
          p.arm[nrow(p.arm),4] <- chrominfo[i,2]
        }
      }
      if(chr.arm != 'no'){
        q.arm <- sample.chrom.seg.new[sample.chrom.seg.new[,chr.arm] == 'q',,drop=F]
        q.arm<- shrink.seg.ai(q.arm)
        q.arm[1,3] <- q.min
        if(any(sample.chrom.seg.new[,chr.arm] == 'p')){
          p.arm <- sample.chrom.seg.new[sample.chrom.seg.new[,chr.arm] == 'p',,drop=F] # split into chromosome arms
          p.arm<- shrink.seg.ai(p.arm)
          p.arm[nrow(p.arm),4] <- p.max
        }
      }
      n.3mb <- which((p.arm[,4] - p.arm[,3]) < 3e6)
      while(length(n.3mb) > 0){
        p.arm <- p.arm[-(n.3mb[1]),]
        p.arm <- shrink.seg.ai(p.arm)
        n.3mb <- which((p.arm[,4] - p.arm[,3]) < 3e6)
      }
      if(nrow(p.arm) >= 2){
        p.arm <- cbind(p.arm[,1:7], c(0,1)[match((p.arm[,4]-p.arm[,3]) >= 10e6, c('FALSE','TRUE'))])
        for(k in 2:nrow(p.arm)){
          if(p.arm[k,8] == 1 & p.arm[(k-1),8]==1 & (p.arm[k,3]-p.arm[(k-1),4]) < 3e6){
            sample.lst <- c(sample.lst, 1)
          }
        }
      }
      n.3mb <- which((q.arm[,4] - q.arm[,3]) < 3e6)
      while(length(n.3mb) > 0){
        q.arm <- q.arm[-(n.3mb[1]),]
        q.arm <- shrink.seg.ai(q.arm)
        n.3mb <- which((q.arm[,4] - q.arm[,3]) < 3e6)
      }
      if(nrow(q.arm) >= 2){
        q.arm <- cbind(q.arm[,1:7], c(0,1)[match((q.arm[,4]-q.arm[,3]) >= 10e6, c('FALSE','TRUE'))])
        for(k in 2:nrow(q.arm)){
          if(q.arm[k,8] == 1 & q.arm[(k-1),8]==1 & (q.arm[k,3]-q.arm[(k-1),4]) < 3e6){
            sample.lst <- c(sample.lst, 1)
            
          }
        }
      }
    }
    output[j] <- sum(sample.lst)
  }
  return(output)
}

#Functional
shrink.seg.ai<-function(chr.seg)
  {
  new.chr <- chr.seg
  if(nrow(chr.seg) > 1){
    new.chr <- matrix(0,0,ncol(chr.seg))
    colnames(new.chr) <- colnames(chr.seg)
    new.chr <- chr.seg
    seg.class <- c(1)
    for(j in 2:nrow(new.chr)){
      seg_test <- new.chr[(j-1),6] == new.chr[j,6] & new.chr[(j-1),7] == new.chr[j,7]
      if(seg_test){
        seg.class <- c(seg.class, seg.class[j-1])
      }
      if(!seg_test){
        seg.class <- c(seg.class, seg.class[j-1]+1)
      }
    }
    for(j in unique(seg.class)){
      new.chr[seg.class %in% j,4] <- max(new.chr[seg.class %in% j,4])
      new.chr[seg.class %in% j,5] <- sum(new.chr[seg.class %in% j,5])
    }
    new.chr<- new.chr[!duplicated(seg.class),]
  }
  if(nrow(chr.seg) == 1){
    new.chr <- chr.seg
  }
  return(new.chr)
}

shrink.seg.ai.wrapper<-function(seg)
  {
  new.seg <- seg[1,]
  for(j in unique(seg[,1])){
    sample.seg <- seg[seg[,1] %in% j,]
    new.sample.seg <- seg[1,]
    for(i in unique(sample.seg[,2])){
      sample.chrom.seg <- sample.seg[sample.seg[,2] %in% i,,drop=F]
      if(nrow(sample.chrom.seg) > 1){
        sample.chrom.seg <- shrink.seg.ai(sample.chrom.seg)
      }
      new.sample.seg <- rbind(new.sample.seg, sample.chrom.seg)
    }
    new.seg <- rbind(new.seg, new.sample.seg[-1,])
  }
  seg <- new.seg[-1,]
  return(seg)
}

#Fixed for columns
preprocess.hrd<-function(seg)
  {
  seg <- seg[!seg[,2] %in% c(paste('chr',c('X','Y','x','y',23,24),sep=''),c('X','Y','x','y',23,24)),]
  seg[,1] <- as.character(seg[,1])
  
  if(! all(seg[,7] <= seg[,6]) ){
    tmp <- seg
    seg[tmp[,7] > tmp[,6],6]  <- tmp[tmp[,7] > tmp[,6],7]
    seg[tmp[,7] > tmp[,6],7]  <- tmp[tmp[,7] > tmp[,6],6]
  }
  seg <- shrink.seg.ai.wrapper(seg)
  return(seg)
  
}  

chrominfo_hg19 <- function()
  {
  # Get chromInfo table from UCSC
  chrom <- GetGzFromUrl("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz", header = FALSE)
  chrom <- subset(chrom, grepl("^chr[0-9XY]{1,2}$", chrom[,1]))
  # Get gap table from UCSC
  gaps <- GetGzFromUrl("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz", header = FALSE)
  centro <- subset(gaps, gaps[,8] == "centromere")
  # Merge the relevant info from the two tables
  chrominfo <- merge(chrom[,1:2], centro[,2:4], by.x = 1, by.y = 1) # merge size and centromere location
  chrominfo$centromere <- rowMeans(chrominfo[,3:4])# convert centromere start and end into one location (the mean)
  chrominfo <- chrominfo[,c(1,2,5,3,4)] # keep chromosome, size and centromere location
  colnames(chrominfo) <- c("chr", "size", "centromere", "centstart", "centend")
  #chrominfo[,1] <- as.character(chrominfo[,1])
  #chrominfo$chr <- sub("chr", "", chrominfo$chr)
  #chrominfo$chr <- sub("X", "23", chrominfo$chr)
  #chrominfo$chr <- sub("Y", "24", chrominfo$chr)
  #chrominfo[,1] <- as.numeric(chrominfo[,1])
  #chrominfo <- chrominfo[order(chrominfo$chr), ]  # order by chromosome number
  rownames(chrominfo) <- as.character(chrominfo[,1])
  chrominfo <- as.data.frame(chrominfo[,c(1,4,5)])
  return(chrominfo)  
}

GetGzFromUrl <- function(url, ...)
  {
  # http://stackoverflow.com/questions/9548630/read-gzipped-csv-directly-from-a-url-in-r
  con <- gzcon(url(url))
  txt <- readLines(con)
  dat <- read.delim(textConnection(txt), ...)  
  return(dat)
}