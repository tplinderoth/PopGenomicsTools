#!/usr/bin/env Rscript

# betaAFOutlier.R
# version 1.2
# This is the script to update and maintain

### LIBRARIES
library(rmutil, quietly=TRUE, warn.conflicts=FALSE) # betabinom

### FUNCTIONS

parseArgs <- function(args) {

	userin = list(dipn=NA, fstbar=NA, obsfreq=NA, outnames=NA, ancmethod=0, seed=NULL, qq=NULL, ancbin=0)

	# parse required input
	userin$dipn = as.numeric(strsplit(args[1],',')[[1]])
	if (userin$dipn[1] <= 0) stop("Pop1 size must be > 0")
	if (userin$dipn[2] <= 0) stop("Pop2 size must be > 0")
	userin$dipn = c(userin$dipn, sum(userin$dipn))

	userin$fstbar = as.numeric(args[2])
	if (userin$fstbar != 9 & (userin$fstbar < 0 || userin$fstbar > 1)) stop("Geome-wide FST out of range [0,1]")

	userin$obsfreq <- read.table(args[3], head=TRUE)
	colnames(userin$obsfreq) = c("chr", "pos", "pop1f", "pop2f")
	nsites = nrow(userin$obsfreq)
	cat(paste0("Observed frequencies file contains ", nsites, " sites\n"))
	# Discard nonvariable sites
	fixedidx = which((userin$obsfreq$pop1f == 0 & userin$obsfreq$pop2f == 0) | (userin$obsfreq$pop1f == 1 & userin$obsfreq$pop2f == 1))
	nfix = length(fixedidx)
	if (nfix > 0) {
		userin$obsfreq = userin$obsfreq[-fixedidx,]
		cat("Pruned ", nfix, "novariable sites\n")
		if (!nrow(userin$obsfreq)) stop("No input sites are variable")
	}

	suffix=c(".dist", ".ancp", ".sigtest", ".pdf", "_qqplot.png")
	userin$outnames = paste0(args[4],suffix)

	# parse optional input
	nargs = length(args)
	i=5
	while (i <= nargs) {
		if (args[i] == "--ancmethod") {
			if (i+1 <= nargs) userin$ancmethod = as.numeric(args[i+1]) else stop("--ancmethod requires a 0/1 argument")
			if (userin$ancmethod) {
				if (userin$ancmethod != 1) stop(paste0("Invalid ancestral allele frequency method: ",userin$ancmethod))
			}
		} else if (args[i] == "--seed") {
			if (i+1 <= nargs) userin$seed = as.numeric(args[i+1]) else stop("--seed requires an integer argument")
		} else if (args[i] == "--plotqq") {
			userin$qq = 1
			i = i-1
		} else if (args[i] == "--ancbin") {
			if (i+1 <= nargs) userin$ancbin = as.numeric(args[i+1]) else stop("--ancbin requires a value in [0,1]")
			if (userin$ancbin < 0 || userin$ancbin > 1) stop("--ancbin argument out of range [0,1]")
		} else {
			stop(paste("Unknown argument,", args[i]))
		}
		i = i+2
	}

	return(userin)
}

expected_af_prob <- function(ndip, binwidth=0) {
	# Calculates the expected probability of randomly sampling a site
	# from binned MAF categories. Based on neutral 1/x prior where x is 
	# the allele frequency.
 
	# ndip: diploid sample size
	# binwidth: size of allele frequency bins, 0 means no binning

	# unfolded SFS probabilities
	n = 2*ndip
	treelen = 0
	sfs = rep(NA,n-1)
	for (k in 1:(n-1)) {
		sfs[k] = 1/k
		treelen = treelen + sfs[k]
	}
	sfs = sfs/treelen
	
	# folded sfs (MAF) probabilities
	sfs.fold = sapply(1:ndip, function(x,sfs,n){ifelse(x<ndip, sfs[x]+sfs[n-x], sfs[x])}, sfs=sfs, n=n)
	
	step = ifelse(binwidth, binwidth, 1/n)
	probs = data.frame(lower=seq(from=0,to=(0.5-step),by=step), upper=seq(from=step,to=0.5,by=step), prob=0)
	if (binwidth) {
		# this is for using binned ancestral probabilities
		j=1
		for(i in 1:length(sfs.fold)) {
			if (i/n > probs$upper[j]) j = j+1
			probs$prob[j] = probs$prob[j] + sfs.fold[i]
		}
	} else probs$prob = sfs.fold

	return(probs)
}

empiric_af_prob <- function(af, binwidth=0, pop1n, pop2n, weight=TRUE) {
	# Calculates the probability of randomly sampling a site
	# from binned MAF category. Ancestral allele frequency
	# considered to be average among descedant pops.

	# af: a matrix where each row is a site and columns are allele frequencies
	# binwidth: size of allele frequency bins, 0 means no binning
	# pop1n: diploid size for pop1
	# pop2n: diploid size for pop2
	# weight descedant population allele frequencies by sample size
	
	ndip = pop1n + pop2n
	if (weight) avgf = (pop1n*af[,1] + pop2n*af[,2])/ndip else avgf = (af[,1] + af[,2])/2
	avgf[avgf > 0.5] = 1-avgf[avgf > 0.5] # convert to MAF
	if (!binwidth) {
		if (is.null(ndip)) stop("Unbinned, empiric ancestral prior requires total population sample size")
		binwidth = 1/(2*ndip)
	}
	probs = data.frame(lower=seq(from=0,to=(0.5-binwidth),by=binwidth), upper=seq(from=binwidth,to=0.5,by=binwidth), prob=NA)
	probs$prob=unname(table(cut(avgf, breaks=seq(from=0,to=0.5,by=binwidth),include.lowest=FALSE, right=TRUE)))
	probs$prob = probs$prob/sum(probs$prob)

	return(probs)
}

empiric_sfs_prob <- function(af, pop1n, pop2n, fold=TRUE, bins=TRUE) {
	# af: 2-column matrix where each row is a site and columns are allele frequencies
	# pop1n: population 1 diploid size
	# pop2n: population 2 diploid size
	# fold: TRUE to fold the spectrum, false for unfolded spectrum
	# bins: TRUE to return low & upper frequency bounds, FALSE to return unbinned frequencies

	sfs = calcAncSFS(af, pop1n, pop2n, fold)
	nchr = 2*(pop1n+pop2n)
	binwidth = 1/nchr
	if (fold) {
		sfs = sfs[-1] # remove fixed counts
		ub = 0.5
	} else {
		sfs = sfs[-c(1,length(sfs))]
		ub = (nchr-1)/nchr
	}

	if (bins) {
		return(data.frame(lower=seq(from=0,to=(ub-binwidth),by=binwidth),upper=seq(from=binwidth, to=ub, by=binwidth), prob=sfs/sum(sfs)))
	} else return(data.frame(freq=seq(from=binwidth,to=ub,by=binwidth), prob=sfs/sum(sfs)))
}


BaldingNicholsSamp <- function(nsamp, fst, ancf, samp=0) {
	# Samples pairs of allele frequencies for each ancestral allele
	# frequency bin and a given FST under the Balding-Nichols Beta
	# distribution model.

	# n: number allele frequency pairs to draw for each ancestral allele frequency bin
	# fst: genome-wide background FST between populations
	# ancf: data.frame of ancestral allele frequency bin probabilities
	# samp: 0 to set p0 to upper bin bound, 1 to uniformly sample p0 from bin

	naf = nrow(ancf)
	gensamp = list()
	for (i in 1:naf) {
		lb = ancf$lower[i]
		ub = ancf$upper[i]
		cat(paste0("Taking ", nsamp, "frequency pairs for ancestral p0 (", lb, ",", ub, "]\n"))
		if (samp) {
			p0 = runif(nsamp, min=lb, max=ub)
			while (length(which(p0==lb))) {
				p0[p0==lb] = runif(length(which(p0==lb)), min=lb, max=ub)
			}
		} else p0 = rep(ub, nsamp) # using bin upper bound because ancestral freq could never be less than singleton
		alpha = (1-fst)/fst * p0
		beta = (1-fst)/fst * (1-p0)
		gensamp[[i]] = matrix(rbeta(2*nsamp, shape1=alpha, shape2=beta), nrow=nsamp, ncol=2)
	}
	return(gensamp)
}

sampleCountsFreq <- function(af, pop1n, pop2n) {
	# Generates allele frequencies for SNPs in a sample
	# based on binomial sampling of alleles.

	# af: list of nsampx2 matrices containing allele frequency pairs for various ancestral frequency bins
	# pop1n: population 1 diploid sample size
	# pop2n: population 2 diploid sample size

	n1 = 2*pop1n
	n2 = 2*pop2n

	countsamp = list()
	for (i in 1:length(af)) {
		cat(paste0("Sampling from frequency matrix ",i,"\n"))
		n = nrow(af[[i]])
		countsamp[[i]] = matrix(c(rbinom(n, size=n1, prob=af[[i]][,1])/n1, rbinom(n, size=n2, prob=af[[i]][,2])/n2), nrow=n, ncol=2)
		# retain only variable sites
		fixed = which((countsamp[[i]][,1] == 0 & countsamp[[i]][,2] == 0) | (countsamp[[i]][,1] == 1 & countsamp[[i]][,2] == 1))
		nfixed = length(fixed)
		if (nfixed > 0) {
			countsamp[[i]] = countsamp[[i]][-fixed,]
			cat("Excluded", nfixed, "nonvariable sites\n")
		}
	}

	return(countsamp)
}

sampleFreq <- function(nsamp, fst, pop1n, pop2n, ancf) {
	# Generates allele frequencies at SNPs among two populations by first
	# drawing allele frequencies under the Balding-Nichols model 
	# based on binomial sampling of alleles.

	# nsamp: number allele frequency pairs to draw for each ancestral allele frequency bin
	# fst: genome-wide background FST between populations
	# ancf: 2-column matrix of ancestral allele frequency bin (1) upper and (2) lower bounds
	# pop1n: population 1 diploid sample size
	# pop2n: population 2 diploid sample size

	rBaldingNichols <- function(nsamp, fst, freq) {
		# freq: A vector of length 2 indicates lower and upper frequecny bounds, otherwise a scalar gives a single frequency

		if (freq[1] > freq[2]) stop("Ancestral frequency lower bound greater than upper bound")
		p0 = runif(nsamp, min=freq[1], max=freq[2])
		while (length(which(p0==freq[1]))) {
			p0[p0==freq[1]] = runif(length(which(p0==freq[1])), min=freq[1], max=freq[2]) # ensure half open intervals
		}
		alpha = (1-fst)/fst * p0
                beta = (1-fst)/fst * (1-p0)
                matrix(rbeta(2*nsamp, shape1=alpha, shape2=beta), nrow=nsamp, ncol=2)
	}


        n1 = 2*pop1n
        n2 = 2*pop2n
	nsets = nrow(ancf)
	afsamp = list()
	idx = 1:nsamp

        for (i in 1:nsets) {
		cat(paste0("Sampling frequency pairs for ancestral p0 (", ancf[i,1], ",", ancf[i,2], "]\n"))
		afsamp[[i]] = matrix(0, nrow=nsamp, ncol=2)
		fixed = idx
		nfill = nsamp
		while (nfill) {
			cat("Filling in", nfill, "draws\n")
			af = rBaldingNichols(nsamp, fst, freq=as.vector(ancf[i,])) # perform genetic sampling under Balding-Nichols model
			countsamp = matrix(c(rbinom(nsamp, size=n1, prob=af[,1]), rbinom(nsamp, size=n2, prob=af[,2])), nrow=nsamp, ncol=2)
			sampfixed = which((countsamp[,1] == 0 & countsamp[,2] == 0) | (countsamp[,1] == n1 & countsamp[,2] == n2))
			countsamp[,1] = countsamp[,1]/n1
			countsamp[,2] = countsamp[,2]/n2
			nsampfixed = length(sampfixed)
			if (nsampfixed < nsamp) {
				if (nsampfixed == 0) segidx = idx else segidx = idx[-sampfixed]
				nseg = length(segidx)
				nreplace = ifelse(nfill > nseg, nseg, nfill)
				afsamp[[i]][fixed[1:nreplace],] = countsamp[segidx[1:nreplace],]
				fixed = which(afsamp[[i]][,1] == 0 & afsamp[[i]][,2] == 0)
				nfill = length(fixed)
			}
		}
        }

        return(afsamp)
}

afDiffProbs <- function(af, diffwidth, ancwidth) {
	# Calculates the probability of observing binned allele frequency differences
	# in a sample given different ancestral population allele frequencies.

	# af: list of nsampx2 matrices containing sample allele frequencies for different ancestral frequency bins
	# diffwidth: width of absolute allele freqeuncy difference bins
	# ancwidth: width of ancestral allele frequency bins

	# initialize matrix with lower and upper bounds for allele frequency bins
	ancf_bins = paste0("(",seq(from=0, to=(0.5-ancwidth), by=ancwidth),",",seq(from=ancwidth, to=0.5, by=ancwidth),"]") # ancestral frequency bin names
	diffp = data.frame(lower=seq(from=0, to=(1-diffwidth), by=diffwidth), upper=seq(from=diffwidth, to=1, by=diffwidth))
	diffp = cbind(diffp, matrix(nrow=nrow(diffp), ncol=length(af)))
	colnames(diffp) = c("lower", "upper", ancf_bins)
	
	# insert bin probabilities
	for (i in 1:length(af)) {
		j = i+2
		if (nrow(af[[i]]) > 0) {
			diffp[,j] = table(cut(abs(af[[i]][,1] - af[[i]][,2]),breaks=seq(from=0,to=1,by=diffwidth),include.lowest=TRUE, right=TRUE))
			diffp[,j] = diffp[,j]/sum(diffp[,j])
		} else diffp[,j] = 0
	}

	return(diffp)
}

weightedDiffProbs <- function(afprobs, ancprobs) {
	# Calculates probability of observing binned allele frequency differences
	# in a sample by summing over the possible ancestral frequencies

	# afprobs: data.frame of lower & upper bin bounds for allele frequency differences and their associated probabilities
	# under different ancestral allele frequency bins
	# ancprobs: data.frame containing the probabilities of sampling sites in different ancestral frequency bins

	distdf = data.frame(lower=afprobs$lower, upper=afprobs$upper, prob=NA, cdf=NA)
	distdf$prob = sapply(1:nrow(afprobs), function(i, probs, prior){sum(probs[i,]*prior)}, probs=afprobs[,3:ncol(afprobs)], prior=ancprobs$prob)
	distdf$cdf[1] = distdf$prob[1]
	for (i in 2:nrow(distdf)) distdf$cdf[i] = distdf$cdf[i-1]+distdf$prob[i]
	
	return(distdf)
}

testObserved <- function(obs, null) {
	# Calculate p and q values for observed allele frequency differences

	# obs: data.frame with columns (1) chr, (2) pos, (3) pop1 allele frequency, (4) pop2 allele frequency
	# null: data.frame specifying the PMS and CDF for allele frequencyc difference bins

	# format bin cutoffs to avoid precision issues
	null$lower = as.numeric(sprintf(null$lower,fmt='%#.3f'))
	null$upper = as.numeric(sprintf(null$upper,fmt='%#.3f'))

	# calculate p-values
	pval = 1-null$cdf
	obs$diff = abs(obs[,3]-obs[,4])
	obs$pval = NA

	for (chr in unique(obs[,1])) {
		cat(paste0("Finding p-values for ",chr,"\n"))
		chridx = which(obs[,1] == chr)
		afdiff = obs$diff[chridx]
		afdiff[afdiff == 0] = 1e-16 # recode allele frequency differences of zero so that they are included in first bin
		p <- rep(NA,length(afdiff))
		for (i in 1:nrow(null)) p[which(afdiff > null$lower[i] & afdiff <= null$upper[i])] = pval[i]
		obs$pval[chridx] = p
	}
	obs$pval[obs$pval < 0] = 0 # prevent negative values due to imprecision

	# calculate q-values
	obs$qval = p.adjust(obs$pval, method="fdr") # Benjamini & Hochberg adjustment

	return(obs)
}

plotFit <- function(obs, null) {
	# Generates plots related to how observed compare to expected allele frequency differences

	# obs: data.frame containing observed allele frequency differences
	# null: data.frame specifying the probabilities and CDF for allele frequency difference bins

	# calculate observed bin probabilities
	obsbins = null[,1:2]
	obsbins$prob = unname(table(cut(obs$diff,breaks=seq(from=0,to=1,by=(null$upper[1]-null$lower[1])),include.lowest=TRUE, right=TRUE)))
	obsbins$prob = as.numeric(obsbins$prob/sum(obsbins$prob))

	# plot observed vs expected bin probabilities
	ymax = max(c(null$prob,obsbins$prob))
	diffmid = (obsbins[,1]+obsbins[,2])/2 # vector of frequency bin midpoints
	plot(x=diffmid, y=obsbins$prob, type="l",col="green",lwd=2,xlab="Allele frequency difference", ylab="Frequency",ylim=c(0,ymax),xlim=c(0,1),main="Allele frequency difference distribution")
	lines(x=diffmid, y=null$prob, type="l",lwd=2,col="blue")
	legend('topright', legend=c("observed", "expected"), col=c("green","blue"), lty=1, bty='n',lwd=2)

	# zoom in on distribution upper tail
	zoomquant = quantile(obs$diff,0.95) # allele frequency difference 0.95 quantile
	binidx = which(obsbins$lower < zoomquant & obsbins$upper >= zoomquant)
	zoom.ymax = max(c(null$prob[binidx], obsbins$prob[binidx]))
	plot(x=diffmid, y=obsbins$prob, type="l", col="green",lwd=2,xlab="Allele frequency difference", ylab="Frequency",ylim=c(0,zoom.ymax),xlim=c(zoomquant,1), main="Observed 95 percentile zoom")
	lines(x=diffmid, y=null$prob, type="l",lwd=2,col="blue")
	legend('topright', legend=c("observed", "expected"), col=c("green","blue"), lty=1, bty='n',lwd=2)

	# plot section of distribution with FDR of being an outlier < 0.05
	sigidx = which(obs$qval < 0.05)
	if (length(sigidx) > 0) {
		mindiff = min(obs$diff[sigidx])
		binidx = which(obsbins$lower < mindiff & obsbins$upper >= mindiff)
		zoom.ymax = max(c(null$prob[binidx], obsbins$prob[binidx]))
		zoom.xmin = obsbins$lower[binidx]
		plot(x=diffmid, y=obsbins$prob, col="green",xlab="Allele frequency difference", ylab="Frequency",ylim=c(0,zoom.ymax),xlim=c(zoom.xmin,1), main="FDR < 0.05 zoom")
		lines(x=diffmid, y=obsbins$prob,col="green",lwd=2)
		points(x=diffmid, y=null$prob,col="blue")
		lines(x=diffmid, y=null$prob,lwd=2,col="blue")
		legend('topright', legend=c("observed", "expected"), col=c("green","blue"), lty=1, bty='n',lwd=2)
	}

}

pvalqq <- function(pobs) {
	# makes a qq-plot of observed vs expected p-values
	# pobs: vector of observed p-values

	zeroidx = which(pobs <= 0)
	if (length(zeroidx) > 0) {
		pobs = pobs[-zeroidx]
		cat(length(zeroidx),"sites with p-value 0 excluded from qq-plot\n")
		if (length(pobs) < 1) {
			cat("Skipping qq-plot because no sites with nonzero p-value\n")
			return(0)
		}
	}
	p.expect <- sort(-log10(runif(length(pobs)))) # generate expected p-values
	plot(p.expect, sort(-log10(pobs)), xlab='-log10(unif[0,1])', ylab='-log10(observed p-values)', main='p-value qqplot')
	abline(0,1, col="red",lty=2)
}

WCFst <- function(f, n1, n2) {
	# f: marix of allele frequencies in [,1] pop1 and [,2] pop2
	# n1: pop1 diploid sample size
	# n2: pop2 diploid sample size

	reynoldsVar <- function(f, n1, n2) {
		npool = n1+n2
		fpool = n1/npool*f[1] + n2/npool*f[2]
		alpha1 = 2*f[1]*(1-f[1])
		alpha2 = 2*f[2]*(1-f[2])
		b = (n1*alpha1 + n2*alpha2)/(npool-1)
		a = (4*n1*(f[1]-fpool)^2 + 4*n2*(f[2]-fpool)^2 - b)/(4*n1*n2/npool)
		return(c(a,b))
	}

	varcomp = t(apply(f,1,reynoldsVar,n1=n1,n2=n2))
	varcomp[,2] = varcomp[,2]+varcomp[,1]
	return(varcomp)
}

wrightFst <- function(f) {
	# f: vector of allele frequencies in [,1] pop1 and [,2] pop2

	# FST = (Ht - Hs)/Ht
	# Assuming both subpops and total pop in HWE & equal subpop sizes

	ftotal = (f[,1]+f[,2])/2 # average allele frequency in pooled population
	ht = 2*ftotal*(1-ftotal)
	hs = f[,1]*(1-f[,1]) + f[,2]*(1-f[,2])
	
	# remove sites with zero heterozygosity in total population
	zerohet = which(ht == 0)
	if (length(zerohet) >0) {
		ht = ht[-zerohet]
		hs = hs[-zerohet]
	}

	if (length(ht)>0) (ht-hs)/ht else return(NA)
}

genomeFst <- function(f, n1, n2) {
	# f: marix of allele frequencies in [,1] pop1 and [,2] pop2
	# n1: pop1 diploid sample size
	# n2: pop2 diploid sample size
	variances = WCFst(f, n1, n2)
	sum(variances[,1])/sum(variances[,2])
}

sampleJointBins <- function(p, n) {
	# p: vector of probabilities
	# n: number of bins to sample

	counts = matrix(nrow=length(p),ncol=3)
	counts[1,] = c(0,p[1],0)
	for (i in 2:nrow(counts)) counts[i,] = c(counts[i-1,2], counts[i-1,2]+p[i],0)
	
	x = runif(n)	
	while (length(which(x==0))) {
			x[x==0] = runif(length(which(x==0)))
	}

	for (i in 1:length(x)) {
		idx = which(x[i] > counts[,1] & x[i] <= counts[,2])
		counts[idx,3] = counts[idx,3] + 1	
	}

	return(counts[,3])
}

calcSFS <- function(f, n, fold=TRUE) {
	# f: vector of allele frequencies
	# n: diploid sample size
	# fold: FALSE=unfolded SFS, TRUE=folded SFS
	
	nhap = 2*n
	sfs = unname(table(cut(round(f*nhap), breaks=seq(from=0,to=nhap,by=1),include.lowest=FALSE, right=TRUE)))
	sfs = c(length(f)-sum(sfs),sfs)
	if (fold) sfs = sapply(1:(n+1), function(x,sfs,n){ifelse(x <= nhap/2, sfs[x]+sfs[n-x+2], sfs[x])}, sfs=sfs, n=nhap)
	names(sfs) = 0:(length(sfs)-1)
	return(sfs)
}

calcAncSFS <- function(fmat, n1, n2, fold=TRUE) {
	# fmat: 2-column matrix of allele frequencies for [,1] pop1 and [,2] pop2. Rows are sites.
	# n1: population 1 diploid sample size
	# n2: population 2 diploid sample size
	# fold: FALSE = unfolded SFS, TRUE = folded SFS

	ntotal = n1+n2
	p0 = (n1*fmat[,1] + n2*fmat[,2])/ntotal
	sfs = calcSFS(p0, ntotal, fold=fold)
	
}

simExpectedSample <- function(obsf, n1, n2, ancf, nsamp, afsamp=NULL) {
	# obsf: 2-column matrix of observed allele frequencies in [,1] pop1 and [,2] pop2. Rows are sites.
	# n1: population1 diploid sample size
	# n2: population2 diploid sample size
	# ancf: data.frame of ancestral allele frequency bin probabilities
	# nsamp: number of allele frequency pairs to generate for each ancestral frequency

	nsites = nrow(obsf)
	if (is.null(afsamp)) {
		obsvar = WCFst(obsf, n1, n2)
		obsfst.avg = sum(obsvar[,1])/sum(obsvar[,2])
		afsamp = sampleFreq(nsamp, obsfst.avg, n1, n2, as.matrix(ancf[,c(1,2)]))
	}
	bincounts = sampleJointBins(ancf$prob, n=nsites)
	expect.freq = NULL
	while (is.null(expect.freq) || nrow(expect.freq) < nsites) {
		for (i in 1:length(afsamp)) expect.freq = rbind(expect.freq, afsamp[[i]][sample(1:nsamp, bincounts[i]),])
	}
	return(expect.freq)
}

betaBinomDist <- function(n1, n2, fst, p, binsize = NULL) {
	# calculates the probability for all possible allele frequency differences
	# n1: population 1 diploid sample size
	# n2: population 2 diploid sample size
	# fst: average FST between pops 1 and 2
	# p: 2-column matrix of [,1] allele frequency and [,2] probability of allele frequency
	# binsize: bin width for allele frequency differences

# THINK I NEED TO USE UNFOLDED W.R.T. REF ALLELE

	bbprob <- function(k, n, fst, p0, p0.prob) {
		if (sum(p0.prob) != 1) stop("p0 probabilities must sum to 1")
		a = (1-fst)/fst * p0 # standard alpha (shape1) for beta
		b = (1-fst)/fst * (1-p0) # standard beta (shape2) for beta
		sum(dbetabinom(k, size=n, m=a/(a+b), s=a+b) * p0.prob)
	}

	hapn1 = 2*n1
	hapn2 = 2*n2

	# allele frequency difference probabilities represented by each joint SFS category (these sum to 1 minus the fixed categories)
	diff = NULL
	for(i in 0:hapn1) {
		p1.prob = bbprob(i, hapn1, fst, p[,1], p[,2])
		f1 = i/hapn1
		for (j in 0:hapn2) {
			if ((i == 0 && j == 0) || (i == hapn1 && j == hapn2)) next # consider only SNPs
			diff = rbind(diff, c(abs(f1 - j/hapn2), p1.prob * bbprob(j, hapn2, fst, p[,1], p[,2])))
		}
	}

	# collapse joint SFS difference probabilities by absolute frequency difference
	diffprob = NULL	
	if (is.null(binsize)) {
		diffprob = data.frame(diff=sort(unique(diff[,1])), prob=NA)
		for (i in 1:nrow(diffprob)) diffprob[i,2] = sum(diff[which(diff[,1] == diffprob[i,1]),2])
	} else {
		diffprob = data.frame(lower=seq(from=0,to=1-binsize,by=binsize), upper=seq(from=binsize,to=1,by=binsize), prob=0)
		diff[,1][diff[,1] == 0] = 1e-16
		for (i in 1:nrow(diff)) {
			idx = which(diffprob$lower < diff[,1][i] & diffprob$upper >= diff[,1][i])
			diffprob$prob[idx] = diffprob$prob[idx] + diff[,2][i]
		}
	}

	diffprob$prob = diffprob$prob/sum(diffprob$prob) # ensure sum to 1 without fixed sites
	
	# calculate CDF
	diffprob$cdf[1] = diffprob$prob[1]
	for (i in 2:nrow(diffprob)) diffprob$cdf[i] = diffprob$cdf[i-1]+diffprob$prob[i]
	
	return(diffprob)
}

### MAIN

## parse user input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1 || args[1] == "--help" || args[1] == "-h") {
cat("\nbetaAFOutlier.R <pop1n,pop2n> <fst> <allele frequency file> <out prefix> [options]\n
pop1n,pop2n: Comma-separated diploid sample sizes for populations 1 and 2.
fst: Genome-wide average FST value between groups being compared or '9' to estimate FST from input.
allele frequency file: TSV file with columns (1) chromosome, (2) position, (3) pop1 allele frequency (4) pop2 allele frequency. Assumes header.
out prefix: Output file name prefix to which '.dist', '.ancf', '.sigtest', and '.pdf' are appended
\nOptional Input
--ancmethod <0|1>: Ancestral allele frequency method where 0 (default) uses empirical ancestral frequency prior, 1 uses expected frequency prior
--ancbin <FLOAT in range [0,1]>: Ancestral allele frequency prior bin width [default: 0 (no binning)]
--plotqq: Generate qq-plot of p-values
--seed <INT>: Set a specific seed\n\n")

options(show.error.messages=FALSE)
stop()}

dat = parseArgs(args)
binwidth = dat$ancbin # width of ancestral allele frequency prior bins
if (is.null(dat$seed)) dat$seed = round(runif(1,min=1,max=10^9))
set.seed(dat$seed) # THINK THIS ONLY APPLIES TO METHOD 1
cat(paste0("\nSeed set to ",dat$seed,"\n"))

## hard-coded parameters
nsamp = 10^6 # number of allele frequency pairs to draw for each ancestral allele frequency
diffwidth = max(1/(2*dat$dipn[1]), 1/(2*dat$dipn[2])) # width of allele frequency difference bin

## Calculate background FST
if (dat$fstbar == 9) {
	cat("\nCalculating Genome-wide FST\n")
	dat$fstbar = genomeFst(dat$obsfreq[,3:4],dat$dipn[1],dat$dipn[2])
	cat("FST:",dat$fstbar,"\n")
} else {
	cat("User set FST:",dat$fstbar,"\n")
}

## METHOD 2
# calculate folded nonfixed ancestral allele frequency probabilities - TRY UNFOLDED WRT TO REF ALLELE
#p0 = empiric_sfs_prob(as.matrix(dat$obsfreq[,3:4]), dat$dipn[1], dat$dipn[2], fold=TRUE, bins=FALSE)

# calculate expected probabilties for all possible allele frequency differences
#nulldist = betaBinomDist(dat$dipn[1], dat$dipn[2], dat$fstbar,p0)

### END METHOD 2

## calculate ancestral allele frequency probabilities
if (binwidth) cat("\nAncestral frequency bin width set to ",binwidth,"\n") else cat("\nNot binning ancestral frequencies\n")
ancf = NULL
if (dat$ancmethod) {
	cat("Using expected ancestral allele frequencies\n\n")
	ancf = expected_af_prob(dat$dipn[3], binwidth)
} else {
	cat("Estimating ancestral allele frequencies from descendant pop frequencies\n\n")
	if (binwidth) {
		ancf = empiric_af_prob(dat$obsfreq[,3:4], binwidth, dat$dipn[1], dat$dipn[2], weight=TRUE) 
	} else ancf = empiric_sfs_prob(as.matrix(dat$obsfreq[,3:4]), dat$dipn[1], dat$dipn[2], fold=TRUE, bins=TRUE)
}

## sample allele frequencies for each ancestral frequency bin
afsamp = sampleFreq(nsamp, dat$fstbar, dat$dipn[1], dat$dipn[2], as.matrix(ancf[,c(1,2)]))

## calculate weighted absolute allele frequency difference probabilities
diffprobs = afDiffProbs(afsamp, diffwidth, ancwidth=(ancf[1,2]-ancf[1,1])) # allele difference probabilities for each ancestral frequency
nulldist = weightedDiffProbs(diffprobs, ancf)

## calculate p and q values for observed differences
cat("\nComparing observed to expected allele frequency differences\n")
difftest = testObserved(dat$obsfreq, nulldist)

## write results
cat("\nWriting results and plotting\n")
write.table(nulldist,file=dat$outnames[1],col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
write.table(ancf,file=dat$outnames[2],col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
write.table(difftest,file=dat$outnames[3],col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

## make plots

# plot fit
pdf(file=dat$outnames[4])
plotFit(difftest, nulldist)
invisible(dev.off())

# plot p-value qq
if (!is.null(dat$qq)) {
	cat("Making p-value qq plot, which may take a while...\n")
	png(file=dat$outnames[5])
	pvalqq(difftest$pval)
	invisible(dev.off())
}
