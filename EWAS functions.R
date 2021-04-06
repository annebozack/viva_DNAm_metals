
# libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(minfi)
library(stringr)
library(limma)
library(ENmix)

# functions
# manhattan polt
manhattan = function(probe, region, FDR = FALSE, annotate = NULL, title = NULL){
	# chromosome as numeric
	probe$chr = as.numeric(gsub("chr", "", probe$chr))

	# dataframes with common columns
	probe = data.frame(name = probe$cpg, p = probe$P.Value, fdr = probe$adj.P.Val, bonf = probe$adj.P.Val.bonf, chr = probe$chr, pos = probe$pos, type = 'position', dmp_sig = NA, color = NA)
	region_df = data.frame(name = as.character(seq(1:nrow(region))), p = NA, fdr = NA, bonf = NA, chr = region$chr, pos = region$start, type = 'region', dmp_sig = NA, color = NA, size = NA)

	# indicating probes within regions
	for (i in 1:length(region$chr)){
		probe$dmp_sig[(probe$chr == region$chr[i]) & (probe$pos >= region$start[i]) & (probe$pos <= region$end[i])] = 1
	}
	# variable for point color
	probe$color[probe$dmp_sig == 1] = 50
	probe$color[probe$fdr < 0.05] = 100
	probe$color[is.na(probe$color)] = probe$chr[is.na(probe$color)]
	# variable for point size
	probe$size[probe$color == 50 | probe$color == 100] = 1
	probe$size[is.na(probe$size)] = 0.5

	# combine dataframes
	df_comb = rbind(probe, region_df)

	# dataset for plotting
		don = df_comb %>% 
	  		# Compute chromosome size
				group_by(chr) %>% summarise(chr_len=as.numeric(max(pos))) %>% 
	  		# Calculate cumulative position of each chromosome
	 			mutate(tot=cumsum(chr_len)-chr_len) %>% dplyr::select(-chr_len) %>%
	  		# Add this info to the initial dataset
	  			left_join(df_comb, ., by=c("chr"="chr")) %>%
	  		# Add a cumulative position of each site
				arrange(chr, pos) %>% mutate(poscum=pos+tot) # %>%
	 		# Prepare X axis
				axisdf = don %>% group_by(chr) %>% summarize(center=(max(poscum) + min(poscum))/2)
		don = merge(don, df_comb, on='name', all.x=T)
		don = don[order(don$size),]
		don_position = don[don$type == 'position',]
		don_region = don[don$type == 'region',]
	
		colors = c("#969696", "#737373", "#969696", "#737373", "#969696", "#737373", "#969696", "#737373", "#969696", "#737373", "#969696", "#737373", "#969696", "#737373", "#969696", "#737373", "#969696", "#737373", "#969696", "#737373", "#969696", "#737373", '#2166ac', 'black')
	
	manhattan = ggplot(don_position, aes(x=poscum, y=-log10(p))) +
	geom_point(aes(color=as.factor(color)), size= don_position$size, alpha = don_position$size) + scale_color_manual(values = colors) +
    # p-value cutoffs
	geom_hline(yintercept=-log10(0.05/nrow(don_position)), colour = '#AB3428', size=.2, alpha = 0.5) +
	geom_vline(xintercept= don_region$poscum, colour = '#4393c3', size=.2) +
	# custom axes:
	scale_x_continuous(expand = c(0.005, 0.005), limits = c(min(don_position$poscum), max(don_position$poscum)), label = axisdf$chr, breaks= axisdf$center) +
	scale_y_continuous(expand = c(0, 0), limits = c(0, (max(-log10(don_position$p)) + 0.5)), breaks = seq(from = 0, to = (max(-log10(don_position$p)) + 0.5), by = 1)) +
	# Custom theme:
    theme_minimal() + theme( 
	legend.position="none", panel.border = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), text = element_text(size = 7.5)) + 
    labs(y=expression(-log[10](italic(p))), x='Chromosome') 
    if (!is.null(title)){
    	manhattan = manhattan + labs(title = title)
    }
    if (FDR == TRUE){
    	manhattan = manhattan + geom_hline(yintercept=-log10(max(don_position$p[don_position$fdr < 0.05])), colour='#AB3428', size=.2, alpha = 0.5, linetype = "dashed")
    } 
    if (!is.null(annotate)){
		manhattan = manhattan + annotate("text", x = max(don_position$poscum)*0.05, y = max(-log10(don_position$p)), label = annotate, size = 4)
	}
	return(manhattan)
}

# lambda
lambda = function(p) median(qchisq(p, df=1, lower.tail=FALSE), na.rm=TRUE) / qchisq(0.5, df=1)

# QQ plot
gg_qqplot = function(pvector){
	l = round(lambda(pvector), 3)
	o = -log10(sort(pvector, decreasing = FALSE))
	e = -log10(ppoints(length(pvector)))
	df = data.frame(o = o, e = e)
	ggplot(df, aes(e, o)) + geom_point(alpha = 0.5, size = 1) + geom_abline(intercept = 0, slope = 1, color = '#AB3428') + labs(y = expression(Observed ~ ~-log[10](italic(p))), x = expression(Expected ~ ~-log[10](italic(p)))) + theme_classic() + annotate("text", x = 1, y = 5, label = paste0('lambda = ', l))
}

# volcano
volcano = function(probe, FDR = FALSE){
		volcano = ggplot(probe, aes(x=logFC, y  = -log10(P.Value))) + 
		geom_point(size = 0.8, alpha=0.4) + 
		geom_hline(aes(yintercept = -log10(0.05/nrow(probe))), color = "#AB3428", size = 0.5) +  
		scale_linetype_manual(name = '', values = c(1,2), guide = guide_legend(override.aes = list(color = c("#AB3428", "#AB3428")))) + theme_minimal() + 
		labs(y=expression(-log[10]*"(P-value)"), x='Effect estimate') + theme(panel.grid.minor.y = element_blank()) + theme(text = element_text(size=8)) + 
		theme(legend.position="none") + scale_y_continuous(expand = c(0, 0)) + theme(panel.grid.minor.x = element_blank(), panel.grid.major.y = element_line(size = 0.2, color = 'gray65'), panel.grid.major.x = element_line(size = 0.2, color = 'gray65'))
		if (FDR == TRUE){
			volcano = volcano + geom_hline(aes(yintercept = -log10(max(P.Value[adj.P.Val < 0.05]))), color = "#AB3428", size = 0.5, linetype = "dashed")
		}
		return(volcano)
}


# DMP analysis
run_DMP <- function(mvals, design){
  # fit model
  l_fit <- limma::lmFit(object = mvals, design = design)
  
  # extract standard errors
  std_err <- l_fit$stdev.unscaled[,2]*l_fit$sigma
  std_err_df <- data.frame(std_err)
  std_err_df$cpg <- rownames(std_err_df)
  
  e_fit <- limma::eBayes(l_fit, robust = TRUE)
  
  # extract results and add Bonferroni correction
  p_top <- limma::topTable(e_fit, adjust = "BH", coef = 2, num = Inf, confint = TRUE)
  p_top <- p_top[order(p_top$P.Value), , drop = FALSE]
  p_top$adj.P.Val.bonf <- topTable(e_fit, adjust="bonferroni", coef=2, number = Inf)$adj.P.Val
  
  # merge results and standard errors
  p_top$cpg <- rownames(p_top)
  p_top <- merge(p_top, std_err_df, by = 'cpg')
  rownames(p_top) <- p_top$cpg
  
  return(p_top)
}

# Combp
acf.table<-function(x,loc,dist.cutoff){
  flag=TRUE; lag=1; result=NULL
  while(flag){
    x1=head(x,-lag); x2=tail(x,-lag); dist=diff(loc,lag=lag)
    index=(dist<dist.cutoff)  
    if(all(!index)){flag=FALSE}else{
      result=rbind(result,data.frame(x1=x1[index],x2=x2[index],dist=dist[index]))
    lag=lag+1
    }
  }
  return(result)  
}

get.acf<-function(data,dist.cutoff,bin.size){
  temp<-NULL
  for (chr in unique(as.vector(data$chr))){
    y<-data[as.vector(data$chr)==chr,]; y<-y[order(y$end),]
    temp<-rbind(temp,acf.table(y$p,y$end,dist.cutoff))
  }
  bin.label<-findInterval(temp$dist,seq(bin.size,dist.cutoff,bin.size))
  temp.stouffer<-by(temp,bin.label,FUN=function(x){cor.test(qnorm(x$x1),
               qnorm(x$x2),alternative="greater")},simplify=FALSE)

  cor.stouffer<-sapply(temp.stouffer,function(x){x$estimate})
  p.stouffer<-sapply(temp.stouffer,function(x){x$p.value})

  if (any(p.stouffer>0.05)){
    index=min(which(p.stouffer>0.05))
    cor.stouffer[index:length(cor.stouffer)]=0
  }
  return(cor.stouffer)
}

# comb_p-like method
combp2<-function(data,dist.cutoff=1000,bin.size=310,seed=0.01,
               region_plot=TRUE,mht_plot=TRUE,nCores=10,verbose=TRUE){
  if(nCores>detectCores()){nCores=detectCores()}
  data=as.data.frame(data)
  data$start=as.numeric(as.vector(data$start))
  data$end=as.numeric(as.vector(data$end))
  data=data[!is.na(data$start) & !is.na(data$end),]
  data$p=as.numeric(as.vector(data$p))

  acf<-get.acf(data,dist.cutoff,bin.size)
  if(verbose){
    cat("P value correlations:\n")
    bin=seq(bin.size,dist.cutoff,bin.size)
    if(!(dist.cutoff%%bin.size==0)){bin=c(bin,dist.cutoff)}
    print(data.frame(bin=bin,acf=acf))
  }

  result<-mclapply(unique(as.vector(data$chr)), function(chr){
    y=data[as.vector(data$chr)==chr,]; y=y[order(y$end),]
    pos=y$end; p=qnorm(y$p)

    temp=sapply(pos,function(i){
      index.i=(abs(pos-i)<bin.size);
      if (sum(index.i)>1){  
        int<-findInterval(c(dist(pos[index.i])),c(bin.size,2*bin.size))
        sd<-sqrt(sum(acf[int+1])*2+sum(index.i))
        return(pnorm(sum(p[index.i]),mean=0,sd=sd))
      }else{return(y$p[index.i])}
    })

    return(data.frame(chr,start=pos,end=pos,s.p=temp))
  },mc.cores=nCores)

  result <- do.call("rbind", result)
  names(result)=c("chr","start","end","s.p")

  result=result[p.adjust(result$s.p,method="fdr")<seed,]

  result.fdr=NULL
  if (nrow(result)>0){
    for (chr in unique(result$"chr")){
      y=data[as.vector(data$chr)==chr,]; y=y[order(y$end),]
      pos=y$end; p=qnorm(y$p)

      result.chr=result[result$"chr"==chr,]
      a=IRanges::IRanges(start=result.chr$start,end=result.chr$end)
      b=IRanges::reduce(a,min.gapwidth=dist.cutoff)

      start=IRanges::start(b); end=IRanges::end(b)
      region.max<-max(IRanges::width(b))
      temp=sapply(1:length(b),function(i){
        index.i=(pos>=start[i] & pos<=end[i]);
        if (sum(index.i)>1){  
          int<-findInterval(c(dist(pos[index.i])),
              seq(bin.size,region.max+bin.size,bin.size))
          sd<-sqrt(sum(ifelse(int<length(acf),
              acf[int+1],0))*2+sum(index.i))
          return(pnorm(sum(p[index.i]),mean=0,sd=sd))
        }else{return(y$p[index.i])}
      })
      result.fdr=rbind(result.fdr,data.frame(chr,start,end,p=temp))
    }
    
    result.fdr$length = (result.fdr$end - result.fdr$start) + 1
    result.fdr = result.fdr[result.fdr$length > 1,]

	if (nrow(result.fdr)>0){
		
		##### BH FDR correction and Sidak correction
    	result.fdr$fdr=p.adjust(result.fdr$p,method="fdr")
    	result.fdr$sidak=(1-(1-result.fdr$p)^(nrow(data)/(result.fdr$end-result.fdr$start+1)))
    	result.fdr<-result.fdr[order(result.fdr$p),]

    	##### use 0-coordinate
    	result.fdr$start=(result.fdr$start-1)
	} else {
		result.fdr = NULL
	}
  }

  if(is.null(result.fdr)){cat("Number of identified DMR:  0\n")}else{
    ndmr=nrow(result.fdr)
  result.fdr$start=as.numeric(as.vector(result.fdr$start))
  result.fdr$end=as.numeric(as.vector(result.fdr$end))
  result.fdr$chr=factor(result.fdr$chr)

    cat("Number of DMRs identified:  ",ndmr, "\n")
    if(region_plot){
      cat("Drawing regional plot: region_plot.pdf ...\n")
      sig=result.fdr
      regplot(ref=data,sig)
    }
  if(mht_plot){
    cat("Drawing manhattan plot: mht.jpg ...\n")
    set2=NULL
    for(i in 1:ndmr){
        set2=c(set2,as.vector(data$probe[as.vector(data$chr)==as.vector(result.fdr$chr[i])
           & data$start>=result.fdr$start[i] & data$start<=result.fdr$end[i]]))
    }
  mhtplot(probe=data$probe,chr=as.vector(data$chr),pos=data$start,p=data$p,color="gray",markprobe=set2)
  }
  #number of probes within eath DMR

  result.fdr$nprobe=NA
  for(i in 1:nrow(result.fdr)){
result.fdr$nprobe[i]=nrow(data[as.vector(data$chr)==as.vector(result.fdr$chr[i])
& data$start>=result.fdr$start[i] & data$end<=result.fdr$end[i],])
}

  write.table(result.fdr,"resu_combp.csv",row.names=FALSE,sep=",")
  }
}

run_EWAS = function(DNAm, pheno, var, covar, anno, path, save_adj = TRUE){
	cat('###########\n')
	cat('#', var, '#\n')
	cat('###########\n', '\n')
	
	# set wd
		if(!dir.exists(path)){
			dir.create(path)
		}
		setwd(path)
		
	# unadjusted EWAS
		modUnadj = model.matrix(formula(paste('~', var)), pheno)
		DNAm_comp = DNAm[,colnames(DNAm) %in% rownames(modUnadj)]
		DNAm_comp = DNAm_comp[,match(rownames(modUnadj), colnames(DNAm_comp))]
		cat('Unadjusted, N = ', ncol(DNAm_comp), '\n')
		DMP_unadj <- run_DMP(mvals = DNAm_comp, design = modUnadj)
		cat('Unadjusted, p<0.05: ', sum(DMP_unadj$P.Value < 0.05), '\n')
		cat('Unadjusted, FDR<0.05: ', sum(DMP_unadj$adj.P.Val < 0.05), '\n')
		cat('Unadjusted, pBonf<0.05: ', sum(DMP_unadj$adj.P.Val.bonf < 0.05), '\n', '\n')
		
	# adjusted EWAS
		varMod = paste(c(var, covar), collapse = ' + ')
		modAdj = model.matrix(formula(paste('~', varMod)), pheno)
		DNAm_comp = DNAm[,colnames(DNAm) %in% rownames(modAdj)]
		DNAm_comp = DNAm_comp[,match(rownames(modAdj), colnames(DNAm_comp))]
		cat('Adjusted, N = ', ncol(DNAm_comp), '\n')
		DMP_adj <- run_DMP(mvals = DNAm_comp, design = modAdj)
		cat('Adjusted, p<0.05: ', sum(DMP_adj$P.Value < 0.05), '\n')
		cat('Adjusted, FDR<0.05: ', sum(DMP_adj$adj.P.Val < 0.05), '\n')
		cat('Adjusted, pBonf<0.05: ', sum(DMP_adj$adj.P.Val.bonf < 0.05), '\n')
		cat('Adjusted, lambda: ', lambda(DMP_adj$P.Value), '\n', '\n')
		
		DMP_adj = cbind(DMP_adj, anno[,'chr'][match(rownames(DMP_adj), rownames(anno))])
		DMP_adj = cbind(DMP_adj, anno[,'pos'][match(rownames(DMP_adj), rownames(anno))])
		DMP_adj = cbind(DMP_adj, anno[,'UCSC_RefGene_Name'][match(rownames(DMP_adj), rownames(anno))])
		colnames(DMP_adj)[c(12,13,14)] = c('chr', 'pos', 'gene')

		# significant CpGs
			if (sum(DMP_adj$adj.P.Val < 0.05) > 0){
				DMP_adj_sig = DMP_adj[DMP_adj$adj.P.Val < 0.05,]
				DMP_adj_sig$logFC_CI = paste0(round(DMP_adj_sig$logFC, 2), ' (', round(DMP_adj_sig$CI.L, 2), ' ,', round(DMP_adj_sig$CI.R, 2), ')')
				DMP_adj_sig$AveExpr = round(DMP_adj_sig$AveExpr,2)
				DMP_adj_sig$adj.P.Val.bonf = round(DMP_adj_sig$adj.P.Val.bonf,4)
				DMP_adj_sig = DMP_adj_sig[,c('logFC_CI', 'AveExpr', 'P.Value', 'adj.P.Val', 'adj.P.Val.bonf', 'chr', 'pos', 'gene')]
			} else {
				DMP_adj_sig = NA
			}
			FDR = F
			if (sum(DMP_adj$adj.P.Val < 0.05) > 1){
				FDR = T
			}	
		
		# combp
			DMP_adj$end = DMP_adj[,13]
			DMP_adj_p = DMP_adj[,c(12,13,15,7,1)]
			colnames(DMP_adj_p) = c( "chr", "start", "end","p", "probe")
			DMP_adj_p$chr = as.numeric(gsub("chr", "", DMP_adj_p$chr))

			if(!dir.exists(paste0(path, '/', 'combp'))){
				dir.create(paste0(path, '/', 'combp'))
			}
			setwd(paste0(path, '/', 'combp'))
			combp2(DMP_adj_p, region_plot = F, mht_plot = F, seed = 0.001, verbose = T)
			
			if ("resu_combp.csv" %in% list.files()){
				region = read.csv(paste0(path, '/', 'combp', '/', 'resu_combp.csv'))
			}
			setwd(path)
		
		# plotting
			qq = gg_qqplot(DMP_adj$P.Value)
			ggsave(filename = paste0(path, '/', var, '_QQ_DMP_adj.png'), plot=qq, device = 'png', width = 4, height = 4, units = 'in', dpi = 300)
			
			vol = volcano(DMP_adj, FDR)
			ggsave(filename = paste0(path, '/', var, '_volcano_DMP_adj.png'), plot=vol, device = 'png', width = 4, height = 4, units = 'in', dpi = 300)


			if (sum('region' %in% objects()) == 0){
				region = data.frame(chr = NA, start = NA, end = NA, p = NA, length = NA, fdr = NA, sidak = NA, nprobe = NA)
			} 
			probe = DMP_adj
			man = manhattan(probe, region, FDR)
			ggsave(filename = paste0(path, '/', var, '_manhattan_DMP_adj.png'), plot=man, device = 'png', width = 5.5, height = 3, units = 'in', dpi = 300)
			
		# saving results
			if (save_adj == TRUE){
				write.csv(DMP_adj, paste0(path, '/', var, '_DMP_adj.csv'))
			}
			
	rm(modUnadj, modAdj, DNAm_comp, DMP_unadj, DMP_adj_p);gc()	
	
	return(list(DMP_adj_sig, region, DMP_adj))
}



# Global DNAm

TestDensities2 = function (Z, y, X = NULL, outcomeType = "C", histBreaks = 200, 
    lambdas = c(0, exp(-10:10)), kernel = "linear", hideProgress = FALSE, 
    adjustmentDichot = FALSE, knotLocs = c(seq(0, 0.3, by = 0.02), 
        seq(0.32, 0.62, by = 0.1), seq(0.72, 1, by = 0.02))) 
{
    require(fda)
    require(SKAT)
    if (!("matrix" %in% class(Z))) 
        stop("Z must be a matrix!")
    if (ncol(Z) != length(y)) 
        stop("Dimensions of Z and y do not match!")
    if (!is.null(X)) {
        if (!("matrix" %in% class(X))) 
            stop("X must be a matrix!")
        if (nrow(X) != length(y)) 
            stop("Dimensions of X and y do not match!")
    }
    getHist = function(x, breaks) {
        ints = hist(x, breaks = breaks, plot = F)$density
        return(ints)
    }
    getGCV = function(lambda, H1, basis, argvals) {
        fdP = fdPar(basis, lambda = lambda)
        F.h1 = smooth.basis(argvals, H1, fdP)
        return(F.h1$gcv)
    }
    if (!hideProgress) 
        print("1. Estimating Histograms")
    breaks.h = seq(0, 1, length = histBreaks)
    H = apply(Z, 2, getHist, breaks.h)
    if (!hideProgress) 
        print("2. Computing Bases")
    breaks.b = knotLocs
    basis = create.bspline.basis(norder = 4, breaks = breaks.b)
    vals = 0.5 * (breaks.h[-1] + breaks.h[-histBreaks])
    if (!hideProgress) 
        print("3. Finding Optimal Lambda")
    if (length(lambdas) > 1) {
        gcvs = sapply(lambdas, getGCV, H, basis, vals)
        avgGCV = apply(gcvs, 2, mean)
        lambda = lambdas[which.min(avgGCV)]
    }
    else {
        lambda = lambdas
    }
    if (!hideProgress) 
        print("4. Obtaining B-spline coefficients")
    fdP = fdPar(basis, lambda = lambda)
    F.h = smooth.basis(vals, H, fdP)$fd
    Z1 = coef(F.h)
    if (!hideProgress) 
        print("5. Testing")
    if (!is.null(X)) {
        nullMod = SKAT_Null_Model(y ~ X, out_type = outcomeType, 
            Adjustment = adjustmentDichot)
    }
    else {
        nullMod = SKAT_Null_Model(y ~ 1, out_type = outcomeType, 
            Adjustment = adjustmentDichot)
    }
    res = SKAT(t(Z1), nullMod, is_check_genotype = F, kernel = kernel, 
        method = "liu")
    return(res$p.value)
}


TestCDFs2 = function (Z, y, X = NULL, outcomeType = "C", histBreaks = 1000, 
    numBases = min(histBreaks/4, 35), lambdas = c(0, exp(-10:10)), 
    kernel = "linear", hideProgress = FALSE, adjustmentDichot = FALSE, 
    knotLocs = seq(0, 1, length = numBases)) 
{
    require(fda)
    require(SKAT)
    if (!("matrix" %in% class(Z))) 
        stop("Z must be a matrix!")
    if (ncol(Z) != length(y)) 
        stop("Dimensions of Z and y do not match!")
    if (!("matrix" %in% class(X))) {
        if (class(X) != "matrix") 
            stop("X must be a matrix!")
        if (nrow(X) != length(y)) 
            stop("Dimensions of X and y do not match!")
    }
    getDist = function(x, breaks) {
        ints = (ecdf(x))(breaks)
        return(ints)
    }
    getGCV = function(lambda, H1, basis, argvals) {
        fdP = fdPar(basis, lambda = lambda)
        F.h1 = smooth.basis(argvals, H1, fdP)
        return(F.h1$gcv)
    }
    if (!hideProgress) 
        print("1. Calculating Empirical CDFs")
    breaks.h = seq(0, 1, length = histBreaks)
    vals = 0.5 * (breaks.h[-1] + breaks.h[-histBreaks])
    H = apply(Z, 2, getDist, vals)
    if (!hideProgress) 
        print("2. Computing Bases")
    breaks.b = knotLocs
    basis = create.bspline.basis(norder = 4, breaks = breaks.b)
    if (!hideProgress) 
        print("3. Finding Optimal Lambda")
    if (length(lambdas) > 1) {
        gcvs = sapply(lambdas, getGCV, H, basis, vals)
        avgGCV = apply(gcvs, 2, mean)
        lambda = lambdas[which.min(avgGCV)]
    }
    else {
        lambda = lambdas
    }
    if (!hideProgress) 
        print("4. Obtaining B-spline coefficients")
    fdP = fdPar(basis, lambda = lambda)
    F.h = smooth.basis(vals, H, fdP)$fd
    Z1 = coef(F.h)
    if (!hideProgress) 
        print("5. Testing")
    if (!is.null(X)) {
        nullMod = SKAT_Null_Model(y ~ X, out_type = outcomeType, 
            Adjustment = adjustmentDichot)
    }
    else {
        nullMod = SKAT_Null_Model(y ~ 1, out_type = outcomeType, 
            Adjustment = adjustmentDichot)
    }
    res = SKAT(t(Z1), nullMod, is_check_genotype = F, kernel = kernel, 
        method = "liu")
    return(res$p.value)
}


