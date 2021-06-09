library(SummarizedExperiment)
library(PharmacoGx)

load('gCSI.sigs.RData')
load('GDSCv2.sigs.RData')
load('CTRPv2.sigs.RData')

stopifnot(length(unique(c(nrow(gCSI.sig.aac), nrow(GDSCv2.sig.aac),nrow(CTRPv2.sig.aac))))==1)

Kvec <- 10^seq(1,4,0.1)


ngene <- nrow(gCSI.sig.aac)

drugs <- c("Bortezomib", "Entinostat", "Sirolimus","Docetaxel","Gemcitabine", "Crizotinib",  "Lapatinib","Vorinostat","Erlotinib","Paclitaxel","Pictilisib")


jaccard_nway <- function(sig_list, drug, K, rank_on="estimate",abs=TRUE,decreasing=FALSE){

	gene.list <- lapply(sig_list, function(sig){
		if(abs){
			if(decreasing){
				rownames(sig)[rank(-abs(sig[,drug,rank_on]))<=K]
			} else {
				rownames(sig)[rank(abs(sig[,drug,rank_on]))<=K]
			}
		}else{
			if(decreasing){
			rownames(sig)[rank(-(sig[,drug,rank_on]))<=K]

			} else {
			rownames(sig)[rank((sig[,drug,rank_on]))<=K]

			}
		}
	})
	return(length(.intersectList(gene.list))/length(.unionList(gene.list)))
}


estimate.aucs <- list()
first.k <- list()

pdf("3wayJaccard.pdf", onefile=TRUE)
for(drug in drugs){


	aac.jac <- sapply(Kvec, function(k)
		jaccard_nway(list(gCSI.sig.aac, GDSCv2.sig.aac, CTRPv2.sig.aac)
			, drug, k, decreasing=TRUE))
	ic50_recomputed.jac <- sapply(Kvec, function(k)
		jaccard_nway(list(gCSI.sig.ic50_recomputed, GDSCv2.sig.ic50_recomputed, CTRPv2.sig.ic50_recomputed)
			, drug, k, decreasing=TRUE))
	ic50_log.jac <- sapply(Kvec, function(k)
		jaccard_nway(list(gCSI.sig.ic50_log, GDSCv2.sig.ic50_log, CTRPv2.sig.ic50_log)
			, drug, k, decreasing=TRUE))
	ic50_logtunc.jac <- sapply(Kvec, function(k)
		jaccard_nway(list(gCSI.sig.ic50_logtunc, GDSCv2.sig.ic50_logtunc, CTRPv2.sig.ic50_logtunc)
			, drug, k, decreasing=TRUE))
	ic50_trunc.jac <- sapply(Kvec, function(k)
		jaccard_nway(list(gCSI.sig.ic50_trunc, GDSCv2.sig.ic50_trunc, CTRPv2.sig.ic50_trunc)
			, drug, k, decreasing=TRUE))
	null_jac <- (Kvec^3/ngene^2)/(3*Kvec - 3*Kvec^2/ngene + Kvec^3/ngene^2) #one power of nrow factored out


	ylim <- c(min(c(min(aac.jac),
		min(ic50_recomputed.jac),
		min(ic50_log.jac),
		min(ic50_logtunc.jac),
		min(ic50_trunc.jac),
	    min(null_jac))),
		max(c(max(aac.jac),
		max(ic50_recomputed.jac),
		max(ic50_log.jac),
		max(ic50_logtunc.jac),
		# max(ic50_trunc.jac),
	    max(null_jac))))
	if(ylim[1]==0) ylim[1] <- min(null_jac)
	plot(Kvec, null_jac, ylim=ylim, col="gray",cex=1, pch=16,log="xy",ylab="Jaccard", xlab="K", main=paste0(drug, " - Ranking by estimate"))
	lines(Kvec, null_jac, ylim=ylim, col="gray")
	points(Kvec, aac.jac, ylim=ylim, col="blue",cex=1, pch=16)
	lines(Kvec, aac.jac, ylim=ylim, col="blue")
	points(Kvec, ic50_recomputed.jac, ylim=ylim, col="red",cex=1, pch=16)
	lines(Kvec, ic50_recomputed.jac, ylim=ylim, col="red")
	points(Kvec, ic50_log.jac, ylim=ylim, col="cyan",cex=1, pch=16)
	lines(Kvec, ic50_log.jac, ylim=ylim, col="cyan")	
	points(Kvec, ic50_logtunc.jac, ylim=ylim, col="black",cex=1, pch=16)
	lines(Kvec, ic50_logtunc.jac, ylim=ylim, col="black")
	points(Kvec, ic50_trunc.jac, ylim=ylim, col="magenta",cex=1, pch=16)
	lines(Kvec, ic50_trunc.jac, ylim=ylim, col="magenta")
	legend("topleft", legend=c("AAC","IC50", "Log IC50", "Trunc IC50", "Log Trunc IC50", "Null"),
						fill=c("blue", "red", "cyan", "magenta", "black", "gray"))


	aac.first.k <- min(Kvec[which( aac.jac!=0)])
	ic50_recomputed.first.k <- min(Kvec[which( ic50_recomputed.jac!=0)])
	ic50_log.first.k <- min(Kvec[which( ic50_log.jac!=0)])
	ic50_logtunc.first.k <- min(Kvec[which( ic50_logtunc.jac!=0)])
	ic50_trunc.first.k <- min(Kvec[which( ic50_trunc.jac!=0)])

	first.k[[drug]] <- c("AAC"=aac.first.k,
	"IC50" = ic50_recomputed.first.k,
	"Log IC50"=ic50_log.first.k,
	"Log Trunc IC50"=ic50_logtunc.first.k,
	'Trunc IC50' = ic50_trunc.first.k)


	# Calculating Area 

	aac.jac[aac.jac==0] <- null_jac[aac.jac==0]
	ic50_recomputed.jac[ic50_recomputed.jac==0] <- null_jac[ic50_recomputed.jac==0]
	ic50_log.jac[ic50_log.jac==0] <- null_jac[ic50_log.jac==0]
	ic50_logtunc.jac[ic50_logtunc.jac==0] <- null_jac[ic50_logtunc.jac==0]
	ic50_trunc.jac[ic50_trunc.jac==0] <- null_jac[ic50_trunc.jac==0]

	estimate.aucs[[drug]] <- c("AAC"=sum(log(aac.jac/null_jac)),
	"IC50" = sum(log(ic50_recomputed.jac / null_jac)),
	"Log IC50"=sum(log(ic50_log.jac / null_jac)),
	"Log Trunc IC50"=sum(log(ic50_logtunc.jac / null_jac)),
	'Trunc IC50' = sum(log(ic50_trunc.jac / null_jac)))

}


estimate.aucs <- do.call(rbind, estimate.aucs)


first.k <- do.call(rbind, first.k)


estimate.aucs <- estimate.aucs*diff(log10(range(Kvec)))/length(Kvec)

library(reshape2)
library(patchwork)

library(ggplot2)

estimate.aucs.m <- melt(estimate.aucs)

colnames(estimate.aucs.m) <- c("Drug", "Sensitivity Measure", "AUC")

ggplot(estimate.aucs.m, aes(Drug, AUC, fill=`Sensitivity Measure`)) + geom_bar(position="dodge", stat="identity") +
 theme(axis.text.x=element_text(angle=45,hjust=1)) + scale_fill_manual(values=c("blue", "red", "cyan", "black", "magenta"))


dev.off()


cbbPalette <- c("#CC79A7", "#E69F00", "#56B4E9","#000000", "#009E73")
pal=c("#E64B35FF", "#8491B4FF","#4DBBD5FF",  "#3C5488FF","#00A087FF"  )
p1 <- ggplot(estimate.aucs.m, aes(Drug, AUC, fill=`Sensitivity Measure`)) + geom_bar(position="dodge", stat="identity") + theme_classic() +
 theme(axis.text.x=element_text(angle=45,hjust=1)) + scale_fill_manual(values=pal) + labs(y = expression("Integral of J"[obs]*"/J"[expected]))

first.k.m <- melt(log10(first.k))

colnames(first.k.m) <- c("Drug", "Sensitivity Measure", "Log10 K")

p2 <- ggplot(first.k.m, aes(Drug, `Log10 K`, fill=`Sensitivity Measure`)) + geom_bar(position="dodge", stat="identity") + theme_classic() +
 theme(axis.text.x=element_text(angle=45,hjust=1)) + scale_fill_manual(values=pal)+ labs(y = expression("log"[10]*"(K"[min]*")"))

png("supp_figure_uni_stability.png", height=11, width=8.5, units = 'in', res=600)
plot_grid(p2, p1, labels = c('A)', 'B)'), label_size = 12, ncol=1)
dev.off()
## calculating statistics for paper
rowSums(apply(apply(-estimate.aucs, 1, rank), 2, function(x)return(x==min(x))))
rowSums(apply(apply(first.k, 1, rank), 2, function(x)return(x==min(x))))

rowSums(apply(apply(-estimate.aucs[,-1], 1, rank), 2, function(x)return(x==min(x))))
rowSums(apply(apply(first.k[,-1], 1, rank), 2, function(x)return(x==min(x))))







pval.aucs <- list()
pdf("3wayJaccardPVal.pdf", onefile=TRUE)
for(drug in drugs){


	aac.jac <- sapply(Kvec, function(k)
		jaccard_nway(list(gCSI.sig.aac, GDSCv2.sig.aac, CTRPv2.sig.aac)
			, drug, k, rank_on="pvalue"))
	ic50_recomputed.jac <- sapply(Kvec, function(k)
		jaccard_nway(list(gCSI.sig.ic50_recomputed, GDSCv2.sig.ic50_recomputed, CTRPv2.sig.ic50_recomputed)
			, drug, k, rank_on="pvalue"))
	ic50_log.jac <- sapply(Kvec, function(k)
		jaccard_nway(list(gCSI.sig.ic50_log, GDSCv2.sig.ic50_log, CTRPv2.sig.ic50_log)
			, drug, k, rank_on="pvalue"))
	ic50_logtunc.jac <- sapply(Kvec, function(k)
		jaccard_nway(list(gCSI.sig.ic50_logtunc, GDSCv2.sig.ic50_logtunc, CTRPv2.sig.ic50_logtunc)
			, drug, k, rank_on="pvalue"))
	ic50_trunc.jac <- sapply(Kvec, function(k)
		jaccard_nway(list(gCSI.sig.ic50_trunc, GDSCv2.sig.ic50_trunc, CTRPv2.sig.ic50_trunc)
			, drug, k, rank_on="pvalue"))
	null_jac <- (Kvec^3/ngene^2)/(3*Kvec - 3*Kvec^2/ngene + Kvec^3/ngene^2) #one power of nrow factored out


	ylim <- c(min(c(min(aac.jac),
		min(ic50_recomputed.jac),
		min(ic50_log.jac),
		min(ic50_logtunc.jac),
		min(ic50_trunc.jac),
	    min(null_jac))),
		max(c(max(aac.jac),
		max(ic50_recomputed.jac),
		max(ic50_log.jac),
		max(ic50_logtunc.jac),
		# max(ic50_trunc.jac),
	    max(null_jac))))
	if(ylim[1]==0) ylim[1] <- min(null_jac)
	plot(Kvec, null_jac, ylim=ylim, col="gray",cex=1, pch=16,log="xy",ylab="Jaccard", xlab="K", main=paste0(drug, " - Ranking by pvalue"))
	lines(Kvec, null_jac, ylim=ylim, col="gray")
	points(Kvec, aac.jac, ylim=ylim, col="blue",cex=1, pch=16)
	lines(Kvec, aac.jac, ylim=ylim, col="blue")
	points(Kvec, ic50_recomputed.jac, ylim=ylim, col="red",cex=1, pch=16)
	lines(Kvec, ic50_recomputed.jac, ylim=ylim, col="red")
	points(Kvec, ic50_log.jac, ylim=ylim, col="cyan",cex=1, pch=16)
	lines(Kvec, ic50_log.jac, ylim=ylim, col="cyan")	
	points(Kvec, ic50_logtunc.jac, ylim=ylim, col="black",cex=1, pch=16)
	lines(Kvec, ic50_logtunc.jac, ylim=ylim, col="black")
	points(Kvec, ic50_trunc.jac, ylim=ylim, col="magenta",cex=1, pch=16)
	lines(Kvec, ic50_trunc.jac, ylim=ylim, col="magenta")
	legend("topleft", legend=c("AAC","IC50", "Log IC50", "Trunc IC50", "Log Trunc IC50", "Null"),
						fill=c("blue", "red", "cyan", "magenta", "black", "gray"))

	# Calculating Area 

	aac.jac[aac.jac==0] <- null_jac[aac.jac==0]
	ic50_recomputed.jac[ic50_recomputed.jac==0] <- null_jac[ic50_recomputed.jac==0]
	ic50_log.jac[ic50_log.jac==0] <- null_jac[ic50_log.jac==0]
	ic50_logtunc.jac[ic50_logtunc.jac==0] <- null_jac[ic50_logtunc.jac==0]
	ic50_trunc.jac[ic50_trunc.jac==0] <- null_jac[ic50_trunc.jac==0]

	pval.aucs[[drug]] <- c("AAC"=sum(log(aac.jac/null_jac)),
	"IC50" = sum(log(ic50_recomputed.jac / null_jac)),
	"Log IC50"=sum(log(ic50_log.jac / null_jac)),
	"Log Trunc IC50"=sum(log(ic50_logtunc.jac / null_jac)),
	'Trunc IC50' = sum(log(ic50_trunc.jac / null_jac)))


}
dev.off()

pval.aucs <- do.call(rbind, pval.aucs)


pval.aucs <- pval.aucs*diff(log10(range(Kvec)))/length(Kvec)


