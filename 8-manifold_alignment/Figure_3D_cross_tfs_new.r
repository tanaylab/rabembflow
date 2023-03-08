
library("metacell")
shades = colorRampPalette(c("darkblue","white","darkred"))(1000)
scdb_init("scdb_embexe/", force_reinit=T)
mm_mc = scdb_mc("mm_embexe")
oc_mc = scdb_mc("oc_embexe")

g_orto = intersect(rownames(mm_mc@e_gc), rownames(oc_mc@e_gc))

rbom = read.csv("../oc_mm_gastru/emb_rbom_match.csv",h=T, stringsAsFactors=F)

mm_rbom_legc = log2(mm_mc@e_gc[g_orto, rbom$mm_embexe_mcid]+1e-5)
oc_rbom_legc = log2(oc_mc@e_gc[g_orto, rbom$oc_embexe_mcid]+1e-5)

g_nz3 = names(which(rowSums(mm_rbom_legc > -15)>2 & rowSums(oc_rbom_legc > -15)>2))

gcor = tgstat::tgs_cor(t(mm_rbom_legc), t(oc_rbom_legc))

oc_gcor_max = apply(gcor[g_nz3,g_nz3], 2, max)

oc_max_a = apply(gcor[g_nz3,g_nz3], 2, which.max)
mm_max_a = apply(gcor[g_nz3,g_nz3], 1, which.max)

all_tfs = toupper(read.table("all_tf_nms.txt", stringsAsFactors=F)$V1)
recip = data.frame(g_oc =names(oc_max_a), 
		g_mm = names(mm_max_a)[oc_max_a], 
		g_oc_back =names(oc_max_a)[mm_max_a[oc_max_a]], 
		gc = round(oc_gcor_max,3))
recip$is_tf = recip$g_mm %in% all_tfs

self_recip = recip[recip[,1]==recip[,2] & recip[,1]==recip[,3],]

library(tgstat)

core_tfs_all = self_recip[self_recip$is_tf,1]
core_tfs = self_recip[self_recip$gc > 0.7 & self_recip$is_tf,1]
core_cross = tgs_cor(t(mm_rbom_legc[core_tfs,]), t(oc_rbom_legc[core_tfs,]))

hc_core = hclust(tgs_dist(core_cross))
n_g = length(core_tfs)
png("core_tfs.png", w=1000, h=1000)
par(mar=c(8,2,2,8))
image(core_cross[hc_core$order, hc_core$order], 
			col=shades,breaks=seq(-1,1,l=1001),
			xaxt='n', yaxt='n')
mtext(core_tfs[hc_core$order], at=seq(0,1,l=n_g),las=2,side=4, cex=1)
mtext(core_tfs[hc_core$order], at=seq(0,1,l=n_g),las=2,side=1, cex=1,las=2)

dev.off()


length(core_tfs[hc_core$order])

core_tfs_all


