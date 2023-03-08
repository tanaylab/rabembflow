do_init=T
do_endo=T
do_troph=T

do_pdf = T
r_pixel = 0.01

fshades = colorRampPalette(c("darkblue","white","darkred"))(1000)
rshades = colorRampPalette(c("white","red", "darkred"))(1000)

# +
mm_color_ord = c(
'#FACB12', #PGC
'#f7f79e', #Surface ectoderm
'#C3C388', #Neural crest
'#CDE088', #Neural tube/Floor plate
'#647A4F', #Forebrain/Midbrain/Hindbrain
'#354E23', #Caudal neural plate
'#649146', #Rostral neural plate
'#879E86', #Neural plate boundary
'#90BF75', #Definitive ectoderm
'#635547', #Epiblast
'#DABE99', #Primitive streak
'#9e6762', #Caudal epiblast
'#8e7e59', #Tail bud - neural
'#703C37', #Tail bud - mesoderm
'#C594BF', #Early nascent mesoderm
'#DFCDE4', #Late nascent mesoderm
'#1a3f52', #Caudal mesoderm
'#45d1c5', #Lateral & intermediate mesoderm
'#408DA1', #Paraxial mesoderm
'#A8DBF7', #Presomitic mesoderm
'#1AA2EB', #Somitic mesoderm
'#8DB5CE', #Rostral mesoderm
'#53f1fc', #Cardiopharyngeal mesoderm
'#B51D8D', #Cardiomyocytes
'#8870ad', #ExE mesoderm
'#cc7818', #Amnion/Chorion progenitor
'#824A09', #Amnion/Chorion
'#532C8A', #Allantois
'#ff891c', #Endothelial cells
'#FBBE92', #Haematoendothelial progenitors
'#C9A997', #Blood projenitors
'#C72228', #Erythroid 1
'#FF5608', #Erythroid 3
'#D96B2B', #Erythroid 2
'#c19f70', #Anterior Primitive Streak
'#0F4A9C', #Node/Notochord
'#F397C0', #Definitive endoderm
'#EF5A9D', #Gut
'#F25CD4', #Hindgut
'#BA9BA1', #Visceral endoderm - anterior
'#F6BFCB', #Visceral endoderm
'#7F6874', #ExE visceral endoderm
'#1A1A1A',  #Parietal endoderm
'#33CC33', # 4.5
'#ffad1e', #Chorion
'#ffd854', #Chorion progenitors
'#b2b2b2', #TSC2
'#19ce5b', #SpT-Gly
'#d6616b', #EPC progenitors
'#67000d', #TGC progenitors
'#2e7ebc', #p-TGC
'#ff78fa' #SpA-TGC
)

oc_color_ord = c( 
'#FACB12', #PGC
'#E2F700', #Median TFAP
'#f7f79e', #Surface ectoderm
'#CDE087', #Floor plate
'#C3C388', #Neural crest
'#879E86', #Neural plate boundary
# '#CDE089', #Neural tube Posterior
'#CDE088', #Neural tube Anterior
'#647A4F', #Forebrain/Midbrain/Hindbrain
'#354E23', #Caudal neurectoderm
'#649146', #rostral neurectoderm
'#90BF75', #Ectoderm - Definitive
'#635547', #epiblast
'#DABE99', #Primitive streak
'#9e6762', #caudal epiblast
'#C594BF', #nascent mesoderm
'#B6A8ED', #Nascent mesoderm - extraembryonic
'#1a3f52', #caudal mesoderm
'#45d1c5', #Lateral & intermediate mesoderm
'#408DA1', #Paraxial Mesoderm
'#A8DBF7', #Presomitic mesoderm
'#1AA2EB', #Somitic mesoderm
'#8DB5CE', #rostral mesoderm
'#53f1fc', #Cardiopharyngeal mesoderm
'#B51D8D', #Caridac
'#8870ad', #ExE mesoderm
'#cc7818', #Amnion
'#A69369', #chorionic mesothelium
'#532C8A', #Allantois
'#966F1B', #Endothelial2
'#ff891c', #Endothelial
'#FBBE92', #Hematoendothelial
'#c9a997', #Blood progenitors
'#C72228', #Erythroid1
'#FF5608', #Erythroid2
'#c19f70', #Anterior primitive streak
'#6666ED', #Node
'#0F4A9C', #notochord
'#F397C0', #definitive endoderm
'#EF5A9E', #Gut
'#EF5A9D', #Foregut
'#F25CD4', #Hindgut
'#CF6191', #Hypoblast-Anterior (Hhex/Lhx1/Cer1)
'#F6BFCB', #hypoblast
'#572E48', #Yolk sac
'#050505', #Parietal Endoderm
'#FF00B7', #DPPA3 TGM3
'#7B00FF', #Dppa3 Hand1
'#00FF15', #Dppa3 low Dusp6 neg
'#D4EB63', #Dusp6 not other
'#54DCE3', #Wfdc2 Aldoc
'#FFFF00', #Gata2 Hand1
'#FFA600' #Lgals3
)
# -

if(do_init) {
    library("metacell")
    scdb_init("scdb_embexe/", force_reinit=T)
    mm_mc = scdb_mc("mm_embexe")
    oc_mc = scdb_mc("oc_embexe")
    mm_legc = log2(1e-5+mm_mc@e_gc)
    oc_legc = log2(1e-5+oc_mc@e_gc)
    mm_lfp = mm_legc - rowMeans(mm_legc)
    oc_lfp = oc_legc - rowMeans(oc_legc)
    mm_color_type = mm_mc@color_key$group
    names(mm_color_type) = mm_mc@color_key$color
    mm_mc_type = mm_color_type[mm_mc@colors]
    oc_color_type = oc_mc@color_key$group
    names(oc_color_type) = oc_mc@color_key$color
    oc_mc_type = oc_color_type[oc_mc@colors]
    oc_mat = scdb_mat("oc_embexe")
    oc_mc_time = tapply(oc_mat@cell_metadata[names(oc_mc@mc),"developmental_time"], oc_mc@mc, mean, na.rm=T)
}

mm_mat = scdb_mat("mm_embexe")

mm_mc_time = tapply(mm_mat@cell_metadata[names(mm_mc@mc),"developmental_time"], mm_mc@mc, mean, na.rm=T)

f_mm_endo = mm_mc_type %in% c("Definitive endoderm","ExE visceral endoderm","Hindgut","Parietal endoderm","Visceral endoderm","Visceral endoderm - anterior")

fG_oc_endo = oc_mc_type %in% c("definitive endoderm","Foregut","Gut","Hindgut","hypoblast","Hypoblast-Anterior (Hhex/Lhx1/Cer1)","Parietal Endoderm","Yolk sac")

f_mm_eendo = mm_mc_type %in% c("ExE visceral endoderm","Parietal endoderm","Visceral endoderm","Visceral endoderm - anterior")

f_oc_eendo = oc_mc_type %in% c("hypoblast","Hypoblast-Anterior (Hhex/Lhx1/Cer1)","Parietal Endoderm","Yolk sac")

f_mm_troph = mm_mc_type %in% c("Chorion", "Chorion progenitors", "EPC progenitors", "p-TGC", "SpA-TGC", "SpT-Gly", "TGC progenitors", "TSC2")

f_oc_troph = oc_mc_type %in% c("Dppa3 Hand1", "Dppa3 low Dusp6 neg", "DPPA3 TGM3", "Dusp6 not other", "Gata2 Hand1", "Lgals3", "Wfdc2 Aldoc")

f_oc_ps = oc_mc_type %in% c("epiblast", "Primitive streak", "Anterior primitive streak")
f_mm_ps = mm_mc_type %in% c("Epiblast", "Primitive streak", "Anterior Primitive Streak")

plot_restric_gg = function(g1, g2, mode = "ee", png_dir=NULL, oc_g1 = toupper(g1), oc_g2 = toupper(g2))
{
	if(mode == "tr") {
		f_oc = f_oc_troph
		f_mm = f_mm_troph
		col_oc = oc_mc@colors[f_oc]
		col_mm = mm_mc@colors[f_mm]
	} else if(mode == "ee") {
		f_oc = f_oc_eendo
		f_mm = f_mm_eendo
		col_oc = oc_mc@colors[f_oc]
		col_mm = mm_mc@colors[f_mm]
	} else {
		f_oc = f_oc_ps
		f_mm = f_mm_ps
		col_oc = oc_mc@colors[f_oc]
		col_oc[col_oc == "#c19f70"] = "#A07AA1"
		col_mm = mm_mc@colors[f_mm]
		col_mm[col_mm == "#C19F70"] = "#A07AA1"
	}
	if(!is.null(png_dir)) {
		png_fn = sprintf("%s/%s.%s.%s.png", png_dir, mode, g1, g2)
		png(png_fn, w=800, h=400);
	}
	layout(matrix(1:2,nrow=1))
	xmin = min(c(oc_legc[oc_g1,f_oc],mm_legc[g1,f_mm]))
	xmax = max(c(oc_legc[oc_g1,f_oc],mm_legc[g1,f_mm]))
	ymin = min(c(oc_legc[oc_g2,f_oc],mm_legc[g2,f_mm]))
	ymax = max(c(oc_legc[oc_g2,f_oc],mm_legc[g2,f_mm]))
	plot(oc_legc[oc_g1,f_oc], oc_legc[oc_g2,f_oc], pch=21, bg=col_oc, cex=1.8, xlim=c(xmin, xmax), ylim=c(ymin, ymax), main="Rabbit", xlab=g1, ylab=g2)
	plot(mm_legc[g1,f_mm], mm_legc[g2,f_mm], pch=21, bg=col_mm, cex=1.8, xlim=c(xmin, xmax), ylim=c(ymin, ymax), main="Mouse", xlab=g1, ylab=g2)
	if(!is.null(png_dir)) {
		dev.off()
	}
}

if(0) {
#     eendo_gs = c("EPAS1", "SULT1E1", "AFP","APOA1", "CER1", "TDGF1", "HHEX","LHX1","CUBN", "DAB2", "FOXA1","FOXA2","EOMES","SOX2", "GSC", "CITED2", "TPM4", "LAMA1", "LAMB1", "MYC", "KLF5","FST", "STRA6", "VIM")
    eendo_gs = c("NOG", "WNT3", "LHX1", "CER1", "HHEX", "LAMB1", "DKK1", "SRGN", "EPAS1", "OTC", "IGF1", "LAMA1", "STMN2", "KLF5", "MYC", "ESRRB", "FOXA1", "CPS1", "SLC3A1", "VIM", "CUBN", "APOM", "APOE", "APOC4", "APOA1", "APOA2", "APOA4", "TTR", "RBP4", "APOB", "AFP", "DAB2", "FN1", "FOXA2", "EOMES", "TRH")
    for(g in eendo_gs) {
        plot_eendo_gg("SOX17", g, png_dir="figs/eendo/")
    }
}

library(tgstat)

grep("APOC", rownames(ee_oc_legc), v=T, ignore.case = T)

max(ee_oc_legc["APOC2",])
max(ee_mm_legc["APOC4",])

bad_names = c('APOC4' )
good_names = c("APOC2" )
names(good_names) = bad_names
mm_genes = rownames(mm_legc)
mm_genes[mm_genes %in% names(good_names)] = good_names[mm_genes[mm_genes %in% names(good_names)]]
rownames(mm_legc) = mm_genes

ee_mm_lfp = ee_mm_legc - rowMeans(ee_mm_legc)

if(do_endo) { #eendoderm
    mm_ee_gmax = apply(mm_legc[,f_mm_eendo], 1, max)
    mm_ee_gmin = apply(mm_legc[,f_mm_eendo], 1, min)
    oc_ee_gmax = apply(oc_legc[,f_oc_eendo], 1, max)
    oc_ee_gmin = apply(oc_legc[,f_oc_eendo], 1, min)
    mm_var_g = names(which(mm_ee_gmax > -13 & (mm_ee_gmax-mm_ee_gmin) > log2(3)))
    mm_var_g_orto = intersect(mm_var_g, rownames(oc_legc))
    oc_var_g = names(which(oc_ee_gmax > -13 & (oc_ee_gmax-oc_ee_gmin) > log2(3)))
    oc_var_g_orto = intersect(oc_var_g, rownames(mm_legc))
    g_ee_union_anchors = union(oc_var_g_orto, mm_var_g_orto)
    g_ee_anchors = intersect(oc_var_g_orto, mm_var_g_orto)
    ee_mm_legc = mm_legc[,f_mm_eendo]
    ee_mm_lfp = ee_mm_legc - rowMeans(ee_mm_legc)
    ee_oc_legc = oc_legc[,f_oc_eendo]
    ee_oc_lfp = ee_oc_legc - rowMeans(ee_oc_legc)
    ee_mm_gcor = tgs_cor(t(ee_mm_lfp[g_ee_anchors,]))
    ee_oc_gcor = tgs_cor(t(ee_oc_lfp[g_ee_anchors,]))
    ee_cross = tgs_cor(ee_oc_lfp[g_ee_anchors,], ee_mm_lfp[g_ee_anchors,])
    hc_ee_oc = hclust(tgs_dist(ee_cross), "ward.D2")
    hc_ee_mm = hclust(tgs_dist(t(ee_cross)), "ward.D2")
    mm_type_rank = as.numeric(as.factor(mm_mc_type[f_mm_eendo]))
    mm_reord = c(3,4,2,1)
    names(mm_reord) = 1:4
    mm_type_rank = mm_reord[mm_type_rank]
    ord_ee_mm = order(1000*mm_type_rank +
                            pmax(colMeans(ee_mm_lfp[c("LAMB1", "KLF5","LAMA1"),])*10,0) +
                            pmax(ee_mm_lfp["EOMES",],0) - pmax(colMeans(ee_mm_lfp[c("DAB2","RBP4"),]) ,0))
    oc_type_rank = as.numeric(as.factor(oc_mc_type[f_oc_eendo]))
    oc_reord = c(1,4,3,2)
    names(oc_reord) = 1:4
    oc_type_rank = oc_reord[oc_type_rank]
    ord_ee_oc = order(1000*oc_type_rank +
                            pmax(colMeans(ee_oc_lfp[c("LAMB1","KLF5","LAMA1"),])*10,0)+ 
                            pmax(ee_oc_lfp["EOMES",],0) - 
                            pmax(colMeans(ee_oc_lfp[c("DAB2","RBP4"),]),0) - 
                            pmax(ifelse(ee_oc_lfp["LAMA1",] < 0, ee_oc_lfp["MYC",]+0.5,0),0)*10)
    n_ee_oc = sum(f_oc_eendo)
    n_ee_mm = sum(f_mm_eendo)
    if(do_pdf) {
        pdf("figs/cross_eendo.pdf", w=10, h=10)
    } else {
        png("figs/cross_eendo.png", w=1000, h=1000)
    }
    layout(matrix(1:4,nrow=2),widths=c(0.9,0.1),heights=c(0.9,0.1))
    par(mar=c(0,2*r_pixel,2*r_pixel,0))
    image(ee_cross[ord_ee_oc, ord_ee_mm], breaks=c(-1,seq(-0.4,0.4,l=999),1),col=fshades, xaxt='n', yaxt='n')
    par(mar=c(2*r_pixel,2*r_pixel,0,0))
    image(matrix((1:n_ee_oc)[ord_ee_oc], ncol=1), breaks=(0:n_ee_oc)+0.5, col=oc_mc@colors[f_oc_eendo], xaxt='n', yaxt='n')
    par(mar=c(0,0,2*r_pixel,2*r_pixel))
    image(matrix((1:n_ee_mm)[ord_ee_mm], nrow=1), breaks=(0:n_ee_mm)+0.5, col=mm_mc@colors[f_mm_eendo], xaxt='n', yaxt='n')
    dev.off()
    eendo_gs = c("NOG", "WNT3", "LHX1", "CER1", "HHEX", "LAMB1", "DKK1", "SRGN", "EPAS1", "OTC", "IGF1", "LAMA1", "STMN2", "KLF5", "MYC", "ESRRB", "FOXA1", "CPS1", "SLC3A1", "VIM", "CUBN", "APOM", "APOE", "APOC2", "APOA1", "APOA2", "APOA4", "TTR", "RBP4", "APOB", "AFP", "DAB2", "FN1", "FOXA2", "EOMES", "TRH")
    #eendo_gs = c("EPAS1","AFP","APOA1", "CER1", "TDGF1", "HHEX","LHX1","CUBN", "DAB2", "FOXA1","FOXA2","EOMES","GSC", "CITED2", "TPM4", "LAMA1", "LAMB1", "MYC", "KLF5","FST", "STRA6", "RBP4","ESRRB", "TTR", "VIM")
    n_g = length(eendo_gs)
    gs_hc = hclust(tgs_dist(tgs_cor(t(cbind(ee_oc_lfp[eendo_gs,],ee_mm_lfp[eendo_gs,]))),"ward.D2"))
    eendo_gs = eendo_gs[gs_hc$order]
    if(do_pdf) {
        pdf("figs/oc_eendo_marks.pdf", w=12, h=8)
    } else {
        png("figs/oc_eendo_marks.png", w=1200, h=800)
    }
    layout(matrix(1:2,nrow=2),heights=c(0.9,0.1))
    par(mar=c(0,8,2,8))
    #image(t(ee_oc_lfp[eendo_gs, ord_ee_oc]), col=rshades, breaks=c(-10,seq(-2,2,l=999),10), xaxt='n', yaxt='n')
    image(t(ee_oc_legc[eendo_gs, ord_ee_oc]), col=rshades, breaks=c(-17,seq(-15,-9,l=999),0), xaxt='n', yaxt='n')
    mtext(eendo_gs, at=seq(0,1,l=n_g),las=2,side=2, font=6)
    mtext(eendo_gs, at=seq(0,1,l=n_g),las=2,side=4, font=6)
    par(mar=c(2,8,0,8))
    image(matrix((1:n_ee_oc)[ord_ee_oc], ncol=1), breaks=(0:n_ee_oc)+0.5, col=oc_mc@colors[f_oc_eendo], xaxt='n', yaxt='n')
    dev.off()
    if(do_pdf) {
        pdf("figs/mm_eendo_marks.pdf", w=12, h=8)
    } else {
        png("figs/mm_eendo_marks.png", w=1200, h=800)
    }
    layout(matrix(1:2,nrow=2),heights=c(0.9,0.1))
    par(mar=c(0,8,2,8))
    #image(t(ee_mm_lfp[eendo_gs, ord_ee_mm]), col=rshades, breaks=c(-10,seq(-2,2,l=999),10), xaxt='n', yaxt='n')
    image(t(ee_mm_legc[eendo_gs, ord_ee_mm]), col=rshades, breaks=c(-17,seq(-15, -9,l=999),0), xaxt='n', yaxt='n')
    mtext(eendo_gs, at=seq(0,1,l=n_g),las=2,side=2, font=6)
    mtext(eendo_gs, at=seq(0,1,l=n_g),las=2,side=4, font=6)
    par(mar=c(2,8,0,8))
    image(matrix((1:n_ee_mm)[ord_ee_mm], ncol=1), breaks=(0:n_ee_mm)+0.5, col=mm_mc@colors[f_mm_eendo], xaxt='n', yaxt='n')
    dev.off()
}

library(matrixStats)
zscore = oc_mc@e_gc[eendo_gs, ord_ee_oc]
zscore = zscore - rowMeans(zscore)
zscore = zscore / rowSds(zscore)


n_ee_oc = length(genes_order)
n_ee_mm = length(genes_order)

layout(matrix(1:2,nrow=2),heights=c(0.9,0.1))
par(mar=c(0,8,2,8))
# image(t(zscore[eendo_gs, ord_ee_oc]), col=rshades, breaks=c(-10,seq(-2,2,l=999),10), xaxt='n', yaxt='n')
image(t(ee_oc_legc[genes_order, ord_ee_oc]), col=rshades, breaks=c(-17,seq(-15,-9,l=999),0), xaxt='n', yaxt='n')
mtext(eendo_gs, at=seq(0,1,l=n_g),las=2,side=2, font=6)
mtext(eendo_gs, at=seq(0,1,l=n_g),las=2,side=4, font=6)
par(mar=c(2,8,0,8))
image(matrix((1:n_ee_oc)[ord_ee_oc], ncol=1), breaks=(0:n_ee_oc)+0.5, col=oc_mc@colors[f_oc_eendo], xaxt='n', yaxt='n')

    eendo_gs = genes_order
    if(do_pdf) {
        pdf("figs/oc_eendo_marks.pdf", w=12, h=8)
    } else {
        png("figs/oc_eendo_marks.png", w=1200, h=800)
    }
    layout(matrix(1:2,nrow=2),heights=c(0.9,0.1))
    par(mar=c(0,8,2,8))
#     image(t(zscore[eendo_gs, ord_ee_oc]), col=rshades, breaks=c(-10,seq(-2,2,l=999),10), xaxt='n', yaxt='n')
    image(t(ee_oc_legc[genes_order, ord_ee_oc]), col=rshades, breaks=c(-17,seq(-15,-9,l=999),0), xaxt='n', yaxt='n')
    mtext(eendo_gs, at=seq(0,1,l=n_g),las=2,side=2, font=6)
    mtext(eendo_gs, at=seq(0,1,l=n_g),las=2,side=4, font=6)
    par(mar=c(2,8,0,8))
    image(matrix((1:n_ee_oc)[ord_ee_oc], ncol=1), breaks=(0:n_ee_oc)+0.5, col=oc_mc@colors[f_oc_eendo], xaxt='n', yaxt='n')
    dev.off()
    if(do_pdf) {
        pdf("figs/mm_eendo_marks.pdf", w=12, h=8)
    } else {
        png("figs/mm_eendo_marks.png", w=1200, h=800)
    }
    layout(matrix(1:2,nrow=2),heights=c(0.9,0.1))
    par(mar=c(0,8,2,8))
#     image(t(ee_mm_lfp[eendo_gs, ord_ee_mm]), col=rshades, breaks=c(-10,seq(-2,2,l=999),10), xaxt='n', yaxt='n')
    image(t(ee_mm_legc[genes_order, ord_ee_mm]), col=rshades, breaks=c(-17,seq(-15, -9,l=999),0), xaxt='n', yaxt='n')
    mtext(eendo_gs, at=seq(0,1,l=n_g),las=2,side=2, font=6)
    mtext(eendo_gs, at=seq(0,1,l=n_g),las=2,side=4, font=6)
    par(mar=c(2,8,0,8))
    image(matrix((1:n_ee_mm)[ord_ee_mm], ncol=1), breaks=(0:n_ee_mm)+0.5, col=mm_mc@colors[f_mm_eendo], xaxt='n', yaxt='n')
    dev.off()

if(do_troph) { #e ectoderm
    mm_tr_gmax = apply(mm_legc[,f_mm_troph], 1, max)
    mm_tr_gmin = apply(mm_legc[,f_mm_troph], 1, min)
    oc_tr_gmax = apply(oc_legc[,f_oc_troph], 1, max)
    oc_tr_gmin = apply(oc_legc[,f_oc_troph], 1, min)
    mm_var_g = names(which(mm_tr_gmax > -13 & (mm_tr_gmax-mm_tr_gmin) > log2(4)))
    mm_var_g_orto = intersect(mm_var_g, rownames(oc_legc))
    oc_var_g = names(which(oc_tr_gmax > -13 & (oc_tr_gmax-oc_tr_gmin) > log2(4)))
    oc_var_g_orto = intersect(oc_var_g, rownames(mm_legc))
    g_tr_union_anchors = union(oc_var_g_orto, mm_var_g_orto)
    g_tr_anchors = intersect(oc_var_g_orto, mm_var_g_orto)
    tr_mm_lfp = mm_legc[,f_mm_troph]
    tr_mm_lfp = tr_mm_lfp - rowMeans(tr_mm_lfp)
    tr_oc_lfp = oc_legc[,f_oc_troph]
    tr_oc_lfp = tr_oc_lfp - rowMeans(tr_oc_lfp)
    tr_mm_gcor = tgs_cor(t(tr_mm_lfp[g_tr_anchors,]))
    tr_oc_gcor = tgs_cor(t(tr_oc_lfp[g_tr_anchors,]))
    tr_cross = tgs_cor(tr_oc_lfp[g_tr_anchors,], tr_mm_lfp[g_tr_anchors,])
#     hc_tr_oc = hclust(tgs_dist(tr_cross), "ward.D2")
#     hc_tr_mm = hclust(tgs_dist(t(tr_cross)), "ward.D2")
    
    mm_type_rank = as.numeric(as.factor(mm_mc_type[f_mm_troph]))
    mm_reord = c(2,1,5,6,4,7,3)
    names(mm_reord) = 1:7
    mm_type_rank = mm_reord[mm_type_rank]
#     mm_type_rank[mm_type_rank == 2] = 0
    ord_ee_mm = order(1000*mm_type_rank +
                            mm_mc_time[f_oc_troph])
    oc_type_rank = as.numeric(as.factor(oc_mc_type[f_oc_troph]))
#     oc_type_rank[oc_type_rank == 2] = 5
    ord_ee_oc = order(1000*oc_type_rank +
                            oc_mc_time[f_mm_troph])
    n_tr_oc = sum(f_oc_troph)
    n_tr_mm = sum(f_mm_troph)
    if(do_pdf) {
        pdf("figs/cross_troph.pdf", w=10, h=10)
    } else {
        png("figs/cross_troph.png", w=1000, h=1000)
    }
    layout(matrix(1:4,nrow=2),widths=c(0.9,0.1),heights=c(0.9,0.1))
    par(mar=c(0,2,2,0)*r_pixel)
    image(tr_cross[ord_ee_oc, ord_ee_mm], breaks=c(-1,seq(-0.4,0.4,l=999),1),col=fshades, xaxt='n', yaxt='n')
    par(mar=c(2,2,0,0)*r_pixel)
    image(matrix((1:n_tr_oc)[ord_ee_oc], ncol=1), breaks=(0:n_tr_oc)+0.5, col=oc_mc@colors[f_oc_troph], xaxt='n', yaxt='n')
    par(mar=c(0,0,2,2)*r_pixel)
    image(matrix((1:n_tr_mm)[ord_ee_mm], nrow=1), breaks=(0:n_tr_mm)+0.5, col=mm_mc@colors[f_mm_troph], xaxt='n', yaxt='n')
    dev.off()
    troph_gs = c("SOX2", "ESRRB", "DPPA3", "EOMES", "GATA2", "HAND1", "FOSL1", "ADM","BMP4", "DUSP6", "TGM3", "MSX1", "DLX3", "CAV3","AQP3", "CDH5")
    n_g = length(troph_gs)
    gs_hc = hclust(tgs_dist(tgs_cor(t(cbind(tr_oc_lfp[troph_gs,],tr_mm_lfp[troph_gs,]))),"ward.D2"))
    troph_gs = troph_gs[gs_hc$order]
    if(do_pdf) {
        pdf("figs/oc_troph_marks.pdf", w=12, h=8)
    } else {
        png("figs/oc_troph_marks.png", w=1200, h=800)
    }
    layout(matrix(1:2,nrow=2),heights=c(0.9,0.1))
    par(mar=c(0,8,2,8)*r_pixel)
    image(t(tr_oc_legc[troph_gs, ord_ee_oc]), col=rshades, breaks=c(-17,seq(-16.5, -10,l=999),0), xaxt='n', yaxt='n')
    mtext(troph_gs, at=seq(0,1,l=n_g),las=2,side=2, font=6)
    mtext(troph_gs, at=seq(0,1,l=n_g),las=2,side=4, font=6)
    par(mar=c(2,8,0,8)*r_pixel)
    image(matrix((1:n_tr_oc)[ord_ee_oc], ncol=1), breaks=(0:n_tr_oc)+0.5, col=oc_mc@colors[f_oc_troph], xaxt='n', yaxt='n')
    dev.off()
    if(do_pdf) {
    pdf("figs/mm_troph_marks.pdf", w=12, h=8)
    } else {
    png("figs/mm_troph_marks.png", w=1200, h=800)
    }
    layout(matrix(1:2,nrow=2),heights=c(0.9,0.1))
    par(mar=c(0,8,2,8)*r_pixel)     
    tr_mm_legc = mm_legc[,f_mm_troph]
    image(t(tr_mm_legc[troph_gs, ord_ee_mm]), col=rshades, breaks=c(-17,seq(-16.5, -10,l=999),0), xaxt='n', yaxt='n')
    mtext(troph_gs, at=seq(0,1,l=n_g),las=2,side=2, font=6)
    mtext(troph_gs, at=seq(0,1,l=n_g),las=2,side=4, font=6)
    par(mar=c(2,8,0,8)*r_pixel)
    image(matrix((1:n_tr_mm)[ord_ee_mm], ncol=1), breaks=(0:n_tr_mm)+0.5, col=mm_mc@colors[f_mm_troph], xaxt='n', yaxt='n')
    dev.off()
    tr_oc_legc = oc_legc[,f_oc_troph]
    gcor_time_oc_tr = apply(tr_oc_lfp, 1, cor, oc_mc_time[f_oc_troph], use="pairwise.complete.obs")
    oc_tr_g_time_h = names(tail(sort(gcor_time_oc_tr),30))
    oc_tr_g_time_l = names(head(sort(gcor_time_oc_tr),30))
    for(g in c(oc_tr_g_time_h)) {
        if(do_pdf) {
        pdf(sprintf("figs/troph/time_h/%s.pdf", g),w=4, h=4)
        } else {
        png(sprintf("figs/troph/time_h/%s.png", g),w=400,h=400)
        }
        plot(oc_mc_time[f_oc_troph], tr_oc_legc[g,], pch=19, col=oc_mc@colors[f_oc_troph], main=g)
        dev.off()
    }
    for(g in c(oc_tr_g_time_l, "KLF5", "GDF3","DUSP4", "DPPA4")) {
        if(do_pdf) {
        pdf(sprintf("figs/troph/time_l/%s.pdf", g),w=4,h=4)
        } else {
        png(sprintf("figs/troph/time_l/%s.png", g),w=400,h=400)
        }
        plot(oc_mc_time[f_oc_troph], tr_oc_legc[g,], pch=19, col=oc_mc@colors[f_oc_troph], main=g)
        dev.off()
    }
    f_oc_tr_late = oc_mc_time[f_oc_troph]>7.7
    gcor_gata2_oc_late_tr = apply(tr_oc_lfp[,f_oc_tr_late], 1, cor, 
                        tr_oc_lfp["GATA2",f_oc_tr_late], 
                        use="pairwise.complete.obs")
    oc_tr_g_late_h = names(tail(sort(gcor_gata2_oc_late_tr),30))
    oc_tr_g_late_l = names(head(sort(gcor_gata2_oc_late_tr),30))
    for(g in c(oc_tr_g_late_h)) {
        if(do_pdf) {
        pdf(sprintf("figs/troph/grad_h/%s.pdf", g),w=4,h=4)
        } else {
        png(sprintf("figs/troph/grad_h/%s.png", g),w=400,h=400)
        }
        plot(oc_mc_time[f_oc_troph], tr_oc_legc[g,], pch=19, col=oc_mc@colors[f_oc_troph], main=g)
        dev.off()
    }
    for(g in c(oc_tr_g_late_l)) {
        if(do_pdf) {
        pdf(sprintf("figs/troph/grad_l/%s.pdf", g),w=4,h=4)
        } else {
        png(sprintf("figs/troph/grad_l/%s.png", g),w=400,h=400)
        }
        plot(oc_mc_time[f_oc_troph], tr_oc_legc[g,], pch=19, col=oc_mc@colors[f_oc_troph], main=g)
        dev.off()
    }
    f_tr_late_gata2_h = f_oc_troph & oc_legc["GATA2",]> -11 & oc_mc_time > 7.7
    f_tr_late_gata2_m = f_oc_troph & oc_legc["GATA2",] > -12.8 &
                                                oc_mc_time > 7.7 & !f_tr_late_gata2_h
    f_tr_late_gata2_l = f_oc_troph & oc_legc["GATA2",] < -12.8 & oc_mc_time > 7.7 & !is.na(oc_mc_time)
    lfp_tr_lgata2_h = rowMeans(oc_lfp[,f_tr_late_gata2_h])
    lfp_tr_lgata2_m = rowMeans(oc_lfp[,f_tr_late_gata2_m])
    lfp_tr_lgata2_l = rowMeans(oc_lfp[,f_tr_late_gata2_l])
    lfp_tr_mm_chorion= rowMeans(mm_lfp[,mm_mc_type=="Chorion"])
    lfp_tr_mm_tgcp= rowMeans(mm_lfp[,mm_mc_type=="TGC progenitors"])
    lfp_tr_mm_epc = rowMeans(mm_lfp[,mm_mc_type=="EPC progenitors"])
    lfp_tr_mm_spt = rowMeans(mm_lfp[,mm_mc_type=="SpT-Gly"])
    g_orto = intersect(names(lfp_tr_mm_chorion), names(lfp_tr_lgata2_h))
    top_n = 30
    top20_chorion = names(tail(sort(lfp_tr_mm_chorion[g_orto]),top_n))
    top20_epc= names(tail(sort(lfp_tr_mm_epc[g_orto]),top_n))
    top20_tgcp = names(tail(sort(lfp_tr_mm_tgcp[g_orto]),top_n))
    top20_spt= names(tail(sort(lfp_tr_mm_spt[g_orto]),top_n))
    top20_lgata2_h= names(tail(sort(lfp_tr_lgata2_h[g_orto]),top_n))
    top20_lgata2_m= names(tail(sort(lfp_tr_lgata2_m[g_orto]),top_n))
    top20_lgata2_l= names(tail(sort(lfp_tr_lgata2_l[g_orto]),top_n))
    df = data.frame(lgata2_h = lfp_tr_lgata2_h[g_orto], 
                         lgata2_m = lfp_tr_lgata2_m[g_orto], 
                        lgata2_l = lfp_tr_lgata2_l[g_orto], 
                         chorion = lfp_tr_mm_chorion[g_orto], 
                         tgcp = lfp_tr_mm_tgcp[g_orto], 
                         spt = lfp_tr_mm_spt[g_orto], 
                         epc = lfp_tr_mm_epc[g_orto])
    rownames(df) = g_orto
    df$dlt_chor = df$chorion - rowMeans(df[,1:3])
    gs = union(union(top20_chorion,top20_lgata2_h), union(top20_lgata2_m, top20_lgata2_l))
    gs = union(gs, top20_epc)
    gs = union(gs, top20_tgcp)
    gs = union(gs, top20_spt)
    pheatmap::pheatmap(df[gs,1:7], col=rshades, 
                breaks=c(-6,seq(-3,3,l=999),6), fontsize=8, w=5,h=20, filename="figs/troph_marks.pdf")
}

for(g in c("IGF1")) {
    if(do_pdf) {
    pdf(sprintf("figs/troph/%s.pdf", g),w=4, h=4)
    } else {
    png(sprintf("figs/troph/%s.png", g),w=400,h=400)
    }
    plot(oc_mc_time[f_oc_troph], tr_oc_legc[g,], pch=19, col=oc_mc@colors[f_oc_troph], main=g)
    dev.off()
}

    mm_type_rank = as.numeric(as.factor(mm_mc_type[f_mm_troph]))
    mm_reord = c(2,1,5,6,4,7,3)
    names(mm_reord) = 1:7
    mm_type_rank = mm_reord[mm_type_rank]
#     mm_type_rank[mm_type_rank == 2] = 0
    ord_ee_mm = order(1000*mm_type_rank +
                            mm_mc_time[f_mm_troph])
    oc_type_rank = as.numeric(as.factor(oc_mc_type[f_oc_troph]))
#     oc_type_rank[oc_type_rank == 2] = 5
    ord_ee_oc = order(1000*oc_type_rank +
                            oc_mc_time[f_oc_troph])
    n_tr_oc = sum(f_oc_troph)
    n_tr_mm = sum(f_mm_troph)
    if(do_pdf) {
        pdf("figs/cross_troph.pdf", w=10, h=10)
    } else {
        png("figs/cross_troph.png", w=1000, h=1000)
    }
    layout(matrix(1:4,nrow=2),widths=c(0.9,0.1),heights=c(0.9,0.1))
    par(mar=c(0,2,2,0)*r_pixel)
    image(tr_cross[ord_ee_oc, ord_ee_mm], breaks=c(-1,seq(-0.4,0.4,l=999),1),col=fshades, xaxt='n', yaxt='n')
    par(mar=c(2,2,0,0)*r_pixel)
    image(matrix((1:n_tr_oc)[ord_ee_oc], ncol=1), breaks=(0:n_tr_oc)+0.5, col=oc_mc@colors[f_oc_troph], xaxt='n', yaxt='n')
    par(mar=c(0,0,2,2)*r_pixel)
    image(matrix((1:n_tr_mm)[ord_ee_mm], nrow=1), breaks=(0:n_tr_mm)+0.5, col=mm_mc@colors[f_mm_troph], xaxt='n', yaxt='n')
    dev.off()
    troph_gs = c("SOX2", "ESRRB", "DPPA3", "EOMES", "GATA2", "HAND1", "FOSL1", "ADM","BMP4", "DUSP6", "TGM3", "MSX1", "DLX3", "CAV3","AQP3", "CDH5")
    n_g = length(troph_gs)
    gs_hc = hclust(tgs_dist(tgs_cor(t(cbind(tr_oc_lfp[troph_gs,],tr_mm_lfp[troph_gs,]))),"ward.D2"))
    troph_gs = troph_gs[gs_hc$order]
    if(do_pdf) {
    pdf("figs/oc_troph_marks.pdf", w=12, h=8)
    } else {
    png("figs/oc_troph_marks.png", w=1200, h=800)
    }
    layout(matrix(1:2,nrow=2),heights=c(0.9,0.1))
    par(mar=c(0,8,2,8)*r_pixel)
    image(t(tr_oc_legc[troph_gs, ord_ee_oc]), col=rshades, breaks=c(-17,seq(-16.5, -10,l=999),0), xaxt='n', yaxt='n')
    mtext(troph_gs, at=seq(0,1,l=n_g),las=2,side=2, font=6)
    mtext(troph_gs, at=seq(0,1,l=n_g),las=2,side=4, font=6)
    par(mar=c(2,8,0,8)*r_pixel)
    image(matrix((1:n_tr_oc)[ord_ee_oc], ncol=1), breaks=(0:n_tr_oc)+0.5, col=oc_mc@colors[f_oc_troph], xaxt='n', yaxt='n')
    dev.off()
    if(do_pdf) {
    pdf("figs/mm_troph_marks.pdf", w=12, h=8)
    } else {
    png("figs/mm_troph_marks.png", w=1200, h=800)
    }
    layout(matrix(1:2,nrow=2),heights=c(0.9,0.1))
    par(mar=c(0,8,2,8)*r_pixel)
    tr_mm_legc = mm_legc[,f_mm_troph]
    image(t(tr_mm_legc[troph_gs, ord_ee_mm]), col=rshades, breaks=c(-17,seq(-16.5, -10,l=999),0), xaxt='n', yaxt='n')
    mtext(troph_gs, at=seq(0,1,l=n_g),las=2,side=2, font=6)
    mtext(troph_gs, at=seq(0,1,l=n_g),las=2,side=4, font=6)
    par(mar=c(2,8,0,8)*r_pixel)
    image(matrix((1:n_tr_mm)[ord_ee_mm], ncol=1), breaks=(0:n_tr_mm)+0.5, col=mm_mc@colors[f_mm_troph], xaxt='n', yaxt='n')
    dev.off()
    tr_oc_legc = oc_legc[,f_oc_troph]
    gcor_time_oc_tr = apply(tr_oc_lfp, 1, cor, oc_mc_time[f_oc_troph], use="pairwise.complete.obs")
    oc_tr_g_time_h = names(tail(sort(gcor_time_oc_tr),30))
    oc_tr_g_time_l = names(head(sort(gcor_time_oc_tr),30))
    for(g in c(oc_tr_g_time_h)) {
        if(do_pdf) {
        pdf(sprintf("figs/troph/time_h/%s.pdf", g),w=4, h=4)
        } else {
        png(sprintf("figs/troph/time_h/%s.png", g),w=400,h=400)
        }
        plot(oc_mc_time[f_oc_troph], tr_oc_legc[g,], pch=19, col=oc_mc@colors[f_oc_troph], main=g)
        dev.off()
    }
    for(g in c(oc_tr_g_time_l, "KLF5", "GDF3","DUSP4", "DPPA4")) {
        if(do_pdf) {
        pdf(sprintf("figs/troph/time_l/%s.pdf", g),w=4,h=4)
        } else {
        png(sprintf("figs/troph/time_l/%s.png", g),w=400,h=400)
        }
        plot(oc_mc_time[f_oc_troph], tr_oc_legc[g,], pch=19, col=oc_mc@colors[f_oc_troph], main=g)
        dev.off()
    }
    f_oc_tr_late = oc_mc_time[f_oc_troph]>7.7
    gcor_gata2_oc_late_tr = apply(tr_oc_lfp[,f_oc_tr_late], 1, cor, 
                        tr_oc_lfp["GATA2",f_oc_tr_late], 
                        use="pairwise.complete.obs")
    oc_tr_g_late_h = names(tail(sort(gcor_gata2_oc_late_tr),30))
    oc_tr_g_late_l = names(head(sort(gcor_gata2_oc_late_tr),30))
    for(g in c(oc_tr_g_late_h)) {
        if(do_pdf) {
        pdf(sprintf("figs/troph/grad_h/%s.pdf", g),w=4,h=4)
        } else {
        png(sprintf("figs/troph/grad_h/%s.png", g),w=400,h=400)
        }
        plot(oc_mc_time[f_oc_troph], tr_oc_legc[g,], pch=19, col=oc_mc@colors[f_oc_troph], main=g)
        dev.off()
    }
    for(g in c(oc_tr_g_late_l)) {
        if(do_pdf) {
        pdf(sprintf("figs/troph/grad_l/%s.pdf", g),w=4,h=4)
        } else {
        png(sprintf("figs/troph/grad_l/%s.png", g),w=400,h=400)
        }
        plot(oc_mc_time[f_oc_troph], tr_oc_legc[g,], pch=19, col=oc_mc@colors[f_oc_troph], main=g)
        dev.off()
    }
    f_tr_late_gata2_h = f_oc_troph & oc_legc["GATA2",]> -11 & oc_mc_time > 7.7
    f_tr_late_gata2_m = f_oc_troph & oc_legc["GATA2",] > -12.8 &
                                                oc_mc_time > 7.7 & !f_tr_late_gata2_h
    f_tr_late_gata2_l = f_oc_troph & oc_legc["GATA2",] < -12.8 & oc_mc_time > 7.7 & !is.na(oc_mc_time)
    lfp_tr_lgata2_h = rowMeans(oc_lfp[,f_tr_late_gata2_h])
    lfp_tr_lgata2_m = rowMeans(oc_lfp[,f_tr_late_gata2_m])
    lfp_tr_lgata2_l = rowMeans(oc_lfp[,f_tr_late_gata2_l])
    lfp_tr_mm_chorion= rowMeans(mm_lfp[,mm_mc_type=="Chorion"])
    lfp_tr_mm_tgcp= rowMeans(mm_lfp[,mm_mc_type=="TGC progenitors"])
    lfp_tr_mm_epc = rowMeans(mm_lfp[,mm_mc_type=="EPC progenitors"])
    lfp_tr_mm_spt = rowMeans(mm_lfp[,mm_mc_type=="SpT-Gly"])
    g_orto = intersect(names(lfp_tr_mm_chorion), names(lfp_tr_lgata2_h))
    top_n = 30
    top20_chorion = names(tail(sort(lfp_tr_mm_chorion[g_orto]),top_n))
    top20_epc= names(tail(sort(lfp_tr_mm_epc[g_orto]),top_n))
    top20_tgcp = names(tail(sort(lfp_tr_mm_tgcp[g_orto]),top_n))
    top20_spt= names(tail(sort(lfp_tr_mm_spt[g_orto]),top_n))
    top20_lgata2_h= names(tail(sort(lfp_tr_lgata2_h[g_orto]),top_n))
    top20_lgata2_m= names(tail(sort(lfp_tr_lgata2_m[g_orto]),top_n))
    top20_lgata2_l= names(tail(sort(lfp_tr_lgata2_l[g_orto]),top_n))
    df = data.frame(lgata2_h = lfp_tr_lgata2_h[g_orto], 
                         lgata2_m = lfp_tr_lgata2_m[g_orto], 
                        lgata2_l = lfp_tr_lgata2_l[g_orto], 
                         chorion = lfp_tr_mm_chorion[g_orto], 
                         tgcp = lfp_tr_mm_tgcp[g_orto], 
                         spt = lfp_tr_mm_spt[g_orto], 
                         epc = lfp_tr_mm_epc[g_orto])
    rownames(df) = g_orto
    df$dlt_chor = df$chorion - rowMeans(df[,1:3])
    gs = union(union(top20_chorion,top20_lgata2_h), union(top20_lgata2_m, top20_lgata2_l))
    gs = union(gs, top20_epc)
    gs = union(gs, top20_tgcp)
    gs = union(gs, top20_spt)
    pheatmap::pheatmap(df[gs,1:7], col=rshades, 
                breaks=c(-6,seq(-3,3,l=999),6), fontsize=8, w=5,h=20, filename="figs/troph_marks.pdf")

library(plot3D)

pdf("eendo_key.pdf")
colkey(col = rshades, breaks = c(-17,seq(-15, -9,l=999),0))
dev.off()


