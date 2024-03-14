set.seed(1234)
library(dplyr)
library(Seurat)
library(SCENIC)
library(AUCell)
library(circlize)
library(ComplexHeatmap)
library(stringr)
library(ggplot2)
library(dendextend)
library(msigdbr)
library(RColorBrewer)
library('stats')
library(tidyHeatmap)
library(colorspace)
library(randomcoloR)
library(colorspace)
library(dendextend)
library(cowplot)
library(scCustomize)
#library(export)
library(ggpubr)
library(readxl)


# Figure 1A ####
excel = read_xlsx('~/Documents/samir_macbook/TMA/TMA97_MSKCC_v9.xlsx', sheet = 1) %>% as.data.frame()
df <- excel[,colSums(is.na(excel))<nrow(excel)]
df = distinct(df)
rownames(df) = paste0('random_',c(1:131))

meta = df %>% dplyr::select( "SoftTissue_Bone"   , "Morphology...26", "Identifier" )
# meta = df[,1:4]; rownames(meta) = rownames(df)
counts = df[,c(6:23)]; rownames(counts) = rownames(df)
counts_nomki = counts[,-12]
counts_mki = counts[,12] %>% as.data.frame(); rownames(counts_mki) = rownames(df)
counts_nomki = t(counts_nomki)
counts = t(counts)

colnames(meta)[1] = 'Site'
meta$Site %>% unique
# meta$Site [meta$Site %>% is.na()] = 'Information not available'
meta$Site = factor(meta$Site, levels = c('Bone','Soft Tissue'))
colnames(meta)[2] = 'Histology'
meta$Histology %>% unique %>% rev()
colnames(meta)[3] = 'Patient'
rownames(meta)
meta$Patient %>% unique %>% length()

counts_nomki %>% rownames()
matrix = matrix(0, nrow = 2, ncol = dim(counts_nomki)[2]) %>% as.data.frame()

counts_nomki = as.data.frame(counts_nomki)
counts_nomki = counts_nomki[c( "AR"   ,   "NKX3.1"  ,"CK8"   ,  "TP63" ,"SYP","INSM1",'ASCL1','NEUROD1','FOXA2',"YAP1",'POU2F3','MYC','SOX2','EZH2','TFF3','TROP2','DLL3'),]


## Make dendrogram for each patient ####
#### p1 ####
m_p1 = counts[,  rownames(meta)[meta$Patient == "Patient 1"], drop = FALSE]
hc1 = hclust(dist(t(m_p1)), method = "ward.D2")
dend_p1 = as.dendrogram(hc1)
order_p1 = m_p1[,order.dendrogram(dend_p1)] %>% colnames()
saveRDS(dend_p1,'20231206_dend_p1_new.rds')

#### p2 ####
m_p2 = counts[,  rownames(meta)[meta$Patient == "Patient 2"], drop = FALSE]
hc2 = hclust(dist(t(m_p2)), method = "ward.D2")
dend_p2 = as.dendrogram(hc2)
order_p2 = m_p2[,order.dendrogram(dend_p2)] %>% colnames()
saveRDS(dend_p2,'20231206_dend_p2_new.rds')

#### p3 ####
m_p3 = counts[,  rownames(meta)[meta$Patient == "Patient 3"], drop = FALSE]
hc3 = hclust(dist(t(m_p3)), method = "ward.D2")
dend_p3 = as.dendrogram(hc3)
order_p3 = m_p3[,order.dendrogram(dend_p3)] %>% colnames()
saveRDS(dend_p3,'20231206_dend_p3_new.rds')

#### p4 ####
m_p4 = counts[,  rownames(meta)[meta$Patient == "Patient 4"], drop = FALSE]
hc4 = hclust(dist(t(m_p4)), method = "ward.D2")
dend_p4 = as.dendrogram(hc4)
order_p4 = m_p4[,order.dendrogram(dend_p4)] %>% colnames()
saveRDS(dend_p4,'20231206_dend_p4_new.rds')

#### p5 ####
m_p5 = counts[,  rownames(meta)[meta$Patient == "Patient 5"], drop = FALSE]
hc5 = hclust(dist(t(m_p5)), method = "ward.D2")
dend_p5 = as.dendrogram(hc5)
order_p5 = m_p5[,order.dendrogram(dend_p5)] %>% colnames()
saveRDS(dend_p5,'20231206_dend_p5_new.rds')

#### p6 ####
m_p6 = counts[,  rownames(meta)[meta$Patient == "Patient 6"], drop = FALSE]
hc6 = hclust(dist(t(m_p6)), method = "ward.D2")
dend_p6 = as.dendrogram(hc6)
order_p6 = m_p6[,order.dendrogram(dend_p6)] %>% colnames()
saveRDS(dend_p6,'20231206_dend_p6_new.rds')

#### p7 ####
m_p7 = counts[,  rownames(meta)[meta$Patient == "Patient 7"], drop = FALSE]
hc7 = hclust(dist(t(m_p7)), method = "ward.D2")
dend_p7 = as.dendrogram(hc7)
order_p7 = m_p7[,order.dendrogram(dend_p7)] %>% colnames()
saveRDS(dend_p6,'20231206_dend_p7_new.rds')

#### p8 ####
m_p8 = counts[,  rownames(meta)[meta$Patient == "Patient 8"], drop = FALSE]
hc8 = hclust(dist(t(m_p8)), method = "ward.D2")
dend_p8 = as.dendrogram(hc8)
order_p8 = m_p8[,order.dendrogram(dend_p8)] %>% colnames()
saveRDS(dend_p8,'20231206_dend_p8_new.rds')

#### p9 ####
m_p9 = counts[,  rownames(meta)[meta$Patient == "Patient 9"], drop = FALSE]
hc9 = hclust(dist(t(m_p9)), method = "ward.D2")
dend_p9 = as.dendrogram(hc9)
order_p9 = m_p9[,order.dendrogram(dend_p9)] %>% colnames()
saveRDS(dend_p9,'20231206_dend_p9_new.rds')

#### p10 ####
m_p10 = counts[,  rownames(meta)[meta$Patient == "Patient 10"], drop = FALSE]
hc10 = hclust(dist(t(m_p10)), method = "ward.D2")
dend_p10 = as.dendrogram(hc10)
order_p10 = m_p10[,order.dendrogram(dend_p10)] %>% colnames()
saveRDS(dend_p10,'20231206_dend_p10_new.rds')

#### p11 ####
m_p11 = counts[,  rownames(meta)[meta$Patient == "Patient 11"], drop = FALSE]
hc11 = hclust(dist(t(m_p11)), method = "ward.D2")
dend_p11 = as.dendrogram(hc11)
order_p11 = m_p11[,order.dendrogram(dend_p11)] %>% colnames()
saveRDS(dend_p11,'20231206_dend_p11_new.rds')

#### p12 ####
m_p12 = counts[,  rownames(meta)[meta$Patient == "Patient 12"], drop = FALSE]
hc12 = hclust(dist(t(m_p12)), method = "ward.D2")
dend_p12 = as.dendrogram(hc12)
order_p12 = m_p12[,order.dendrogram(dend_p12)] %>% colnames()
saveRDS(dend_p12,'20231206_dend_p12_new.rds')

#### p13 ####
m_p13 = counts[,  rownames(meta)[meta$Patient == "Patient 13"], drop = FALSE]
hc13 = hclust(dist(t(m_p13)), method = "ward.D2")
dend_p13 = as.dendrogram(hc13)
order_p13 = m_p13[,order.dendrogram(dend_p13)] %>% colnames()
saveRDS(dend_p13,'20231206_dend_p13_new.rds')

#### p14 ####
m_p14 = counts[,  rownames(meta)[meta$Patient == "Patient 14"], drop = FALSE]
hc14 = hclust(dist(t(m_p14)), method = "ward.D2")
dend_p14 = as.dendrogram(hc14)
order_p14 = m_p14[,order.dendrogram(dend_p14)] %>% colnames()
saveRDS(dend_p14,'20231206_dend_p14_new.rds')

#### p15 ####
m_p15 = counts[,  rownames(meta)[meta$Patient == "Patient 15"], drop = FALSE]
hc15 = hclust(dist(t(m_p15)), method = "ward.D2")
dend_p15 = as.dendrogram(hc15)
order_p15 = m_p15[,order.dendrogram(dend_p15)] %>% colnames()
saveRDS(dend_p15,'20231206_dend_p15_new.rds')

#### p16 ####
m_p16 = counts[,  rownames(meta)[meta$Patient == "Patient 16"], drop = FALSE]
hc16 = hclust(dist(t(m_p16)), method = "ward.D2")
dend_p16 = as.dendrogram(hc16)
order_p16 = m_p16[,order.dendrogram(dend_p16)] %>% colnames()
saveRDS(dend_p16,'20231206_dend_p16_new.rds')


## Compare between patients ####
cluster_between_groups_ward = function(mat, factor) {
  
  if (!is.factor(factor)) {
    factor = factor(factor, levels = unique(factor))
  }
  
  dend_list = list()
  order_list = list()
  for(le in unique(levels(factor))) {
    m = mat[, factor == le, drop = FALSE]
    if (ncol(m) == 1) {
      order_list[[le]] = which(factor == le)
      dend_list[[le]] = structure(which(factor == le), class = "dendrogram", leaf = TRUE,
                                  height = 0, label = 1, members = 1)
    } else if(ncol(m) > 1) {
      hc1 = hclust(dist(1:ncol(m)),method = "ward.D2")
      dend_list[[le]] = reorder(as.dendrogram(hc1), wts = 1:ncol(m), agglo.FUN = mean)
      order_list[[le]] = which(factor == le)[order.dendrogram(dend_list[[le]])]
      order.dendrogram(dend_list[[le]]) = order_list[[le]]
    }
    attr(dend_list[[le]], ".class_label") = le
  }
  
  parent = as.dendrogram(hclust(dist(t(sapply(order_list, function(x) rowMeans(mat[, x, drop = FALSE]))))))
  dend_list = lapply(dend_list, function(dend) dendrapply(dend, function(node) {
    attr(node, "height") = 0
    node
  }))
  dend = merge_dendrogram(parent, dend_list)
  order.dendrogram(dend) = unlist(order_list[order.dendrogram(parent)])
  return(dend)
}

orders = ls()[ls() %>% str_detect('order_p')]
order_list = c()
for (i in orders){
  order_list= c(order_list, get(i))
}
meta = meta[order_list,]
m = counts
## meta = arrange(meta,Patient)
m = m[,rownames(meta)]
hc_patients = cluster_between_groups_ward(m , meta$Patient)
dend_test = as.dendrogram(hc_patients )
plot(dend_test )
meta[order.dendrogram(hc_patients),]
order = meta[m  %>% .[,order.dendrogram(hc_patients)] %>% colnames(),'Patient'] %>% unique
## "Patient 16" "Patient 1"  "Patient 10" "Patient 5"  "Patient 6"  "Patient 12""Patient 13" "Patient 14" "Patient 9"  "Patient 2"  "Patient 11" "Patient 15" "Patient 8"  "Patient 3"  "Patient 4" 
## "Patient 7"  "Patient 5"  "Patient 6"  "Patient 16" "Patient 1"  "Patient 10" "Patient 12" "Patient 13" "Patient 14" "Patient 9"  "Patient 2"  "Patient 11" "Patient 15" "Patient 8"  "Patient 3"  "Patient 4"    12/06

saveRDS(dend_test,'20231206_cluster_between_pateints_dend_new.rds')


order_col = c(colnames(m_p7),colnames(m_p5), colnames(m_p6), colnames(m_p16),colnames(m_p1),colnames(m_p10),colnames(m_p12),colnames(m_p13),colnames(m_p14),colnames(m_p9),colnames(m_p2),colnames(m_p11),colnames(m_p15),colnames(m_p8),colnames(m_p3),colnames(m_p4))
dend_p = as.dendrogram(hclust(dist(rbind(colMeans(t(m_p7)), rbind(colMeans(t(m_p5)),rbind(colMeans(t(m_p6)), rbind(colMeans(t(m_p16))),rbind(colMeans(t(m_p1)), rbind(colMeans(t(m_p10)), rbind(colMeans(t(m_p12)), rbind(colMeans(t(m_p13)), rbind(colMeans(t(m_p14)), rbind(colMeans(t(m_p9)), rbind(colMeans(t(m_p2)), rbind(colMeans(t(m_p11)), rbind(colMeans(t(m_p15)), rbind(colMeans(t(m_p8)),rbind(colMeans(t(m_p3)), colMeans(t(m_p4)))))))))))))))))))
plot(dend_p)
order.dendrogram(dend_p)
dend_m = merge_dendrogram(dend_p, list(dend_p7, dend_p5, dend_p6, dend_p16 ,dend_p1,dend_p10,dend_p12,dend_p13,dend_p14,dend_p9,dend_p2,dend_p11,dend_p15,dend_p8,dend_p3,dend_p4))
order.dendrogram(dend_m)
plot(dend_m)

metaa = meta[order_col,]
metaa[order.dendrogram(dend_m),] ## this is it!

saveRDS(dend_m,'20231206_dend_patients_new.rds')
saveRDS(order_col ,'20231206_order_col_new.rds')

order_col = readRDS('20231206_order_col_new.rds')
dend_m = readRDS('20231206_dend_patients_new.rds')

metaa = meta[order_col,]
row = metaa  %>% rownames()

meta = meta[order_col,]
counts_nomki = counts_nomki[,rownames(meta)]

meta_raw = meta 

col_fun = colorRamp2(c(0, 200), c(  "#F2F3F4","red3"))
colours <- list( 
  'Site' = setNames(c('gold','#D6CADD','#4B0082','#DF00FF','#A420D0'), levels(meta$Site )),
  'Histology' =setNames(c('#73C2FB','orange','red3'), c('PRAD','HGC','NEPC'))) 
meta$Patient = NULL


colAnn <- HeatmapAnnotation(df = meta,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'),show_annotation_name = FALSE)


hmap = ComplexHeatmap::Heatmap(counts_nomki, name="H-score", cluster_rows = F,cluster_columns = dend_m, row_names_gp = grid::gpar(fontsize = 15,fontface = 'bold'), column_split = 16 , na_col = "#91A3B0", col = col_fun, row_split = c(rep('1',4), rep('2',1),rep('3',1),rep('5',9), rep('6',2)),top_annotation = colAnn,row_gap = unit(c(5), "mm"), show_column_names = F, column_title = " ", row_title = " ")
draw(hmap)

saveRDS(hmap,'20240310_nomki_heat_new.rds')

counts_mki = counts_mki[row,]
counts_mki = counts_mki %>% t() %>% as.data.frame() 
colnames(counts_mki) = row
rownames(counts_mki) = 'KI67'

## "Patient 16" "Patient 1"  "Patient 10" "Patient 5"  "Patient 6"  "Patient 12""Patient 13" "Patient 14" "Patient 9"  "Patient 2"  "Patient 11" "Patient 15" "Patient 8"  "Patient 3"  "Patient 4" 
## "Patient 7"  "Patient 5"  "Patient 6"  "Patient 16" "Patient 1"  "Patient 10" "Patient 12" "Patient 13" "Patient 14" "Patient 9"  "Patient 2"  "Patient 11" "Patient 15" "Patient 8"  "Patient 3"  "Patient 4"    12/06
metaa = meta_raw

table(meta$Patient)
column_split = c(rep('1',length(which(metaa$Patient == 'Patient 7'))), rep('2',length(which(metaa$Patient == 'Patient 5'))), rep('3',length(which(metaa$Patient == 'Patient 6'))), rep('4',length(which(metaa$Patient == 'Patient 16'))), rep('5',length(which(metaa$Patient == 'Patient 1'))), rep('51',length(which(metaa$Patient == 'Patient 10'))), rep('6',length(which(metaa$Patient == 'Patient 12'))), rep('7',length(which(metaa$Patient == 'Patient 13'))), rep('8',length(which(metaa$Patient == 'Patient 14'))),rep('9',length(which(metaa$Patient == 'Patient 9'))), rep('10',length(which(metaa$Patient == 'Patient 2'))), rep('11',length(which(metaa$Patient == 'Patient 11'))), rep('12',length(which(metaa$Patient == 'Patient 15'))), rep('13',length(which(metaa$Patient == 'Patient 8'))),rep('14',length(which(metaa$Patient == 'Patient 3'))), rep('15',length(which(metaa$Patient == 'Patient 4')))  )

col_fun = colorRamp2(c(0, 100), c(  "#F2F3F4","gray4"))
hmapp = ComplexHeatmap::Heatmap(counts_mki, name="KI67", cluster_rows = F,cluster_columns = F, row_names_gp = grid::gpar(fontsize = 15,fontface = 'bold'), column_split = column_split  , na_col = "#91A3B0", col = col_fun,row_gap = unit(c(5), "mm"), show_column_names = F, column_title = " ", row_title = " ")
draw(hmapp)

saveRDS(hmapp,'20231206_mki_heat_new.rds')


## Merge heatmap ####
hmap = readRDS('20240310_nomki_heat_new.rds')
hmapp = readRDS('20231206_mki_heat_new.rds')
ht_list = hmap %v% hmapp
draw(ht_list)
saveRDS(ht_list,'20240310_heatmap_merge_new_tff.rds')

pdf('20240310_heatmap_tma_new.pdf', width = 17, height =6)
draw(ht_list)
dev.off()

# Figure 1B ####
# Spearman's ranked correlation test - break tie rank (zero-values) ####
excel = read_xlsx('~/Documents/samir_macbook/TMA/TMA97_MSKCC_v9.xlsx', sheet = 1) %>% as.data.frame()
df <- excel[,colSums(is.na(excel))<nrow(excel)]
df = distinct(df)
rownames(df) = paste0('random_',c(1:131))

meta = df %>% dplyr::select( "SoftTissue_Bone"   , "Morphology...26", "Identifier" )
rownames(meta) = rownames(df)
counts = df[,c(6:23)]; rownames(counts) = rownames(df)
counts_nomki = counts[,-12]
counts_nomki = counts_nomki[,c('ASCL1','DLL3','YAP1','INSM1','TROP2','SYP','AR','FOXA2','SOX2','EZH2')]

counts_mki = counts[,12] %>% as.data.frame(); rownames(counts_mki) = rownames(df)

counts_nomki = cbind(counts_nomki, counts_mki)

dff = cbind(counts_nomki, meta)

dff$NEPC_score_mean = rowMeans(dff[c("ASCL1",  "INSM1","FOXA2")], na.rm = TRUE)
rank(y, ties.method = "random")

#### DLL3 #####
x1 = rank(dff$DLL3, ties.method = 'random')
y1 = rank(dff$NEPC_score_mean, ties.method = 'random')
p = cor.test(x1,y1,method = 'spearman') # rho = 0.5819089, p-value < 2.2e-16 ( == 0)
p$statistic
plot(x1,y1)

plot(dff$DLL3,dff$NEPC_score_mean)

#### SYP #####
x1 = rank(dff$SYP, ties.method = 'random')
p = cor.test(x1,y1,method = 'spearman') # rho = 0.5360647, p-value =2.843e-12
p$p.value
plot(x1,y1)

#### YAP1 ####
x1 = rank(dff$YAP1, ties.method = 'random')
p = cor.test(x1,y1,method = 'spearman') # rho = -0.3353494, p-value =  9.886e-05
p$p.value
plot(x1,y1)

#### INSM1 ####
x1 = rank(dff$INSM1, ties.method = 'random')
p = cor.test(x1,y1,method = 'spearman') # rho = 0.7844606  , p-value < 2.2e-16
p$p.value
plot(x1,y1)

#### TROP2 ####
x1 = rank(dff$TROP2, ties.method = 'random')
p = cor.test(x1,y1,method = 'spearman') # rho = -0.3653019  , p-value = 2.042e-05
p$p.value
plot(x1,y1)

#### EZH2 #####
x1 = rank(dff$EZH2, ties.method = 'random')
y1 = rank(dff$., ties.method = 'random')
p = cor.test(x1,y1,method = 'spearman') # rho = 0.5819089, p-value < 2.2e-16 ( == 0)
p$statistic
plot(x1,y1)

## Dot plot - correlation ####
excel = read_xlsx('~/Documents/samir_macbook/TMA/TMA97_MSKCC_v9.xlsx', sheet = 1) %>% as.data.frame()
df <- excel[,colSums(is.na(excel))<nrow(excel)]
df = distinct(df)
rownames(df) = paste0('random_',c(1:131))

meta = df %>% dplyr::select( "SoftTissue_Bone"   , "Morphology...26", "Identifier" )
rownames(meta) = rownames(df)
# meta = df[,1:4]; rownames(meta) = rownames(df)
counts = df[,c(6:23)]; rownames(counts) = rownames(df)
counts_nomki = counts[,-12]
counts_nomki = counts_nomki[,c('ASCL1','DLL3','YAP1','INSM1','TROP2','SYP','AR')]

dff = cbind(counts_nomki, meta)

pdf(paste0('/Users/jooyoung/Google Drive/My Drive/Main_Figure_v3/Figure1B_ASCL1_DLL3_cor_Hscore_all.pdf'), width = 4.6*1.3, height = 3.9*1.3)
ggplot(dff, aes(x = ASCL1,y=DLL3)) + geom_point(size= 3.5 , aes( col = Morphology...26, fill = Morphology...26))  + scale_color_manual(values = c('#73C2FB','orange','red3'))  + theme_classic()  + guides(col=guide_legend(title="Histology")) + xlab('ASCL1 H-score')+ ylab('DLL3 H-score')+ stat_cor(label.x.npc = 0) + geom_smooth(method = "lm",size = 0.5,alpha  = 0.1,show_guide = FALSE, color = "black") + guides(fill = FALSE)+  scale_fill_manual(values = c('#73C2FB','orange','red3')) + theme(text = element_text(family = 'Helvetica'))
dev.off()

pdf(paste0('/Users/jooyoung/Google Drive/My Drive/Main_Figure_v3/Figure1B_ASCL1_YAP1_cor_Hscore_all.pdf'), width = 4.6*1.3, height = 3.9*1.3)
ggplot(dff, aes(x = ASCL1,y=YAP1)) + geom_point(size= 3.5 , aes( col = Morphology...26, fill = Morphology...26))  + scale_color_manual(values = c('#73C2FB','orange','red3'))  + theme_classic()  + guides(col=guide_legend(title="Histology")) + xlab('ASCL1 H-score')+ ylab('YAP1 H-score')+ stat_cor(label.x.npc = 0) + geom_smooth(method = "lm",size = 0.5,alpha  = 0.1,show_guide = FALSE, color = "black") + guides(fill = FALSE)+  scale_fill_manual(values = c('#73C2FB','orange','red3')) + theme(text = element_text(family = 'Helvetica'))
dev.off()


pdf(paste0('/Users/jooyoung/Google Drive/My Drive/Main_Figure_v3/Figure1B_ASCL1_SYP_cor_Hscore_all.pdf'), width = 4.6*1.3, height = 3.9*1.3)
ggplot(dff, aes(x = ASCL1,y=SYP)) + geom_point(size= 3.5 , aes( col = Morphology...26, fill = Morphology...26))  + scale_color_manual(values = c('#73C2FB','orange','red3'))  + theme_classic()  + guides(col=guide_legend(title="Histology")) + xlab('ASCL1 H-score')+ ylab('SYP H-score')+ stat_cor(label.x.npc = 0) + geom_smooth(method = "lm",size = 0.5,alpha  = 0.1,show_guide = FALSE, color = "black") + guides(fill = FALSE)+  scale_fill_manual(values = c('#73C2FB','orange','red3')) + theme(text = element_text(family = 'Helvetica'))
dev.off()

pdf(paste0('/Users/jooyoung/Google Drive/My Drive/PNAS_fin/Main/Figure1_EZH2_KI67_cor_Hscore_all.pdf'), width = 4.6*1.3, height = 3.9*1.3)
ggplot(df, aes(x = EZH2,y=Ki67)) + geom_point(size= 3.5 , aes( col =Morphology, fill = Morphology))  + scale_color_manual(values = c('#73C2FB','orange','red3'))  + theme_classic()  + guides(col=guide_legend(title="Histology")) + xlab('EZH2 H-score')+ ylab('Ki67 H-score')+ stat_cor(label.x.npc = 0) + geom_smooth(method = "lm",size = 0.5,alpha  = 0.1,show_guide = FALSE, color = "black") + guides(fill = FALSE)+  scale_fill_manual(values = c('#73C2FB','orange','red3')) + theme(text = element_text(family = 'Helvetica'))
dev.off()


# Supplemetary Fig 1 with zero-inflated modified wilcox test ####
library(rstatix)
library(ZIR)
excel = read_xlsx('~/Documents/samir_macbook/TMA/TMA97_MSKCC_v9.xlsx', sheet = 1) %>% as.data.frame()
df <- excel[,colSums(is.na(excel))<nrow(excel)]
df = distinct(df)
rownames(df) = paste0('random_',c(1:131))

colnames(df)[26] = 'Morphology'
colnames(df)[21] = 'TTF3'
colnames(df)[17] = 'Ki67'

for (y in colnames(df)[6:23] ){
  if ( !  y %in% c('TP63','NEUROD1','POU2F3')){
    formula = as.formula( paste(y, 'Morphology', sep="~") )
    
    tmp = df[,c(y,'Morphology')]
    
    ADC = tmp[tmp$Morphology == "ADC",y][!tmp[tmp$Morphology == "ADC",y] %>% is.na()]
    HGC = tmp[tmp$Morphology == "HGC",y][!tmp[tmp$Morphology == "HGC",y] %>% is.na()]
    NEPC = tmp[tmp$Morphology == "NEPC",y][!tmp[tmp$Morphology == "NEPC",y] %>% is.na()]
    
    stat.test = data.frame(matrix(ncol = 13, nrow = 3)); colnames(stat.test) = c(".y." ,"group1","group2","n1","n2","statistic","p","p.adj","p.adj.signif","y.position","groups","xmin","xmax")
    
    stat.test[,1] = y
    stat.test$group1 = c('ADC','ADC','HGC')
    stat.test$group2 = c('HGC','NEPC','NEPC')
    
    stat.test$n1[1] = length(ADC)
    stat.test$n1[2] = length(ADC)
    stat.test$n1[3] = length(HGC)
    stat.test$n2[1] = length(HGC)
    stat.test$n2[2] = stat.test$n2[3] = length(NEPC)
    
    stat.test$statistic[1] = ziw(ADC, HGC, perm = T)$statistics
    stat.test$statistic[2] = ziw(ADC, NEPC, perm = T)$statistics
    stat.test$statistic[3] = ziw(HGC, NEPC, perm = T)$statistics
    
    stat.test$p[1] = ziw(ADC, HGC, perm = T)$p.value
    stat.test$p[2] = ziw(ADC, NEPC, perm = T)$p.value
    stat.test$p[3] = ziw(HGC, NEPC, perm = T)$p.value
    
    stat.test$p.adj[1] = p.adjust(ziw(ADC, HGC, perm = T)$p.value, method = "bonferroni")
    stat.test$p.adj[2] = p.adjust(ziw(ADC, NEPC, perm = T)$p.value, method = "bonferroni")
    stat.test$p.adj[3] = p.adjust(ziw(HGC, NEPC, perm = T)$p.value, method = "bonferroni")
    
    stat.test = stat.test %>% add_significance() 
    maxval = df[,y][! df[,y] %>% is.na()] %>% max()
    
    stat.test$manual_position <- maxval * c(1.05,1.2,1.11)
    stat.test$label <- stat.test$p.adj.signif
    
    ylab = y
    if (y == 'Ki67'){
      ylab = 'Ki-67'
    }
    
    p = ggplot(df, aes(x= Morphology, y = df[,y], col =Morphology )) + geom_boxplot()   + scale_color_manual(values = c('#73C2FB','orange','red3'))  + theme_classic()  + guides(col=guide_legend(title="Histology")) + theme(axis.title.x = element_blank()) +  ggsignif::geom_signif(data=stat.test,aes(xmin=group1,xmax=group2,annotations=label,y_position=manual_position),manual=TRUE, inherit.aes=FALSE) + ylab(paste0(y, ' H-Score'))
    
    pdf(paste0('/Users/jooyoung/Google Drive/My Drive/PNAS_fin/20240131_rev_v1/Supp_Figure1_ziw_',y,'.pdf'), width = 4, height = 4)
    print(p)
    dev.off()
  }
}


