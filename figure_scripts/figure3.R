library(data.table)
library(dplyr)
library(ggplot2)
library(smoother)
library(splines)
library(R.utils)
library(zoo)
library(ggthemes)
library(forcats)
library(cowplot)
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

annot=fread("metadata/fig3_samples.csv")

# get genomic xlsites
fls = basename(annot[[1]]) %>% gsub(".small.Aligned.out.sorted.cdnacounts",".Aligned.out.sorted.cdnacounts",.)
fil = paste0("../pre_processing/results/genome_xlsites/", fls)
fil

all = list()
for (i in 1:length(fil)){ 
  row <- fil[[i]]
  nm = basename(row) %>% gsub(".Aligned.*","",.)
  tmp = fread(row, stringsAsFactors = FALSE)
  tmp$sample = nm
  all[[i]] = tmp
}
all = do.call(rbind, all)

gxl = all %>% dplyr::group_by(sample) %>% dplyr::summarise(total_Genome_xl = sum(V5))

# merge with annotations
annot$sample = basename(annot$fil) %>% gsub(".small.*","",.)
annot2 = left_join(annot,gxl,by=c("sample"))


all = list()
for (i in 1:nrow(annot2)){ 
  row <- annot2[i,]
  print(row)
  tmp = fread(row$fil, stringsAsFactors = FALSE)
  tmp$ss = row$ss
  tmp$prot = row$prot
  tmp$rep = row$rep
  tmp$sub = row$sub
  tmp$condition = row$condition
  tmp$genome_xls = row$total_Genome_xl
  substrate = tmp$V1[!grepl("^U", tmp$V1)][[1]]
  du_df = data.frame(V1=substrate,V2=seq(-1,1200), V3=seq(-1,1200), V4=".",V5=0,V6="+",ss=row$ss, prot=row$prot,rep=row$rep, sub=row$sub, condition=row$condition,genome_xls=row$total_Genome_xl)
  ne = full_join(tmp,du_df,by=c("V1","V2","V3", "V4","V5","V6","ss", "prot","rep", "sub", "condition", "genome_xls")) %>% arrange(-V2)
  ne = ne %>% ungroup() %>% dplyr::group_by(V1,V2,ss,prot,rep,sub,condition,genome_xls) %>% filter(V5==max(V5))
  all[[i]] = ne
}

all = do.call(rbind, all)
all <- all %>% dplyr::group_by(sub,condition,ss,prot) %>% dplyr::mutate(lib_size=sum(as.numeric(V5)))
ubc4 <- all %>% filter(V1=="WT_long"| V1=="AC_long") %>% mutate(V2=as.numeric(V2)-232)
act1 <- all %>% filter(V1=="WT_actin"| V1=="AC_actin") %>% mutate(V2=as.numeric(V2)-398)

# merge everything together
ff = rbind(ubc4, act1)

# for debranch figure s2
fdb = ff %>% ungroup() %>% filter(V1=="WT_long" & condition=="UV" & rep != "rep3") %>% mutate(norm=V5/lib_size)
ggplot(fdb, aes(x=V2+1,y=norm)) + geom_col() + xlim(-10,50) + theme_few() + facet_wrap(~rep, ncol=1, scales="free_y")
ggsave("DBR1.pdf",width=6, height=4)

# remove the 5 branch point nucleotides
# remember the signal is at -1, despite the fact the actual branch point is at 0
ff = ff %>% filter(V2!=-3 & V2!=-2 & V2!=-1 & V2!=0 & V2!=1)

#################### MAIN FIGURE REPLICATES MERGED ####################

# to merge replicates
ff = ff %>% ungroup() %>% 
  dplyr::group_by(ss,prot,sub,condition) %>%
  dplyr::mutate(lib_size=sum(V5)) %>%
  filter(!grepl("^U", V1)) %>%
  ungroup() %>% dplyr::group_by(V1,V2,ss,prot,sub,condition,lib_size) %>%
  dplyr::summarise(V5=sum(V5), genome_xls=sum(genome_xls)) %>% ungroup()

# to normalise by genome xls
ff$norm = ff$V5/ff$genome_xls
ff[is.na(ff)] <- 0

ff = ff %>% filter(V2 > -100 & V2 < 100)

# to normalise by library size
#ff <- ff %>% dplyr::group_by(sub,ss,prot,condition) %>% dplyr::mutate(norm=as.numeric(V5)/lib_size) 
ff <- ff %>% dplyr::group_by(sub,ss,prot,condition) %>% dplyr::mutate(smoothed= smth.gaussian(norm, window=20))
ff <- data.frame(ff)
ff$smoothed[ff$smoothed<0] <- 0
ff$smoothed[is.na(ff$smoothed)] <- 0
# subtract ctrlcDNAs from UV cDNAs
ff <- ff %>% dplyr::group_by(sub,ss,prot,V2) %>% dplyr::mutate(ctrl_minus=ifelse(smoothed>0,smoothed-smoothed[condition=="ctrl"],0))
ff <- data.frame(ff)
ff$ctrl_minus[ff$ctrl_minus<0] <- 0
ff$ctrl_minus[is.na(ff$ctrl_minus)] <- 0
ff <- ff %>% dplyr::group_by(sub,ss,prot,condition) %>% dplyr::mutate(ctrl_minus= smth.gaussian(ctrl_minus, window=20))
ff <- data.frame(ff)
ff$ctrl_minus[ff$ctrl_minus<0] <- 0
ff$ctrl_minus[is.na(ff$ctrl_minus)] <- 0


# Replicate, lines all on one graph
ff_uv <- ff %>% filter(condition=="UV")

ggplot(ff_uv, aes(x=V2, y=ctrl_minus, color=interaction(ss,prot))) +
  geom_vline(xintercept = c(0), linetype=5) +
  scale_color_manual(values = cbp1) +
  geom_line(aes(y=ctrl_minus),size=1) +
  scale_x_continuous(limits =c(-50,60), breaks=seq(-50,60,5), labels=c(insert(seq(-50,60,10), ats=2:length(seq(-50,60,10)), ""))) +
  ylab("normalised score") +
  xlab("") +
  theme_few() +
  facet_wrap(~sub, scales="free_y") +
  theme(text = element_text(size=20))
ggsave("fig3_prp16_part_a_repCombined_ctrlminus.eps", height=6, width=20)

# snRNA mapping
reps_merg = all
count_reps_merg = reps_merg %>% ungroup() %>% group_by(V1,V2,ss,prot,sub,condition) %>% 
  summarise(V5=sum(V5)) %>% ungroup() %>% group_by(ss,prot,sub,condition) %>% mutate(lib_size=sum(V5)) %>% ungroup() %>%
group_by(V1,sub,condition,ss,prot) %>% summarise(total=sum(V5)) %>% ungroup() %>% 
  mutate(V1=paste0(V1, " snRNA") %>% gsub(".*long snRNA|.*actin snRNA","pre-mRNA substrate", .), 
         ss=gsub("WT","AG", ss),
         prot=paste0(prot, " Prp16"))
# get a lib size and lib size norm
count_reps_merg %>% ungroup() %>% dplyr::group_by(sub, condition, ss, prot) %>% 
  dplyr::mutate(libsize=sum(total)) %>% ungroup() %>%
  dplyr::group_by(sub, ss, prot) %>% 
  dplyr::mutate(log2fc= log2(total/libsize*1000000) - log2(total[condition=="ctrl"]/libsize[condition=="ctrl"]*1000000))


ggplot(count_reps_merg, aes(x=interaction(sub,ss,prot,condition), y=total, fill=V1)) + geom_col() + 
  theme_few() + ylab("number of crosslinks") + coord_flip() + facet_wrap(~condition+sub, ncol=1, scales="free") +
  ggtitle("all reps together") + xlab("") + theme(legend.title = element_blank())
ggsave("fig2_all_reps_together_UBC4-snRNA-substrate-proportions.eps", height=6, width=5)

# JUST snRNA proportion enrichment mapping
reps_merg_sn = all
count_reps_merg_sn = reps_merg_sn %>% ungroup() %>%
  mutate(V1=paste0(V1, " snRNA") %>% gsub(".*long snRNA|.*actin snRNA","pre-mRNA substrate", .), 
         ss=gsub("WT","AG", ss),
         prot=paste0(prot, " Prp16")) %>% filter(V1 != "pre-mRNA substrate") %>% group_by(V1,V2,ss,prot,sub,condition) %>% 
    summarise(V5=sum(V5)) %>% ungroup() %>% group_by(ss,prot,sub,condition) %>% mutate(lib_size=sum(V5)) %>% ungroup() %>%
    group_by(V1,sub,condition,ss,prot) %>% summarise(total=sum(V5)) %>% ungroup()
# get a lib size and lib size norm
snrna_enrichment = count_reps_merg_sn %>% ungroup() %>% dplyr::group_by(sub, condition, ss, prot) %>% 
  dplyr::mutate(libsize=sum(total)) %>% ungroup() %>%
  dplyr::group_by(sub, ss, prot) %>% 
  dplyr::mutate(log2fc= log2(total/libsize*1000000) - log2(total[condition=="ctrl"]/libsize[condition=="ctrl"]*1000000)) %>%
  filter(condition=="UV")
ggplot(snrna_enrichment, aes(x=fct_relevel(as.factor(V1), "U1 snRNA", "U4 snRNA","U6 snRNA","U2 snRNA","U5 snRNA"), y=log2fc, color=interaction(sub,ss,prot,condition), fill=V1)) + 
  geom_bar(stat = "identity",position = "dodge") +
 xlab("") + theme_few() + theme(legend.title = element_blank(), panel.grid.major.y = element_line(color="grey"),
                                panel.grid.minor.y = element_line(color="grey")) +  ylab("log2FC UV/ctrl") + coord_flip() +
  facet_wrap(~sub, ncol=2, scales="free_y") + guides(group=guide_legend())

ggsave("fig2_all_reps_together_UBC4-snRNA-enrichment.eps", height=5, width=8)

#######################################################################
#################### Supplementary, replicates separate ####################
# merge everything together
ff = rbind(ubc4, act1)

# remove the 5 branch point nucleotides
# remember the signal is at -1, despite the fact the actual branch point is at 0
ff = ff %>% filter(V2!=-3 & V2!=-2 & V2!=-1 & V2!=0 & V2!=1)

# to get lib size
ff = ff %>% ungroup() %>% 
  dplyr::group_by(ss,prot,sub,condition, rep) %>%
  dplyr::mutate(lib_size=sum(V5)) %>%
  filter(!grepl("^U", V1)) %>%
  ungroup() 

# to normalise by genome xls
ff$norm = ff$V5/ff$genome_xls
ff[is.na(ff)] <- 0

ff = ff %>% filter(V2 > -100 & V2 < 100)

# to normalise by library size
# ff <- ff %>% dplyr::group_by(sub,ss,prot,condition,rep) %>% dplyr::mutate(norm=as.numeric(V5)/lib_size) 
ff <- ff %>% dplyr::group_by(sub,ss,prot,condition,rep) %>% dplyr::mutate(smoothed= smth.gaussian(norm, window=20))
ff <- data.frame(ff)
ff$smoothed[ff$smoothed<0] <- 0
ff$smoothed[is.na(ff$smoothed)] <- 0
# subtract ctrlcDNAs from UV cDNAs
ff <- ff %>% dplyr::group_by(sub,ss,prot,V2,rep) %>% dplyr::mutate(ctrl_minus=ifelse(smoothed>0,smoothed-smoothed[condition=="ctrl"],0))
ff <- data.frame(ff)
ff$ctrl_minus[ff$ctrl_minus<0] <- 0
ff$ctrl_minus[is.na(ff$ctrl_minus)] <- 0
ff <- ff %>% dplyr::group_by(sub,ss,prot,condition,rep) %>% dplyr::mutate(ctrl_minus= smth.gaussian(ctrl_minus, window=20))
ff <- data.frame(ff)
ff$ctrl_minus[ff$ctrl_minus<0] <- 0
ff$ctrl_minus[is.na(ff$ctrl_minus)] <- 0

ff$condition = factor(ff$condition, levels=c("UV","ctrl"))
ff=ff %>% mutate(ss=gsub("WT","AG", ss),prot=paste0(prot, " Prp16"))
aff = ff %>% filter(sub=="ACT1")
aff_pl = ggplot(aff, aes(x=V2, y=smoothed, colour=rep)) +
  geom_vline(xintercept = c(0), linetype=5) +
  scale_color_manual(values = c("#6383be","#23344f","#2c518c")) +
  geom_line(aes(y=smoothed,linetype=condition),size=1) +
  scale_x_continuous(limits =c(-50,50), breaks=seq(-50,50,10), labels=c(insert(seq(-50,60,20), ats=2:length(seq(-50,60,20)), ""))) +
  ylab("normalised score") +
  xlab("") +
  theme_few() +
  facet_wrap(~ss+prot, scales="free_y",nrow=1) +
  theme(text = element_text(size=18))  + theme(legend.position = "none")

uff = ff %>% filter(sub=="UBC4")
uff_pl = ggplot(uff, aes(x=V2, y=smoothed, colour=rep)) +
  geom_vline(xintercept = c(0), linetype=5) +
  scale_color_manual(values = c("#6383be","#23344f","#2c518c")) +
  geom_line(aes(y=smoothed,linetype=condition),size=1) +
  scale_x_continuous(limits =c(-50,50), breaks=seq(-50,50,10), labels=c(insert(seq(-50,60,20), ats=2:length(seq(-50,60,20)), ""))) +
  ylab("normalised score") +
  xlab("") +
  theme_few() +
  facet_wrap(~ss+prot, scales="free_y",nrow=1) +
  theme(text = element_text(size=18))  + theme(legend.position = "none")

# crosslink graphs
reps_merg = all
count_reps_merg = reps_merg %>% ungroup() %>% group_by(ss,prot,sub,condition,rep) %>% mutate(lib_size=sum(V5)) %>% ungroup() %>%
  group_by(V1,sub,condition,ss,prot,rep) %>% summarise(total=sum(V5)) %>% ungroup() %>% 
  mutate(V1=paste0(V1, " snRNA") %>% gsub(".*long snRNA|.*actin snRNA","pre-mRNA substrate", .), 
         ss=gsub("WT","AG", ss),
         prot=paste0(prot, " Prp16"))
acrm_uv = count_reps_merg %>% filter(condition=="UV") %>% filter(sub=="ACT1") %>% filter(!grepl("^U", V1)) 
acrm_pl = ggplot(acrm_uv, group=rep,aes(x=rep, y=total,fill=rep)) + geom_col() + scale_fill_manual(values = c("#6383be","#23344f","#2c518c")) +
  theme_few() + ylab("# substrate crosslinks") + facet_wrap(~ss+prot, nrow=1,scales="free_x") +
  xlab("") + theme(legend.title = element_blank()) + theme(legend.position = "none")+
  theme(text = element_text(size=18))

ucrm_uv = count_reps_merg %>% filter(condition=="UV") %>% filter(sub=="UBC4") %>% filter(!grepl("^U", V1)) 
ucrm_pl = ggplot(ucrm_uv, group=rep,aes(x=rep, y=total,fill=rep)) + geom_col() + scale_fill_manual(values = c("#6383be","#23344f","#2c518c")) +
  theme_few() + ylab("# substrate crosslinks") + facet_wrap(~ss+prot, nrow=1,scales="free_x") +
  xlab("") + theme(legend.title = element_blank()) + theme(legend.position = "none") +
  theme(text = element_text(size=18))

# put it all together
pdf("genome_norm_fin_s3.pdf",height=12,width=10)
plot_grid(acrm_pl, aff_pl, ucrm_pl, uff_pl, labels = c('a', 'b','c','d'), label_size = 12,ncol=1,align="v", rel_heights = c(1, 2, 1, 2))
dev.off()
############################################################################
