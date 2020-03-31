library(data.table)
library(dplyr)
library(ggplot2)
library(smoother)
library(R.utils)
library(zoo)
library(ggthemes)
library(forcats)
library(cowplot)
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

annot=fread("metadata/fig4_samples.csv")

# get genomic xlsites
fls = basename(annot[[1]]) %>% gsub(".small.Aligned.out.sorted.cdnacounts",".Aligned.out.sorted.cdnacounts",.)
fil = paste0("../results/genome_xlsites/", fls)
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
act1EXT <- all %>% filter(grepl("ACT1",V1)) %>% mutate(V2=as.numeric(V2)-403)
act1 <- all %>% filter(V1=="WT_actin") %>% mutate(V2=as.numeric(V2)-398)

# merge everything together
ff = rbind(act1EXT, act1)


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

ff = ff %>% filter(V2 > -100 & V2 < 150)

# to normalise by library size
#ff <- ff %>% dplyr::group_by(sub,ss,prot,condition) %>% dplyr::mutate(norm=as.numeric(V5)/lib_size) 
ff <- ff %>% dplyr::group_by(sub,ss,prot,condition) %>% dplyr::mutate(smoothed= smth.gaussian(norm, window=20))
ff <- data.frame(ff)
ff$smoothed[ff$smoothed<0] <- 0
ff$smoothed[is.na(ff$smoothed)] <- 0
# subtract ctrlcDNAs from UV cDNAs
ff <- ff %>% dplyr::group_by(sub,ss,prot,V2) %>% dplyr::mutate(ctrl_minus=ifelse(norm>0,norm-norm[condition=="ctrl"],0))
ff <- data.frame(ff)
ff$ctrl_minus[ff$ctrl_minus<0] <- 0
ff$ctrl_minus[is.na(ff$ctrl_minus)] <- 0
ff <- ff %>% dplyr::group_by(sub,ss,prot,condition) %>% dplyr::mutate(ctrl_minus= smth.gaussian(ctrl_minus, window=20))
ff <- data.frame(ff)
ff$ctrl_minus[ff$ctrl_minus<0] <- 0
ff$ctrl_minus[is.na(ff$ctrl_minus)] <- 0

ff$condition = factor(ff$condition, levels=c("UV", "ctrl"))

# Replicate, lines all on one graph
ff_uv <- ff %>% filter(condition=="UV")

ggplot(ff_uv, aes(x=V2, y=ctrl_minus, color=interaction(ss,prot))) +
  geom_vline(xintercept = c(0), linetype=5) +
  scale_color_manual(values = cbp1) +
  geom_line(aes(y=ctrl_minus),size=1) +
  scale_x_continuous(limits =c(-50,100), breaks=seq(-50,100,5), labels=c(insert(seq(-50,100,10), ats=2:length(seq(-50,100,10)), ""))) +
  ylab("normalised score") +
  xlab("") +
  theme_few() +
  facet_wrap(~ss, scales="free_y",ncol=1) +
  theme(text = element_text(size=20))
ggsave("fig4_libsize.eps", height=6, width=20)

ff$ss = factor(ff$ss, levels=c("WT","ext20", "ext40","extSL"))
ggplot(ff, aes(x=V2, y=smoothed, color=sub, linetype=condition)) +
  geom_vline(xintercept = c(0,44,64,84), linetype=5) +
  scale_color_manual(values = c("#6383be","#b8babc")) +
  geom_line(aes(y=smoothed),size=1) +
  scale_x_continuous(limits =c(-50,110), breaks=seq(-50,110,10), labels=c(insert(seq(-50,110,20), ats=2:length(seq(-50,110,20)), ""))) +
  ylab("normalised score") +
  xlab("") +
  theme_few() +
  facet_wrap(~ss, ncol=1,scales="free") +
  theme(text = element_text(size=20))
ggsave("fig4_yeast_genome_seperate_reps_combined.eps", height=8, width=8)

#raw zoom
ggplot(ff %>% filter(ss=="ext40" | ss=="extSL"), aes(x=V2, y=V5, fill=condition)) + geom_col() + theme_few() + xlim(-10,20) +
  facet_wrap(~ss, ncol=1,scales="free") + ylab("# crosslinks") + xlab("distance from brA") +
  scale_fill_manual(values = c("#6383be","#cccccc"))
ggsave("fig4_yeast_genome_stem_loop_zoom.eps", height=4, width=6)

#raw
ggplot(ff, aes(x=V2, y=smoothed, fill=condition)) +
  geom_vline(xintercept = c(0), linetype=5) +
  scale_fill_manual(values = cbp1) +
  geom_col() +
  scale_x_continuous(limits =c(-50,100), breaks=seq(-50,100,5), labels=c(insert(seq(-50,100,10), ats=2:length(seq(-50,100,10)), ""))) +
  ylab("normalised score") +
  xlab("") +
  theme_few() +
  facet_wrap(~ss+condition, scales="free_y",ncol=1) +
  theme(text = element_text(size=20))
ggsave("fig4_raw_repscombined_libsize.eps", height=12, width=6)


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

#######################################################################
#################### Supplementary, replicates separate ####################
# merge everything together
ff = rbind(act1EXT, act1)

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

ff = ff %>% filter(V2 > -100 & V2 < 150)

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
ff$ss = factor(ff$ss, levels=c("WT","ext20", "ext40","extSL"))

ggplot(ff, aes(x=V2, y=smoothed, colour=rep)) +
  geom_vline(xintercept = c(0,44,64,84), linetype=5) +
  scale_color_manual(values = c("#6383be","#23344f","#2c518c")) +
  geom_line(aes(y=smoothed,linetype=condition),size=1) +
  scale_x_continuous(limits =c(-50,110), breaks=seq(-50,110,10), labels=c(insert(seq(-50,110,20), ats=2:length(seq(-50,110,20)), ""))) +
  ylab("normalised score") +
  xlab("") +
  theme_few() +
  facet_wrap(~ss+rep, scales="free",ncol=2) +
  theme(text = element_text(size=18))  + theme(legend.position = "none")
ggsave("fig4_yeast_genome_seperate.eps", height=7, width=8)

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
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                