library(data.table)
library(dplyr)
library(ggplot2)
library(smoother)
library(splines)
library(R.utils)
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

annot=fread("metadata/fig6_samples.csv")
all = list()
for (i in 1:nrow(annot)){ 
  row <- annot[i,]
  print(row)
  tmp = fread(row$fil, stringsAsFactors = FALSE)
  tmp$ss = row$ss
  tmp$prot = row$prot
  tmp$rep = row$rep
  tmp$sub = row$sub
  tmp$condition = row$condition
  tmp$complex = row$complex
  print(tmp)
  all[[i]] = tmp
}
all = do.call(rbind, all)

# correct coordinates - make everything relative to 3'SS
allc <- all %>% mutate(V2 = ifelse(V1=="ACT1_dG", V2-450,
                                   ifelse(V1=="WT_actin", V2-442, V2)))
# remove the 5 branch point nucleotides for C star complex
alld = allc %>% filter(V2!="-46" & V2!="-45" & V2!="-44" & V2!="-43" & V2!="-42")
                             
ff=alld

# to merge replicates
ff <- alld %>% ungroup() %>% group_by(V1,V2,ss,prot,sub,condition,complex) %>% 
  summarise(V5=sum(V5)) %>% ungroup() %>% group_by(ss,prot,sub,condition, complex) %>% mutate(lib_size=sum(V5), norm=V5/lib_size)

# to show crosslinking distributions
reps_merg = ff
count_reps_merg = reps_merg %>% ungroup() %>% 
  group_by(V1,sub,condition,ss,prot, complex) %>% summarise(total=sum(V5)) %>% ungroup() %>% 
  mutate(V1=paste0(V1, " snRNA") %>% gsub(".*dG snRNA|WT_actin snRNA","pre-mRNA substrate", .), 
         prot=paste0(prot, " Prp22"))

ggplot(count_reps_merg %>% filter(condition=="UV"), aes(x=interaction(sub,ss,prot,condition, complex), y=total, fill=V1)) + geom_col() + 
  theme_few() + ylab("number of crosslinks") + coord_flip() + facet_wrap(~complex+sub, ncol=2, scales="free") +
  ggtitle("all reps together") + xlab("") + theme(legend.title = element_blank())
ggsave("prp22_all_reps_together_-snRNA-substrate-proportions.eps", height=5, width=12)

# to show snRNA enrichment
reps_merg_sn = ff 
count_reps_merg_sn = reps_merg_sn %>% ungroup() %>%
  mutate(V1=paste0(V1, " snRNA") %>% gsub(".*dG snRNA|.*actin snRNA|.*MS2 snRNA","pre-mRNA substrate", .), 
         prot=paste0(prot, " Prp22")) %>% filter(V1 != "pre-mRNA substrate") %>% group_by(V1,V2,ss,prot,sub,condition,complex) %>% 
  summarise(V5=sum(V5)) %>% ungroup() %>% group_by(ss,prot,sub,condition,complex) %>% mutate(lib_size=sum(V5)) %>% ungroup() %>%
  group_by(V1,sub,condition,ss,prot,complex) %>% summarise(total=sum(V5)) %>% ungroup()
# get a lib size and lib size norm
snrna_enrichment = count_reps_merg_sn %>% ungroup() %>% dplyr::group_by(sub, condition, ss, prot, complex) %>% 
  dplyr::mutate(libsize=sum(total)) %>% ungroup() %>%
  dplyr::group_by(sub, ss, prot, complex) %>% 
  dplyr::mutate(log2fc= log2(total/libsize*1000000) - log2(total[condition=="ctrl"]/libsize[condition=="ctrl"]*1000000)) %>%
  filter(condition=="UV")
ggplot(snrna_enrichment, aes(x=fct_relevel(as.factor(V1), "U1 snRNA", "U4 snRNA","U6 snRNA","U2 snRNA","U5 snRNA"), y=log2fc, color=interaction(sub,ss,prot,condition,complex), fill=V1)) + 
  geom_bar(stat = "identity",position = "dodge") +
  xlab("") + theme_few() + theme(legend.title = element_blank(), panel.grid.major.y = element_line(color="grey"),
                                 panel.grid.minor.y = element_line(color="grey")) +  ylab("log2FC UV/ctrl") + coord_flip() +
  facet_wrap(~sub, ncol=2, scales="free_y") + guides(group=guide_legend()) 

ggsave("prp22_all_reps_together_-snRNA-enrichment.eps", height=5, width=8)


# to ctrl-normalise the data
ff <- ff %>% dplyr::group_by(sub,ss,prot,condition, complex) %>% dplyr::mutate(smoothed= smth.gaussian(norm, window=10))
ff <- data.frame(ff)
ff$smoothed[ff$smoothed<0] <- 0
ff$smoothed[is.na(ff$smoothed)] <- 0
# subtract ctrlcDNAs from UV cDNAs
ff <- ff %>% dplyr::group_by(sub,ss,prot,V2,complex) %>% dplyr::mutate(ctrl_minus=ifelse(smoothed>0,smoothed-smoothed[condition=="ctrl"],0))
ff <- data.frame(ff)
ff$ctrl_minus[ff$ctrl_minus<0] <- 0
ff$ctrl_minus[is.na(ff$ctrl_minus)] <- 0
ff <- ff %>% dplyr::group_by(sub,ss,prot,condition,complex) %>% dplyr::mutate(ctrl_minus= smth.gaussian(ctrl_minus, window=10))
ff <- data.frame(ff)
ff$ctrl_minus[ff$ctrl_minus<0] <- 0
ff$ctrl_minus[is.na(ff$ctrl_minus)] <- 0

ff_uv <- ff %>% filter(condition=="UV") %>% filter(grepl("dG|actin", V1)) %>%
  filter(V2 > -100 & V2 < 100)

ggplot(ff_uv, aes(x=V2, y=ctrl_minus, group=interaction(ss,sub,prot,complex),color=interaction(prot,complex))) +
  geom_vline(xintercept = c(0, -44), linetype=5) +
  scale_x_continuous(limits =c(-100,100), breaks=seq(-100,100,5), labels=c(insert(seq(-100,100,10), ats=2:length(seq(-100,100,10)), ""))) +
  scale_color_manual(values = cbp1) +
  geom_line(aes(y=ctrl_minus),size=1) +
  ylab("normalised score") +
  xlab("") +
  theme_few() +
  facet_wrap(~sub+ss, scales="free", ncol=1) +
  theme(text = element_text(size=20))
ggsave("PRP22_prp18ksprofile_repCombined_ctrlminus.eps", height=8, width=13)

ggplot(ff_uv, aes(x=V2, y=smoothed, group=interaction(ss,sub,prot,complex),color=interaction(prot,complex))) +
  geom_vline(xintercept = c(0, -44), linetype=5) +
  scale_x_continuous(limits =c(-100,100), breaks=seq(-100,100,5), labels=c(insert(seq(-100,100,10), ats=2:length(seq(-100,100,10)), ""))) +
  scale_color_manual(values = cbp1) +
  geom_line(aes(y=smoothed),size=1) +
  ylab("normalised score") +
  xlab("") +
  theme_few() +
  facet_wrap(~sub+ss, scales="free", ncol=1) +
  theme(text = element_text(size=20))
ggsave("PRP22_prp18ksprofile_repCombined_norm.eps", height=8, width=13)

ggplot(ff_uv, aes(x=V2, y=smth.gaussian(V5, window=10), group=interaction(ss,sub,prot,complex),color=interaction(prot,complex))) +
  geom_vline(xintercept = c(0, -44), linetype=5) +
  scale_x_continuous(limits =c(-100,100), breaks=seq(-100,100,5), labels=c(insert(seq(-100,100,10), ats=2:length(seq(-100,100,10)), ""))) +
  scale_color_manual(values = cbp1) +
  geom_line(aes(y=smth.gaussian(V5, window=10)),size=1) +
  ylab("normalised score") +
  xlab("") +
  theme_few() +
  facet_wrap(~ss, scales="free", ncol=1) +
  theme(text = element_text(size=20))
ggsave("PRP22_prp18ksprofile_repCombined_justSmoothNotnorm.eps", height=8, width=13)








#######################
ff <- ff %>% group_by(ss,prot,rep,genotype) %>% mutate(norm=V5/sum(V5)) %>% filter(V2 != -44 & V2 != -43 & V2 != -45 & V2 != -46 & V2 != -47)

ff_rep1 <- ff %>% filter(rep==1)
ff_rep2 <- ff %>% filter(rep==2)
ggplot(ff_rep2, aes(x=V2, y=V5, color=interaction(prot,rep,genotype))) +
  geom_vline(xintercept = c(0, -44), linetype=5) +
  scale_color_manual(values = cbp1) +
  #geom_col(fill="grey") +
  geom_line(aes(y=smth.gaussian(V5, window=20)), size=1) +
  #geom_smooth(aes(y=V5), method = "loess", span = 0.2, se=FALSE) #+
  #geom_smooth(aes(y=norm, color=interaction(ss,prot)), se=F,method="glm", formula=y~ns(x,10)) +
  scale_x_continuous(limits =c(-100,100), breaks=seq(-100,100,5), labels=c(insert(seq(-100,100,10), ats=2:length(seq(-100,100,10)), ""))) +
  ylab("cDNAs") +
  xlab("") +
  facet_wrap(~genotype +ss) +
  ylim(0,35) +
  theme_few()

ggsave("fig7_prp18kd_prp22_REP2.eps", height=3, width=15)