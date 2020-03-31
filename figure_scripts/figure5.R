library(data.table)
library(dplyr)
library(ggplot2)
library(smoother)
library(R.utils)
library(ggthemes)

cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

annot=fread("metadata.csv")

# get genomic xlsites
fls = basename(annot[[1]]) %>% gsub(".small.Aligned.out.sorted.cdnacounts|.Aligned.out.forward.sorted.bam.intronRemoved.cDNAcounts",".Aligned.out.sorted.cdnacounts",.)
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
annot$sample = basename(annot$fil) %>% gsub(".small.*|.Aligned.*intronRemoved.*","",.)
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
  tmp$complex = row$complex
  tmp$genome_xls = row$total_Genome_xl
  substrate = tmp$V1[!grepl("^U", tmp$V1)]
  if (length(substrate)==0){
    substrate="WTLS54_sub"
  } else {
    substrate=substrate[[1]]
  }
  print(substrate)
  du_df = data.frame(V1=substrate,V2=seq(-1,2000), V3=seq(-1,2000), V4=".",V5=0,V6="+",ss=row$ss, prot=row$prot,rep=row$rep, sub=row$sub, condition=row$condition, complex=row$complex, genome_xls=row$total_Genome_xl )
  ne = full_join(tmp,du_df,by=c("V1","V2","V3", "V4","V5","V6","ss", "prot","rep", "sub", "condition", "complex", "genome_xls")) %>% arrange(-V2)
  ne = ne %>% ungroup() %>% dplyr::group_by(V1,V2,ss,prot,rep,sub,condition,complex,genome_xls) %>% filter(V5==max(V5))
  print(ne)
  all[[i]] = ne
}
all = do.call(rbind, all)

# correct coordinates - make everything relative to 3'SS
allc <- all %>% ungroup() %>% mutate(V2 = ifelse(V1=="ACT1_dG" & complex=="Cstar", V2-450,
                                                 ifelse(V1=="dG_long_sub" & complex=="Cstar", V2-258,
                                                        ifelse(V1=="WTLS54_sub" & complex=="P" & V2>20, V2-96,
                                                               ifelse(V1=="ACT1_WT-MS2" & complex=="P" & V2>27, V2-310, V2))))) #%>%
  #mutate(V2 = ifelse(V1=="WTLS54_sub" & complex=="P", V2-20,
                     ifelse(V1=="ACT1_WT-MS2" & complex=="P", V2-27, V2)))

# remove the 5 branch point nucleotides for C star complex
alld = allc %>% filter(ifelse(V1=="ACT1_dG" & complex=="Cstar", V2!="-47" & V2!="-46" & V2!="-45" & V2!="-44" & V2!="-43",
                              ifelse(V1=="dG_long_sub" & complex=="Cstar", V2!="-29" & V2!="-28" & V2!="-27" & V2!="-26" & V2!="-25", V2 < 100000000000)))

# to merge replicates
ff = alld %>% ungroup() %>% 
  dplyr::group_by(ss,prot,sub,condition,complex) %>%
  dplyr::mutate(lib_size=sum(V5)) %>%
  filter(!grepl("^U", V1)) %>%
  ungroup() %>% dplyr::group_by(V1,V2,ss,prot,sub,condition,complex,lib_size) %>%
  dplyr::summarise(V5=sum(V5), genome_xls=sum(genome_xls)) %>% ungroup()

ff$norm = ff$V5/ff$genome_xls
ff[is.na(ff)] <- 0


# to show crosslinking distributions
reps_merg = ff
count_reps_merg = reps_merg %>% ungroup() %>% 
  group_by(V1,sub,condition,ss,prot, complex) %>% summarise(total=sum(V5)) %>% ungroup() %>% 
  mutate(V1=paste0(V1, " snRNA") %>% gsub(".*dG snRNA|.*sub snRNA|.*MS2 snRNA","pre-mRNA substrate", .), 
         prot=paste0(prot, " Prp22"))

ggplot(count_reps_merg %>% filter(condition=="UV"), aes(x=interaction(sub,ss,prot,condition, complex), y=total, fill=V1)) + geom_col() + 
  theme_few() + ylab("number of crosslinks") + coord_flip() + facet_wrap(~complex+sub, ncol=2, scales="free") +
  ggtitle("all reps together") + xlab("") + theme(legend.title = element_blank())
ggsave("prp22_all_reps_together_-snRNA-substrate-proportions.eps", height=5, width=12)

# to show snRNA enrichment
reps_merg_sn = ff
count_reps_merg_sn = reps_merg_sn %>% ungroup() %>%
  mutate(V1=paste0(V1, " snRNA") %>% gsub(".*dG snRNA|.*sub snRNA|.*MS2 snRNA","pre-mRNA substrate", .), 
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

ff_uv <- ff %>% filter(condition=="UV") %>% filter(grepl("sub|ACT1", V1)) %>%
  filter(V2 > -50 & V2 < 50)

ggplot(ff_uv, aes(x=V2, y=ctrl_minus, group=interaction(ss,sub,prot,complex),color=interaction(ss,prot,complex))) +
  geom_vline(xintercept = c(0, -26, -44), linetype=5) +
  scale_x_continuous(limits =c(-50,50), breaks=seq(-50,50,5), labels=c(insert(seq(-50,50,10), ats=2:length(seq(-50,50,10)), ""))) +
  scale_color_manual(values = cbp1) +
  geom_line(aes(y=ctrl_minus),size=1) +
  ylab("normalised score") +
  xlab("") +
  theme_few() +
  facet_wrap(~sub+complex, scales="free") +
  theme(text = element_text(size=20))
ggsave("PRP22profile_repCombined_ctrlminus_yeastGenomeNorm.eps", height=6, width=13)

###### Replicates separately, also UV and ctrl separately

# to merge replicates
# to merge replicates
ff = alld %>% ungroup() %>% 
  dplyr::group_by(ss,prot,sub,condition,complex,rep) %>%
  dplyr::mutate(lib_size=sum(V5)) %>%
  filter(!grepl("^U", V1)) %>%
  ungroup() %>% dplyr::group_by(V1,V2,ss,prot,sub,condition,complex,lib_size,rep) %>%
  dplyr::summarise(V5=sum(V5), genome_xls=sum(genome_xls)) %>% ungroup()

ff$norm = ff$V5/ff$genome_xls
ff[is.na(ff)] <- 0

ff <- ff %>% dplyr::group_by(sub,ss,prot,condition, complex,rep) %>% dplyr::mutate(smoothed= smth.gaussian(norm, window=10))
ff$condition = factor(ff$condition, levels=c("UV", "ctrl"))
ff_uv <- ff %>% filter(grepl("sub|ACT1", V1)) %>%
  filter(V2 > -50 & V2 < 50)

ggplot(ff_uv, aes(x=V2, y=smoothed, group=interaction(condition,rep,prot,complex),color=interaction(rep,prot,complex), linetype=condition)) +
  geom_vline(xintercept = c(0, -26, -44), linetype=5) +
  scale_x_continuous(limits =c(-50,50), breaks=seq(-50,50,5), labels=c(insert(seq(-50,50,10), ats=2:length(seq(-50,50,10)), ""))) +
  scale_colour_manual(values = c("#036630","#0e8943","#97c681","#cbedb9","#5d598e","#908ece")) +
  geom_line(aes(y=smoothed),size=1) +
  ylab("normalised score") +
  xlab("") +
  theme_few() +
  facet_wrap(~sub+complex, scales="free") +
  theme(text = element_text(size=20))
ggsave("PRP22profile_repSep_yeastGenomeNorm.eps", height=6, width=13)



#FIN
options(scipen = 999)
plot_all = fin %>% filter(grepl("36_dn|42_dn|41|38", exp)) %>% filter(V2 > -75 & V2 < 75)
ggplot(plot_all, aes(x=V2, group=exp, colour=exp)) +
  geom_vline(xintercept = c(0), linetype=5) +
  #geom_col(fill="grey") +
  #geom_line(aes(y=smth.gaussian(V5, window=10)), size=1) #+
  #geom_smooth(aes(y=V5), method = "loess", span = 0.2, se=FALSE) #+
  geom_line(aes(y=smth.gaussian(libsize_norm, window=10)), size=1) +
  #geom_smooth(aes(y=libsize_norm), se=F,method="glm", formula=y~ns(x,15)) +
  scale_x_continuous(limits =c(-50,50), breaks=seq(-50,50,5), labels=c(insert(seq(-50,50,10), ats=2:length(seq(-50,50,10)), ""))) +
  ylab("cDNAs") +
  xlab("") +
  facet_wrap(complex~sub, scales="free_y") + 
  theme_few()
ggsave("fig5_p_to_ctstar_comparsion_combined_gaussian.eps", height=6, width=12)

#REPLICATE
options(scipen = 999)
plot_all = fin %>% filter(grepl("28_dn|40_dn|39|37", exp)) %>% filter(V2 > -75 & V2 < 75)
ggplot(ff_uv, aes(x=V2, colour=exp)) +
  geom_vline(xintercept = c(0), linetype=5) +
  #geom_col(fill="grey") +
  #geom_line(aes(y=smth.gaussian(V5, window=10)), size=1) #+
  #geom_smooth(aes(y=V5), method = "loess", span = 0.2, se=FALSE) #+
  geom_line(aes(y=smth.gaussian(libsize_norm, window=10)), size=1) +
  #geom_smooth(aes(y=libsize_norm), se=F,method="glm", formula=y~ns(x,15)) +
  scale_x_continuous(limits =c(-50,50), breaks=seq(-50,50,5), labels=c(insert(seq(-50,50,10), ats=2:length(seq(-50,50,10)), ""))) +
  ylab("cDNAs") +
  xlab("") +
  facet_wrap(complex~sub, scales="free_y") + 
  theme_few()
ggsave("fig5_p_to_ctstar_comparsion_combined_gaussian_REPLICATE.eps", height=6, width=12)
