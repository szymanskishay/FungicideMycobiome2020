#Alpha Diversity Tests
#This portion uses:
#library(agricolae)
#library(phyloseq)
#library(vegan)
#library(ggplot2)
#library(ggpubr)
###### Alpha Diversity Tests (Whole) ####
sd<-data.frame(sample_data(spraydata)) #Retrieving metadata from phyloseq
sd$shannon<-vegan::diversity(otu_table(spraydata), MARGIN=2, index="shannon") #shannon index for whole data set
sd$simpson<-vegan::diversity(otu_table(spraydata), MARGIN=2, index="simpson") #simpson
sd$fishera<-fisher.alpha(otu_table(spraydata), MARGIN=2) #fisher

shannon_anova_tissue<-aov(shannon~Tissue, sd) #asessing variance and the like
shannon_anova_spray<-aov(shannon~Spray, sd)
simpson_anova_tissue<-aov(simpson~Tissue, sd)
simpson_anova_spray<-aov(simpson~Spray, sd)
fisher_anova_tissue<-aov(fishera~Tissue, sd)
fisher_anova_spray<-aov(fishera~Spray, sd)


tukey_Shannon_Spray<-HSD.test(shannon_anova_spray,  "Spray", group=TRUE, unbalanced=TRUE) #separation of alpha diversity based on spray type
tukey_Simpson_Spray<-HSD.test(simpson_anova_spray, "Spray", group=TRUE, unbalanced=TRUE) #unbalanced=TRUE due to disparate sample counts in a few instances
tukey_fisher_Spray<-HSD.test(fisher_anova_spray, "Spray", group=TRUE, unbalanced=TRUE)

tukey_Shannon_tissue<-HSD.test(shannon_anova_tissue,  "Tissue", group=TRUE, unbalanced=TRUE) #separation of alpha diversity based on tissue type
tukey_Simpson_tissue<-HSD.test(simpson_anova_tissue, "Tissue", group=TRUE, unbalanced=TRUE) #unbalanced=TRUE due to disparate sample counts
tukey_fisher_tissue<-HSD.test(fisher_anova_tissue, "Tissue", group=TRUE, unbalanced=TRUE)


#make into new objects for easier plotting
tk.Shannon.spray<-cbind(tukey_Shannon_Spray$means[,2:8], tukey_Shannon_Spray$groups[order(row.names((tukey_Shannon_Spray$groups))),])
tk.Simpson.spray<-cbind(tukey_Simpson_Spray$means[,2:8], tukey_Simpson_Spray$groups[order(row.names((tukey_Simpson_Spray$groups))),])
tk.fisher.spray<-cbind(tukey_fisher_Spray$means[,2:8], tukey_fisher_Spray$groups[order(row.names((tukey_fisher_Spray$groups))),])

tk.Shannon.tissue<-cbind(tukey_Shannon_tissue$means[,2:8], tukey_Shannon_tissue$groups[order(row.names((tukey_Shannon_tissue$groups))),])
tk.Simpson.tissue<-cbind(tukey_Simpson_tissue$means[,2:8], tukey_Simpson_tissue$groups[order(row.names((tukey_Simpson_tissue$groups))),])
tk.fisher.tissue<-cbind(tukey_fisher_tissue$means[,2:8], tukey_fisher_tissue$groups[order(row.names((tukey_fisher_tissue$groups))),])




##### Alpha Diversity Tests (Tissue Restricted) ####
pulpdata<-subset_samples(spraydata, Tissue%in%c("Pulp")) #split by tissue type
skindata<-subset_samples(spraydata, Tissue%in%c("Skin"))
sd.pulp<-data.frame(sample_data(pulpdata)) #Retrieving metadata from phyloseq
sd.pulp$shannon<-vegan::diversity(otu_table(pulpdata), MARGIN=2, index="shannon")
sd.pulp$simpson<-vegan::diversity(otu_table(pulpdata), MARGIN=2, index="simpson")
sd.pulp$fishera<-fisher.alpha(otu_table(pulpdata), MARGIN=2)

shannon_anova_spray.pulp<-aov(shannon~Spray, sd.pulp)
simpson_anova_spray.pulp<-aov(simpson~Spray, sd.pulp)
fisher_anova_spray.pulp<-aov(fishera~Spray, sd.pulp)


tukey_Shannon_spray.pulp<-HSD.test(shannon_anova_spray.pulp,  "Spray", group=TRUE, unbalanced=TRUE)
tukey_Simpson_spray.pulp<-HSD.test(simpson_anova_spray.pulp, "Spray", group=TRUE, unbalanced=TRUE)
tukey_fisher_spray.pulp<-HSD.test(fisher_anova_spray.pulp, "Spray", group=TRUE, unbalanced=TRUE)


tk.Shannon.spray.pulp<-cbind(tukey_Shannon_spray.pulp$means[,2:8], tukey_Shannon_spray.pulp$groups[order(row.names((tukey_Shannon_spray.pulp$groups))),])
tk.Simpson.spray.pulp<-cbind(tukey_Simpson_spray.pulp$means[,2:8], tukey_Simpson_spray.pulp$groups[order(row.names((tukey_Simpson_spray.pulp$groups))),])
tk.fisher.spray.pulp<-cbind(tukey_fisher_spray.pulp$means[,2:8], tukey_fisher_spray.pulp$groups[order(row.names((tukey_fisher_spray.pulp$groups))),])



skindata<-subset_samples(spraydata, Tissue%in%c("Skin"))
sd.skin<-data.frame(sample_data(skindata)) #Retrieving metadata from phyloseq
sd.skin$shannon<-vegan::diversity(otu_table(skindata), MARGIN=2, index="shannon")
sd.skin$simpson<-vegan::diversity(otu_table(skindata), MARGIN=2, index="simpson")
sd.skin$fishera<-fisher.alpha(otu_table(skindata), MARGIN=2)

shannon_anova_spray.skin<-aov(shannon~Spray, sd.skin)
simpson_anova_spray.skin<-aov(simpson~Spray, sd.skin)
fisher_anova_spray.skin<-aov(fishera~Spray, sd.skin)

tukey_Shannon_spray.skin<-HSD.test(shannon_anova_spray.skin,  "Spray", group=TRUE, unbalanced=TRUE)
tukey_Simpson_spray.skin<-HSD.test(simpson_anova_spray.skin, "Spray", group=TRUE, unbalanced=TRUE)
tukey_fisher_spray.skin<-HSD.test(fisher_anova_spray.skin, "Spray", group=TRUE, unbalanced=TRUE)

tk.Shannon.spray.skin<-cbind(tukey_Shannon_spray.skin$means[,2:8], tukey_Shannon_spray.skin$groups[order(row.names((tukey_Shannon_spray.skin$groups))),])
tk.Simpson.spray.skin<-cbind(tukey_Simpson_spray.skin$means[,2:8], tukey_Simpson_spray.skin$groups[order(row.names((tukey_Simpson_spray.skin$groups))),])
tk.fisher.spray.skin<-cbind(tukey_fisher_spray.skin$means[,2:8], tukey_fisher_spray.skin$groups[order(row.names((tukey_fisher_spray.skin$groups))),])
######MAKING PLOTS #####

shannon.tissue.plot<-ggplot(tk.Shannon.tissue, aes(x=row.names(tk.Shannon.tissue),
                                                     ymin=Min,
                                                     lower=Q25,
                                                     middle=Q50,
                                                     upper=Q75,
                                                     ymax=Max,
                                                     fill=row.names(tk.Shannon.tissue)))+
  theme_classic()+
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Q75, label=groups), size=5, vjust=-1, nudge_x = -0.15)+
  xlab("Tissue")+
  ylab("Richness (H')")+
  labs(tag="A")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=16, face="bold"),axis.title.y=element_text(size=16, face="bold"),legend.title=element_text(size=18, face="bold"),  legend.text=element_text(size=16))+
  guides(fill=guide_legend(title="Tissue"))+
  theme(plot.tag=element_text(size=16, face="bold"))

shannon.spray.plot<-ggplot(tk.Shannon.spray, aes(x=row.names(tk.Shannon.spray),
                                                  ymin=Min,
                                                  lower=Q25,
                                                  middle=Q50,
                                                  upper=Q75,
                                                  ymax=Max,
                                                  fill=row.names(tk.Shannon.spray)))+
  theme_classic()+ scale_fill_manual(values=c("Abound"="#1372C5", "Abound+NuFilm"="#e67e00", "Conventional"="#5C1F94",                                "Double Nickel"="#1CC171","Lifegard"="#B3220F", "Untreated"="#FFFFFF"))+
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Q75, label=groups), size=5, vjust=-1, nudge_x = -0.2)+
  xlab("Treatment")+
  ylab("Richness (H')")+
  labs(tag="A")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=16, face="bold"),axis.title.y=element_text(size=16, face="bold"),legend.title=element_text(size=18, face="bold"),  legend.text=element_text(size=16))+
  guides(fill=guide_legend(title="Spray"))+
  theme(plot.tag=element_text(size=16, face="bold"))


simpson.tissue.plot<-ggplot(tk.Simpson.tissue, aes(x=row.names(tk.Simpson.tissue),
                                                   ymin=Min,
                                                   lower=Q25,
                                                   middle=Q50,
                                                   upper=Q75,
                                                   ymax=Max,
                                                   fill=row.names(tk.Simpson.tissue)))+
  theme_classic()+
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Q75, label=groups), size=5, vjust=-1, nudge_x = -0.15)+
  xlab("Tissue")+
  ylab("Simpson's Diversity")+
  labs(tag="B")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=16, face="bold"),axis.title.y=element_text(size=16, face="bold"),legend.title=element_text(size=18, face="bold"),  legend.text=element_text(size=16))+
  guides(fill=guide_legend(title="Tissue"))+
  theme(plot.tag=element_text(size=16, face="bold"))

simpson.spray.plot<-ggplot(tk.Simpson.spray, aes(x=row.names(tk.Simpson.spray),
                                                 ymin=Min,
                                                 lower=Q25,
                                                 middle=Q50,
                                                 upper=Q75,
                                                 ymax=Max,
                                                 fill=row.names(tk.Simpson.spray)))+
  theme_classic()+ scale_fill_manual(values=c("Abound"="#1372C5", "Abound+NuFilm"="#e67e00", "Conventional"="#5C1F94",                                "Double Nickel"="#1CC171","Lifegard"="#B3220F", "Untreated"="#FFFFFF"))+
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Q75, label=groups), size=5, vjust=-1, nudge_x = -0.2)+
  xlab("Treatment")+
  ylab("Simpson's Diversity")+
  labs(tag="B")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=16, face="bold"),axis.title.y=element_text(size=16, face="bold"),legend.title=element_text(size=18, face="bold"),  legend.text=element_text(size=16))+
  guides(fill=guide_legend(title="Spray"))+
  theme(plot.tag=element_text(size=16, face="bold"))

fisher.tissue.plot<-ggplot(tk.fisher.tissue, aes(x=row.names(tk.fisher.tissue),
                                                   ymin=Min,
                                                   lower=Q25,
                                                   middle=Q50,
                                                   upper=Q75,
                                                   ymax=Max,
                                                   fill=row.names(tk.fisher.tissue)))+
  theme_classic()+ 
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Q75, label=groups), size=5, vjust=-1, nudge_x = -0.15)+
  xlab("Tissue")+
  ylab("Fisher's alpha")+
  labs(tag="C")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=16, face="bold"),axis.title.y=element_text(size=16, face="bold"),legend.title=element_text(size=18, face="bold"),  legend.text=element_text(size=16))+
  guides(fill=guide_legend(title="Tissue"))+
  theme(plot.tag=element_text(size=16, face="bold"))

fisher.spray.plot<-ggplot(tk.fisher.spray, aes(x=row.names(tk.fisher.spray),
                                                 ymin=Min,
                                                 lower=Q25,
                                                 middle=Q50,
                                                 upper=Q75,
                                                 ymax=Max,
                                                 fill=row.names(tk.fisher.spray)))+
  theme_classic()+ scale_fill_manual(values=c("Abound"="#1372C5", "Abound+NuFilm"="#e67e00", "Conventional"="#5C1F94",                                "Double Nickel"="#1CC171","Lifegard"="#B3220F", "Untreated"="#FFFFFF"))+
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Q75, label=groups), size=5, vjust=-1, nudge_x = -0.2)+
  labs(tag="C")+
  xlab("Treatment")+
  ylab("Fisher's alpha")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=16, face="bold"),axis.title.y=element_text(size=16, face="bold"),legend.title=element_text(size=18, face="bold"),  legend.text=element_text(size=16))+
  guides(fill=guide_legend(title="Spray"))+
  theme(plot.tag=element_text(size=16, face="bold"))

tiff(filename="Figure_S1.tiff", width=750, height=500, units="px")
ggarrange(shannon.tissue.plot,simpson.tissue.plot,fisher.tissue.plot, nrow=1, legend="bottom", common.legend = TRUE)
dev.off()
tiff(filename="Figure_S2.tiff", width=1500, height=500, units="px")
ggarrange(shannon.spray.plot, simpson.spray.plot, fisher.spray.plot, nrow=1,legend="bottom", common.legend = TRUE)
dev.off()



shannon.skin.plot<-ggplot(tk.Shannon.spray.skin, aes(x=row.names(tk.Shannon.spray.skin),
                                     ymin=Min,
                                     lower=Q25,
                                     middle=Q50,
                                     upper=Q75,
                                     ymax=Max,
                                     fill=row.names(tk.Shannon.spray.skin)))+
  theme_classic()+ scale_fill_manual(values=c("Abound"="#1372C5", "Abound+NuFilm"="#e67e00", "Conventional"="#5C1F94",                                "Double Nickel"="#1CC171","Lifegard"="#B3220F", "Untreated"="#FFFFFF"))+
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Q75, label=groups), size=5, vjust=-1, nudge_x = -0.15)+
  labs(tag="A")+
  xlab("Treatment")+
  ylab("Richness (H')")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=16, face="bold"),axis.title.y=element_text(size=16, face="bold"),legend.title=element_text(size=18, face="bold"),  legend.text=element_text(size=16))+
  guides(fill=guide_legend(title="Treatment"))+
  theme(plot.tag=element_text(size=16, face="bold"))

shannon.pulp.plot<-ggplot(tk.Shannon.spray.pulp, aes(x=row.names(tk.Shannon.spray.pulp),
                                                     ymin=Min,
                                                     lower=Q25,
                                                     middle=Q50,
                                                     upper=Q75,
                                                     ymax=Max,
                                                     fill=row.names(tk.Shannon.spray.pulp)))+
  theme_classic()+ scale_fill_manual(values=c("Abound"="#1372C5", "Abound+NuFilm"="#e67e00", "Conventional"="#5C1F94",                                "Double Nickel"="#1CC171","Lifegard"="#B3220F", "Untreated"="#FFFFFF"))+
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Q75, label=groups), size=5, vjust=-1, nudge_x = -0.15)+
  labs(tag="D")+
  xlab("Treatment")+
  ylab("Richness (H')")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=16, face="bold"),axis.title.y=element_text(size=16, face="bold"),legend.title=element_text(size=18, face="bold"),  legend.text=element_text(size=16))+
  guides(fill=guide_legend(title="Treatment"))+
  theme(plot.tag=element_text(size=16, face="bold"))

Simpson.skin.plot<-ggplot(tk.Simpson.spray.skin, aes(x=row.names(tk.Simpson.spray.skin),
                                                     ymin=Min,
                                                     lower=Q25,
                                                     middle=Q50,
                                                     upper=Q75,
                                                     ymax=Max,
                                                     fill=row.names(tk.Simpson.spray.skin)))+
  theme_classic()+ scale_fill_manual(values=c("Abound"="#1372C5", "Abound+NuFilm"="#e67e00", "Conventional"="#5C1F94",                                "Double Nickel"="#1CC171","Lifegard"="#B3220F", "Untreated"="#FFFFFF"))+
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Q75, label=groups), size=5, vjust=-1, nudge_x = -0.15)+
  labs(tag="B")+
  xlab("Treatment")+
  ylab("Simpson's Diversity")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=16, face="bold"),axis.title.y=element_text(size=16, face="bold"),legend.title=element_text(size=18, face="bold"),  legend.text=element_text(size=16))+
  guides(fill=guide_legend(title="Treatment"))+
  theme(plot.tag=element_text(size=16, face="bold"))

Simpson.pulp.plot<-ggplot(tk.Simpson.spray.pulp, aes(x=row.names(tk.Simpson.spray.pulp),
                                                     ymin=Min,
                                                     lower=Q25,
                                                     middle=Q50,
                                                     upper=Q75,
                                                     ymax=Max,
                                                     fill=row.names(tk.Simpson.spray.pulp)))+
  theme_classic()+ scale_fill_manual(values=c("Abound"="#1372C5", "Abound+NuFilm"="#e67e00", "Conventional"="#5C1F94",                                "Double Nickel"="#1CC171","Lifegard"="#B3220F", "Untreated"="#FFFFFF"))+
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Q75, label=groups), size=5, vjust=-1, nudge_x = -0.15)+
  labs(tag="E")+
  xlab("Treatment")+
  ylab("Simpson's Diversity")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=16, face="bold"),axis.title.y=element_text(size=16, face="bold"),legend.title=element_text(size=18, face="bold"),  legend.text=element_text(size=16))+
  guides(fill=guide_legend(title="Treatment"))+
  theme(plot.tag=element_text(size=16, face="bold"))


fisher.skin.plot<-ggplot(tk.fisher.spray.skin, aes(x=row.names(tk.fisher.spray.skin),
                                                     ymin=Min,
                                                     lower=Q25,
                                                     middle=Q50,
                                                     upper=Q75,
                                                     ymax=Max,
                                                     fill=row.names(tk.fisher.spray.skin)))+
  theme_classic()+ scale_fill_manual(values=c("Abound"="#1372C5", "Abound+NuFilm"="#e67e00", "Conventional"="#5C1F94",                                "Double Nickel"="#1CC171","Lifegard"="#B3220F", "Untreated"="#FFFFFF"))+
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Q75, label=groups), size=5, vjust=-1, nudge_x = -0.15)+
  labs(tag="C")+
  xlab("Treatment")+
  ylab("Fisher's alpha")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=16, face="bold"),axis.title.y=element_text(size=16, face="bold"),legend.title=element_text(size=18, face="bold"),  legend.text=element_text(size=16))+
  guides(fill=guide_legend(title="Treatment"))+
  theme(plot.tag=element_text(size=16, face="bold"))

fisher.pulp.plot<-ggplot(tk.fisher.spray.pulp, aes(x=row.names(tk.fisher.spray.pulp),
                                                     ymin=Min,
                                                     lower=Q25,
                                                     middle=Q50,
                                                     upper=Q75,
                                                     ymax=Max,
                                                     fill=row.names(tk.fisher.spray.pulp)))+
  theme_classic()+ scale_fill_manual(values=c("Abound"="#1372C5", "Abound+NuFilm"="#e67e00", "Conventional"="#5C1F94",                                "Double Nickel"="#1CC171","Lifegard"="#B3220F", "Untreated"="#FFFFFF"))+
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Q75, label=groups), size=5, vjust=-1, nudge_x = -0.15)+
  labs(tag="F")+
  xlab("Treatment")+
  ylab("Fisher's alpha")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=16, face="bold"),axis.title.y=element_text(size=16, face="bold"),legend.title=element_text(size=18, face="bold"),  legend.text=element_text(size=16))+
  guides(fill=guide_legend(title="Treatment"))+
  theme(plot.tag=element_text(size=16, face="bold"))


tiff(filename="Figure_1.tiff", width=1500, height=750, units="px")
ggarrange(shannon.skin.plot, Simpson.skin.plot, fisher.skin.plot, shannon.pulp.plot, Simpson.pulp.plot, fisher.pulp.plot, nrow=2, ncol=3,legend="bottom", common.legend = TRUE)
dev.off()
###Output section used to help ease table ####


#m<-c("Treatment", "Stat", "STD", "Group")
#output.tk.Simpson.spray<-cbind(row.names(tk.Simpson.spray),tk.Simpson.spray$simpson, tk.Simpson.spray$std, tk.Simpson.spray$groups)
#output.tk.Shannon.spray<-cbind(row.names(tk.Shannon.spray),tk.Shannon.spray$shannon, tk.Shannon.spray$std, tk.Shannon.spray$groups)
#output.tk.fisher.spray<-cbind(row.names(tk.fisher.spray),tk.fisher.spray$fisher, tk.fisher.spray$std, tk.fisher.spray$groups)
#output.spray.whole<-rbind(m, output.tk.Shannon.spray, m, output.tk.Simpson.spray, m, output.tk.fisher.spray)
#write.table(output.spray.whole, "Whole Spray.txt")
#output.tk.Simpson.spray.pulp<-cbind(row.names(tk.Simpson.spray.pulp),tk.Simpson.spray.pulp$simpson, tk.Simpson.spray.pulp$std, tk.Simpson.spray.pulp$groups)
#output.tk.Shannon.spray.pulp<-cbind(row.names(tk.Shannon.spray.pulp),tk.Shannon.spray.pulp$shannon, tk.Shannon.spray.pulp$std, tk.Shannon.spray.pulp$groups)
#output.tk.fisher.spray.pulp<-cbind(row.names(tk.fisher.spray.pulp),tk.fisher.spray.pulp$fisher, tk.fisher.spray.pulp$std, tk.fisher.spray.pulp$groups)
#output.spray.pulp<-rbind(m, output.tk.Shannon.spray.pulp, m, output.tk.Simpson.spray.pulp, m, output.tk.fisher.spray.pulp)
#write.table(output.spray.pulp, "Pulp Spray.txt")
#output.tk.Simpson.spray.skin<-cbind(row.names(tk.Simpson.spray.skin),tk.Simpson.spray.skin$simpson, tk.Simpson.spray.skin$std, tk.Simpson.spray.skin$groups)
#output.tk.Shannon.spray.skin<-cbind(row.names(tk.Shannon.spray.skin),tk.Shannon.spray.skin$shannon, tk.Shannon.spray.skin$std, tk.Shannon.spray.skin$groups)
#output.tk.fisher.spray.skin<-cbind(row.names(tk.fisher.spray.skin),tk.fisher.spray.skin$fisher, tk.fisher.spray.skin$std, tk.fisher.spray.skin$groups)
#output.spray.skin<-rbind(m, output.tk.Shannon.spray.skin, m, output.tk.Simpson.spray.skin, m, output.tk.fisher.spray.skin)
#write.table(output.spray.skin, "Skin Spray.txt")


