
###Treatment Bars
#reset ps.spray if needed
ps.spray<- subset_samples(ps.spray.naut, Extraction%in%c("Kingfisher", "Chelex"))
otu_table(ps.spray)[otu_table(ps.spray) <= 4] <- 0 ### This sets any OTUs whose count in a single sample is less than 4 to 0. The purpose of this is to ensure that if there was low level barcode switching that it gets removed.

otu_table(ps.spray) <- otu_table(ps.spray)[which(rowSums(otu_table(ps.spray)) >= 10),]### PCR Errors
# this removes otus which had less than ten total reads to account for potential pcr relics

ps.spray.pulp<- subset_samples(ps.spray, Tissue%in%c("Pulp"))
ps.spray.skin<- subset_samples(ps.spray, Tissue%in%c("Skin"))
##Skin##----
ps.spray.skin.melt<- ps.spray.skin %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at Genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format                       # Filter out low abundance taxa
  arrange(Genus)           # Sort data frame alphabetically by Genus


dat_ps.spray.skinbar <- data.table(ps.spray.skin.melt)


dat_ps.spray.skinbar[(Abundance <= 0.04), Genus:= "Other"] 
write.csv(dat_ps.spray.skinbar, "genus-table-skin.csv")
ps.spray.skinbarplot= ggplot(dat_ps.spray.skinbar, aes(x = Sample, y = Abundance, fill = Genus)) + 
  facet_wrap(~Spray, strip.position = "bottom", scales="free_x") +
  #theme(axis.text.x = element_text(angle = 90))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Aureobasidium"="#fff000",	
                               "Epicoccum"="#1b1be4",	
                               "Filobasidium"="#795f38",	
                               "Heterocephalacria"="#848dc5",	
                               "Metschnikowia"="#89f170",	
                               "Papiliotrema"="#6a6a6a",	
                               "Paraboeremia"="#facc76",	
                               "Sporobolomyces"="#9e1100",	
                               "Taphrina"="#d7791c",	
                               "Vishniacozyma"="#829E92",	
                               "Unassigned"="#ff87dd",	
                               "Other"="#101010"	
  ))+
  # Remove x axis title
  theme_classic()+
  theme(axis.title.x = element_blank()) + 
  theme(legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.25, "cm")) +
  theme(legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(plot.title = element_text(size = 16, hjust = 0, face="bold")) +
  ggtitle("Skin")+
  theme(axis.text.x= element_blank())+
  #theme(axis.text.x = element_text(size = 0, angle = 0, vjust = 0, hjust = 0)) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.y = element_text(angle = 90, size = 16, face = "bold", vjust=-2)) +
  theme(legend.position="right")+
  theme(axis.text.y=element_text(size=14))+
  theme(strip.text = element_text(size = 12, face = "bold")) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Genera > 4%) \n") +
  xlab("")+
  labs(tag="B")+  theme(plot.tag=element_text(size=16, face="bold"))
plot(ps.spray.skinbarplot)

tiff(filename="Figure_4b.tiff", width=1000, height= 480, units="px") #creates figure 4b 
plot(ps.spray.skinbarplot)
dev.off()
##Pulp##----
ps.spray.pulp.bar<- ps.spray.pulp %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at Genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format                       # Filter out low abundance taxa
  arrange(Genus)           # Sort data frame alphabetically by Genus

dat_ps.spray.pulpbar <- data.table(ps.spray.pulp.bar)

dat_ps.spray.pulpbar[(Abundance <= 0.04), Genus:= "Other"] 
write.csv(dat_ps.spray.pulpbar, "pulpgenus.csv")
##
ps.spray.pulpbarplot= ggplot(dat_ps.spray.pulpbar, aes(x = Sample, y = Abundance, fill = Genus)) + 
  facet_wrap(~Spray, strip.position = "bottom", scales="free_x") +
  #theme(axis.text.x = element_text(angle = 90))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c( "Alternaria"="#155c1a",
                                "Aspergillus"="#d8c494",
                                "Aureobasidium"="#fff000",
                                "Botrytis"="#7c1fa9",
                                "Candida"="#24a4bf",
                                "Cladophialophora"="#aa654e",
                                "Coniothyrium"="#059600",
                                "Curvibasidium"="#7b8ad2",
                                "Didymella"="#f0fb73",
                                "Dioszegia"="#c4fe6b",
                                "Epicoccum"="#1b1be4",
                                "Erythrobasidium"="#dfc3ef",
                                "Exophiala"="#ff7826",
                                "Filobasidium"="#795f38",
                                "Genolevuria"="#14dc70",
                                "Hannaella"="#ccc204",
                                "Heterocephalacria"="#848dc5",
                                "Issatchenkia"="#f081f5",
                                "Mortierella"="#fd8c6e",
                                "Naganishia"="#f8381b",
                                "Papiliotrema"="#6a6a6a",
                                "Paraboeremia"="#facc76",
                                "Peltaster"="#fdb7b7",
                                "Phaeococcomyces"="#22ffdd",
                                "Phoma"="#e67e00",
                                "Pichia"="#9f45fc",
                                "Pithomyces"="#799205",
                                "Podosphaeria"="#ff7b6c",
                                "Pulvinula"="#99ea77",
                                "Pyrenochaetopsis"="#fff0da",
                                "Rhodotorula"="#ed3e7e",
                                "Sphaerosporella"="#aa654e",
                                "Sporobolomyces"="#9e1100",
                                "Taphrina"="#d7791c",
                                "Tetraplosphaeria"="#b8b8b8",
                                "Thanatephorus"="#65c5c2",
                                "Trichosporon"="#281137",
                                "Wilcoxina"="#06bf00",
                                "Unassigned"="#ff87dd",
                                "Other"="#101010"))+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  theme_classic()+
  theme(legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.25, "cm")) +
  theme(legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(plot.title = element_text(size = 16, hjust = 0, face="bold")) +
  ggtitle("Pulp")+
  theme(axis.text.x= element_blank())+
  theme(axis.text.y=element_text(size=14))+
  #theme(axis.text.x = element_text(size = 0, angle = 0, vjust = 0, hjust = 0)) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.y = element_text(angle = 90, size = 18, face = "bold", vjust=-2)) +
  theme(legend.position="right")+
  theme(strip.text = element_text(size = 12, face = "bold")) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Genera > 4%) \n") +
  xlab("")+
  labs(tag="A")+
  theme(plot.tag=element_text(size=16, face="bold"))
plot(ps.spray.pulpbarplot)
tiff(filename="Figure_4a.tiff", width= 1200, height= 500, units="px") #creates figure 4a
plot(ps.spray.pulpbarplot)
dev.off()


f4.0<-ggarrange(ps.spray.pulpbarplot, ps.spray.skinbarplot, nrow=2)
tiff(filename="TESTDIFF.tiff", width=1200, height=1200, units="px")
plot(f4.0)
dev.off()

f4<-ggarrange(ps.spray.pulpbarplot, ps.spray.skinbarplot, nrow=2,legend="bottom", common.legend=TRUE)
tiff(filename="TESTBOTTOM.tiff", width=1200, height= 1200, units="px")
plot(f4)
dev.off()

f4.2<-ggarrange(ps.spray.pulpbarplot, ps.spray.skinbarplot, nrow=2,legend="right", common.legend=TRUE)
tiff(filename="TESTRIGHT.tiff", width=1200, height= 1200, units="px")
plot(f4.2)
dev.off()
