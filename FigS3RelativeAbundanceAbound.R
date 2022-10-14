ps.spray<- subset_samples(ps.spray.naut, Extraction%in%c("Kingfisher", "Chelex")) #reset ps.spray
otu_table(ps.spray)[otu_table(ps.spray) <= 4] <- 0 ### This sets any OTUs whose count in a single sample is less than 4 to 0. The purpose of this is to ensure that if there was low level barcode switching that it gets removed.

otu_table(ps.spray) <- otu_table(ps.spray)[which(rowSums(otu_table(ps.spray)) >= 10),]### PCR Errors
# this removes otus which had less than ten total reads to account for potential pcr relics
ps.abound<-subset_samples(ps.spray, Spray %in% c("Abound", "Abound+NuFilm"))
ps.abound.s<-subset_samples(ps.abound, Tissue %in% c("Skin"))
ps.abound.P<-subset_samples(ps.abound, Tissue %in% c("Pulp"))

##### Figure 5a #####
ps.abound.P.bar<- ps.abound.P %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at Genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format                       # Filter out low abundance taxa
  arrange(Genus)           # Sort data frame alphabetically by Genus

dat_ps.abound.P <- data.table(ps.abound.P.bar)

dat_ps.abound.P[(Abundance <= 0.04), Genus:= "Other"] # merge low abundance genera so the plot isn't too noisy. This is up to you.
write.csv(dat_ps.abound.P$Genus, file="ListOfGenu4s_abound_p.csv")
# Plot---
# manually filled colors
ps.abound.P.barplot= ggplot(dat_ps.abound.P, aes(x = Sample, y = Abundance, fill = Genus)) + 
  facet_wrap(~Spray, strip.position = "bottom", scales="free_x") +
  # facet_wrap(~Tissue, strip.position= "bottom", scales="free_x")+
  #theme(axis.text.x = element_text(angle = 90))+
  geom_bar(stat = "identity") +
  theme_classic()+
  scale_fill_manual(values = c("Alternaria"="#155c1a",
                               "Aureobasidium"="#fff000",
                               "Botrytis"="#7c1fa9",
                               "Cladophialophora"="#aa654e",
                               "Coniothyrium"="#059600",
                               "Didymella"="#f0fb73",
                               "Epicoccum"="#1b1be4",
                               "Filobasidium"="#795f38",
                               "Genolevuria"="#14dc70",
                               "Leptosphaerulina"="#000000",
                               "Papiliotrema"="#6a6a6a",
                               "Paraboeremia"="#facc76",
                               "Parastagonospora"="#000000",
                               "Pichia"="#9f45fc",
                               "Pithomyces"="#799205",
                               "Pulvinula"="#99ea77",
                               "Podosphaeria"="#ff7b6c",
                               "Rhodotorula"="#ed3e7e",
                               "Sporobolomyces"="#9e1100",
                               "Trichosporon"="#281137",
                               "Unassigned"="#ff87dd",
                               "Other"="#101010"))+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  theme(legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.25, "cm")) +
  theme(legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14)) +
  theme(strip.text.x = element_text(size = 14, face = "bold")) +
  theme(axis.text.x = element_blank()) +
  theme(plot.title = element_text(size = 16, hjust = 0, face="bold")) +
  ggtitle("Pulp")+
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.y = element_text(angle = 90, size = 16, face = "bold", vjust=-2)) +
  theme(legend.position="right")+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Genera > 4%) \n") +
  xlab("")+
  labs(tag="A")
plot(ps.abound.P.barplot)

tiff(filename="Abound-Pulp.tiff", width=1000, height=480, units="px") #makes Figure 5a
plot(ps.abound.P.barplot)
dev.off()
#### Abound Skin #### Figure 5b ####
ps.abound.s.bar<- ps.abound.s %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at Genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format                       # Filter out low abundance taxa
  arrange(Genus)           # Sort data frame alphabetically by Genus

dat_ps.abound.s <- data.table(ps.abound.s.bar)

dat_ps.abound.s[(Abundance <= 0.04), Genus:= "Other"] # merge low abundance genera so the plot isn't too noisy. This is up to you.
write.csv(dat_ps.abound.s$Genus, file="ListOfGenus_abound_s.csv")
# Plot---
# manually filled colors
ps.abound.s.barplot= ggplot(dat_ps.abound.s, aes(x = Sample, y = Abundance, fill = Genus)) + 
  facet_wrap(~Spray, strip.position = "bottom", scales="free_x") +
  # facet_wrap(~Tissue, strip.position= "bottom", scales="free_x")+
  #theme(axis.text.x = element_text(angle = 90))+
  geom_bar(stat = "identity") +
  theme_classic()+
  scale_fill_manual(values = c("Aureobasidium"="#fff000",
                               "Epicoccum"="#1b1be4",
                               "Filobasidium"="#795f38",
                               "Papiliotrema"="#6a6a6a",
                               "Paraboeremia"="#facc76",
                               "Sporobolomyces"="#9e1100",
                               "Unassigned"="#ff87dd",
                               "Other"="#101010"))+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  theme(legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.25, "cm")) +
  theme(legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14)) +
  theme(strip.text.x = element_text(size = 14, face = "bold")) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y=element_text(size=12))+
  theme(plot.title = element_text(size = 16, hjust = 0, face="bold")) +
  ggtitle("Skin")+
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.y = element_text(angle = 90, size = 16, face = "bold", vjust=-2)) +
  theme(legend.position="right")+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Genera > 4%) \n") +
  xlab("")+
  labs(tag="B")
plot(ps.abound.s.barplot)

tiff(filename="Abound-skin.tiff", width=1000, height=480, units="px") #makes Figure 5b
plot(ps.abound.s.barplot)
dev.off()

fs3<-ggarrange(ps.abound.P.barplot, ps.abound.s.barplot, nrow=2, legend="bottom", common.legend=TRUE)
tiff(filename="Figure_S3.tiff", width=1000, height=1000, units="px")
plot(fs3)
dev.off()
