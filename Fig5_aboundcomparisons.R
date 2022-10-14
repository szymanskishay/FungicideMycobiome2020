
##Abound Comparisons##
ps.spray<- subset_samples(ps.spray.naut, Extraction%in%c("Kingfisher", "Chelex"))
otu_table(ps.spray)[otu_table(ps.spray) <= 4] <- 0 ### This sets any OTUs whose count in a single sample is less than 4 to 0. The purpose of this is to ensure that if there was low level barcode switching that it gets removed.

otu_table(ps.spray) <- otu_table(ps.spray)[which(rowSums(otu_table(ps.spray)) >= 10),]### PCR Errors
# this removes otus which had less than ten total reads to account for potential pcr relics#shorten name because easy (also resetting it for this script)

spray.pulp<- subset_samples(ps.spray, Tissue%in%c("Pulp"))
spray.skin<- subset_samples(ps.spray, Tissue%in%c("Skin"))
pulp.abound.wc<- subset_samples(spray.pulp, Spray%in%c("Untreated", "Abound", "Abound+NuFilm"))
skin.abound.wc<- subset_samples(spray.skin, Spray%in%c("Untreated", "Abound", "Abound+NuFilm"))
###Pulp+Control####
pulp.abound.wc.n = phyloseq_to_metagenomeSeq(pulp.abound.wc) #n for normalized
p_biom_pulp.abound.wc<-cumNormStat(pulp.abound.wc.n)
biom_quant_pulp.abound.wc<-cumNorm(pulp.abound.wc.n, p=p_biom_pulp.abound.wc)
normFactors(biom_quant_pulp.abound.wc)
pulp.abound.wc.nf <-MRcounts(biom_quant_pulp.abound.wc, norm=T) #nf for normFactors?
#create physeq object with normalized otu table
otu_table(pulp.abound.wc) <- otu_table(pulp.abound.wc.nf, taxa_are_rows = TRUE)

pulp.abound.wc.ord = ordinate(pulp.abound.wc, method ="PCoA", distance="bray")
# this function will visualize the ordination along the two axes that account for the most variation
pulp.abound.wc.pca = plot_ordination(pulp.abound.wc, pulp.abound.wc.ord, color="Spray") + 
  theme_classic()+
  labs(title="Abound (Pulp)", tag="A") + #the previous line also had shape for a factor, which is good.
  #theme(plot.title = element_text(size = 12, face = "bold", hjust = -1)) +
  geom_point(size=3, alpha=0.9) + 
  geom_point(size=3, shape=1, color="Black")+
  scale_color_manual(values=c("Abound"="#1372C5", "Abound+NuFilm"="#e67e00", "Untreated"="#000000"))+
  stat_ellipse(type="norm", linetype=1)+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0), legend.title=element_text(size = 18, face="bold"), legend.text=element_text(size=18), axis.title = element_text(size=20), axis.text=element_text(size=16))+
  theme(plot.tag= element_text(size=20))
plot(pulp.abound.wc.pca)
#Call Tiff
tiff("pulp.abound.wc.4-27-22-pca.tiff", width=750, height=550, units="px") #makes figure 5a
#Run Plot
plot(pulp.abound.wc.pca)
#Close Tiff
dev.off()


#PERMANOVA##
#preparing 3 objects to be passed to different parts for the permanova
otu_pulp.abound.wc <- as.data.frame(otu_table(pulp.abound.wc))
taxa_pulp.abound.wc <- as.data.frame(as.matrix(tax_table(pulp.abound.wc)))
metadata_pulp.abound.wc <- as.data.frame(as.matrix(sample_data(pulp.abound.wc)))

#adonis


adonis(t(otu_pulp.abound.wc) ~ Spray, strata=metadata_pulp.abound.wc$Block,  data=metadata_pulp.abound.wc, permutations=9999)

# betadispersion
vegan::vegdist(t(otu_pulp.abound.wc), method="bray") -> dist_otu_pulp.abound.wc #make base info
permdisp_otu_pulp.abound.wc.trt<- betadisper(dist_otu_pulp.abound.wc, metadata_pulp.abound.wc$Spray)#Bdisp Spray
anova(permdisp_otu_pulp.abound.wc.trt, permutations = 9999)#set to write

###skin+Control####
skin.abound.wc.n = phyloseq_to_metagenomeSeq(skin.abound.wc) #n for normalized
p_biom_skin.abound.wc<-cumNormStat(skin.abound.wc.n)
biom_quant_skin.abound.wc<-cumNorm(skin.abound.wc.n, p=p_biom_skin.abound.wc)
normFactors(biom_quant_skin.abound.wc)
skin.abound.wc.nf <-MRcounts(biom_quant_skin.abound.wc, norm=T) #nf for normFactors?
#create physeq object with normalized otu table
otu_table(skin.abound.wc) <- otu_table(skin.abound.wc.nf, taxa_are_rows = TRUE)

skin.abound.wc.ord = ordinate(skin.abound.wc, method ="PCoA", distance="bray")
# this function will visualize the ordination along the two axes that account for the most variation
skin.abound.wc.pca = plot_ordination(skin.abound.wc, skin.abound.wc.ord, color="Spray") + 
  theme_classic()+  
  labs(title="Abound (Skin)", tag="B") + #the previous line also had shape for a factor, which is good.
  #theme(plot.title = element_text(size = 12, face = "bold", hjust = -1)) +
  geom_point(size=3, alpha=0.9) + 
  geom_point(size=3, shape=1, color="Black")+
  scale_color_manual(values=c("Abound"="#1372C5", "Abound+NuFilm"="#e67e00", "Untreated"="#000000"))+
  stat_ellipse(type="norm", linetype=1)+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0), legend.title=element_text(size = 18, face="bold"), legend.text=element_text(size=18), axis.title = element_text(size=20), axis.text=element_text(size=16))+
  theme(plot.tag= element_text(size=20))

plot(skin.abound.wc.pca)
#Call Tiff
tiff("skin.abound.wc-pca.tiff", width=750, height=550, units="px") #makes figure 5b
#Run Plot
plot(skin.abound.wc.pca)
#Close Tiff
dev.off()
figure5aboundpca<-ggarrange(pulp.abound.wc.pca,skin.abound.wc.pca, widths=c(1,1), legend="bottom", common.legend=TRUE)
tiff("figure_5.tiff", width=1000, height=550, units="px")
figure5aboundpca
dev.off()
#makes unified figure 5
#PERMANOVA##
#preparing 3 objects to be passed to different parts for the permanova
otu_skin.abound.wc <- as.data.frame(otu_table(skin.abound.wc))
taxa_skin.abound.wc <- as.data.frame(as.matrix(tax_table(skin.abound.wc)))
metadata_skin.abound.wc <- as.data.frame(as.matrix(sample_data(skin.abound.wc)))

#adonis


adonis(t(otu_skin.abound.wc) ~ Spray, strata=metadata_skin.abound.wc$Block,  data=metadata_skin.abound.wc, permutations=9999)

# betadispersion
vegan::vegdist(t(otu_skin.abound.wc), method="bray") -> dist_otu_skin.abound.wc #make base info
permdisp_otu_skin.abound.wc.trt<- betadisper(dist_otu_skin.abound.wc, metadata_skin.abound.wc$Spray)#Bdisp Spray
anova(permdisp_otu_skin.abound.wc.trt, permutations = 9999)#set to write