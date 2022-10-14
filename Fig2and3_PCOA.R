#PCA/Permanova plots
#USE CODE "MAKE-PHYLOSEQ-OBJECT.R" FIRST
# Begin Work Proper----
ps.spray<-ps.noncontam_spray #shorten name because easy (also resetting it for this script)
# this line of code removes any otus without reads so that only real otus are being counted
otu_table(ps.spray) <- otu_table(ps.spray)[which(rowSums(otu_table(ps.spray)) >= 1),]
##BIG STEP: SUBSETTING. This doc is here basically for this part.
ps.spray.naut<- subset_samples(ps.spray, origin%in%c("Blueberry")) #removes controls, keeps Fast Xtract
ps.spray<- subset_samples(ps.spray.naut, Extraction%in%c("Kingfisher", "Chelex")) #Kingfisher+Chelex will be used for most data, keeping name shorter for better downstream applications
spray.pulp<- subset_samples(ps.spray, Tissue%in%c("Pulp"))
spray.skin<- subset_samples(ps.spray, Tissue%in%c("Skin"))
#Load in each comparison type----
ps.spray_pulp<- spray.pulp
ps.spray_skin<- spray.skin


###### Whole-Figure 2
ps.spray.n = phyloseq_to_metagenomeSeq(ps.spray) #n for normalized
p_biom_spray<-cumNormStat(ps.spray.n)
biom_quant_spray<-cumNorm(ps.spray.n, p=p_biom_spray)
normFactors(biom_quant_spray)
ps.spray.nf <-MRcounts(biom_quant_spray, norm=T)
#create physeq object with normalized otu table
otu_table(ps.spray) <- otu_table(ps.spray.nf, taxa_are_rows = TRUE)

spray.ord = ordinate(ps.spray, method ="PCoA", distance="bray")
# this function will visualize the ordination along the two axes that account for the most variation
spray.pca = plot_ordination(ps.spray, spray.ord, color="Tissue") + 
  #labs(title="Pulp") + #the previous line also had shape for a factor, which is good.
  theme_classic()+
  theme(legend.title=element_text(size=18, face="bold"), legend.text=element_text(size=16)) +
  geom_point(size=3, alpha=0.9) + 
  geom_point(size=3, shape=1, color="Black")+
  scale_color_manual(values=c("Pulp"="Green", "Skin"="Blue"))+
  theme(axis.title=element_text(size=16), axis.text=element_text(size=14))+
  stat_ellipse(type="norm", linetype=1)
plot(spray.pca)
#Call Tiff
tiff("Figure_2.tiff", width=500, height=400, units="px")
#Run Plot
plot(spray.pca)
#Close Tiff
dev.off()

##Pulp##----
ps.spray_pulp.n = phyloseq_to_metagenomeSeq(ps.spray_pulp) #n for normalized
p_biom_spray_pulp<-cumNormStat(ps.spray_pulp.n)
biom_quant_spray_pulp<-cumNorm(ps.spray_pulp.n, p=p_biom_spray_pulp)
normFactors(biom_quant_spray_pulp)
ps.spray_pulp.nf <-MRcounts(biom_quant_spray_pulp, norm=T) 
#create physeq object with normalized otu table
otu_table(ps.spray_pulp) <- otu_table(ps.spray_pulp.nf, taxa_are_rows = TRUE) #note: This overwrites the OTU table that is currently at ps.spray_pulp

spray_pulp.ord = ordinate(ps.spray_pulp, method ="PCoA", distance="bray")
# this function will visualize the ordination along the two axes that account for the most variation for pulp samples
spray_pulp.pca = plot_ordination(ps.spray_pulp, spray_pulp.ord, color="Spray") + #Color represents the treatment (Spray)
  theme_classic()+
 labs(title="Pulp", tag="A") + 
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0), legend.title=element_text(size=18, face="bold"), legend.text=element_text(size=16)) +
  theme(axis.title=element_text(size=16), axis.text=element_text(size=14))+
   geom_point(size=3, alpha=0.9) + 
  geom_point(size=3, shape=1, color="Black")+ #adds outline for circles
  scale_color_manual(values=c("Abound"="#1372C5", "Abound+NuFilm"="#e67e00", "Conventional"="#5C1F94", 
                              "Double Nickel"="#1CC171","Lifegard"="#B3220F", "Untreated"="#000000"), name="Treatment")+
  stat_ellipse(type="norm", linetype=1)+ #treatment ellipses to better visualize groupings
  theme(plot.tag= element_text(size=20))
plot(spray_pulp.pca) #lets you see this first


##Skin##----
ps.spray_skin.n = phyloseq_to_metagenomeSeq(ps.spray_skin) #n for normalized
p_biom_spray_skin<-cumNormStat(ps.spray_skin.n)
biom_quant_spray_skin<-cumNorm(ps.spray_skin.n, p=p_biom_spray_skin)
normFactors(biom_quant_spray_skin)
ps.spray_skin.nf <-MRcounts(biom_quant_spray_skin, norm=T) 
#create physeq object with normalized otu table
otu_table(ps.spray_skin) <- otu_table(ps.spray_skin.nf, taxa_are_rows = TRUE)

spray_skin.ord = ordinate(ps.spray_skin, method ="PCoA", distance="bray")
# this function will visualize the ordination along the two axes that account for the most variation
spray_skin.pca = plot_ordination(ps.spray_skin, spray_skin.ord, color="Spray") + 
  theme_classic()+
  labs(title="Skin", tag="B")+
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0), legend.title=element_text(size=18, face="bold"), legend.text=element_text(size=16)) +
  theme(axis.title=element_text(size=16), axis.text=element_text(size=14))+
  geom_point(size=3, alpha=0.9) + 
  geom_point(size=3, shape=1, color="Black")+
  scale_color_manual(values=c("Abound"="#1372C5", "Abound+NuFilm"="#e67e00", "Conventional"="#5C1F94", 
                              "Double Nickel"="#1CC171","Lifegard"="#B3220F", "Untreated"="#000000"), name="Treatment")+
  stat_ellipse(type="norm", linetype=1)+
  theme(plot.tag= element_text(size=20))
plot(spray_skin.pca)

figure_3_tissue.pca<-ggarrange(spray_pulp.pca, spray_skin.pca, legend="bottom", common.legend=TRUE) #makes unified Figure 2
tiff("Figure_3.tiff", width=1000, height=500, units="px")
plot(figure_3_tissue.pca)
dev.off()
