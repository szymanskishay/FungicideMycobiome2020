
########################
library(ALDEx2)
Sprame<-as.data.frame(otu_table(spraydata))
Spramet<-t(Sprame)
Spramet2 = Spramet %>% arrange(row.names.data.frame(Spramet))
pulp.frame<-as.data.frame(otu_table(pulpdata))
pulp.frame.1<-pulp.frame
colnames(pulp.frame)<-sort(colnames(pulp.frame))
write.csv(Spramet, "Spramet.csv")
#export to excel to rearrange for being alphabetical
SprametSort<-read.csv("Spramet-sort.csv")
rownames(SprametSort)<-SprametSort[,1]
View(SprametSort)
length(colnames(SprametSort))
SprametSort<-SprametSort[,2:969]
SprameSorted<-t(SprametSort)
condTis<-c(rep("Pulp",70), rep("Skin", 71))
obj.T.clr<-aldex.clr(reads=SprameSorted, conds=condTis, denom="zero")
obj.T.e<-aldex.effect(obj.T.clr)
obj.T.t<-aldex.ttest(obj.T.clr)
obj.T.te<-cbind(obj.T.e, obj.T.t)
obj.T.te.s=obj.T.te[which(obj.T.te$we.eBH<0.05),]
obj.T.em<-merge(obj.T.te, as.data.frame(tax_table(spraydata)), by="row.names")
obj.T.em.h<-obj.T.em[which(obj.T.em$effect>.5),]
obj.T.em.l<-obj.T.em[which(-.5>obj.T.em$effect),]
obj.T.em.whole<-rbind(obj.T.em.h,obj.T.em.l)
write.table(obj.T.em.whole, "TissueDIFFABNeffects.txt")
obj.T.em.tests<-obj.T.em.whole[-which(obj.T.em.whole$Genus == "Unassigned"),]

FigDAT<-ggplot(obj.T.em.tests, aes(x=Genus, y=effect))+
  geom_point(size=2)+  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust=0.1))  +
  geom_line(size=1)+
  geom_hline(yintercept=0)+
  scale_fill_manual()+
  labs(title="Differential Abundance Between Tissues")
tiff(filename="FigurePLACEHOLDER.tiff", width= 750, height= 500, units="px")
FigDAT
dev.off()


#obj.T<-aldex(reads=SprameSorted, conditions=condTis, test="kw")
#obj.T.sig<-obj.T[which(obj.T$kw.eBH<alpha),]
#obj.T.sig=cbind(as(obj.T.sig, "data.frame"), as(tax_table(spraydata)[rownames(obj.T.sig),], "matrix"))
#obj.T.sig.t<-obj.T.sig[,c(1:4,10 )]
#write.table(obj.T.sig.t, "Tissue-Diffential-Abundance-KruskalWallace.txt")





SpramePulp<-SprameSorted[,1:70]
condSpray.P<-c(rep("Untreated",12),rep("Abound",12), rep("AboundNuFilm",12),rep("DoubleNickel",10),rep("Conventional",12),rep("LifeGard",12))
condSpray.S<-c(rep("Untreated",12),rep("Abound",12), rep("AboundNuFilm",12),rep("DoubleNickel",12), rep("Conventional",11),rep("LifeGard",12))
SprameSkin<-SprameSorted[,71:141]

obj.P<-aldex(reads=SpramePulp, conditions = condSpray.P,test="kw")
obj.P.sig<-obj.P[which(obj.P$kw.eBH<alpha),]
obj.P.sig=cbind(as(obj.P.sig, "data.frame"), as(tax_table(pulpdata)[rownames(obj.P.sig),], "matrix"))
obj.P.sig.t<-obj.P.sig[,c(1:4,10 )]

write.table(obj.P.sig.t, "Pulp-Spray-Diff-Abundance-KW.txt")


obj.S<-aldex(reads=SprameSkin, conditions=condSpray.S, test="kw")  
obj.S.sig<-obj.S[which(obj.S$kw.eBH<alpha),]
obj.S.sig=cbind(as(obj.S.sig, "data.frame"), as(tax_table(skindata)[rownames(obj.S.sig),], "matrix"))
obj.S.sig.t<-obj.S.sig[,c(1:4,10 )]
write.table(obj.S.sig.t, "Skin-Spray-Diff-ABN-KW.txt")




#####

#Abound v Control Comparison Skin
skin.UU.Ab.Sprame<-SprameSkin[, c(1:24)]
skin.UU.Ab.clr<-aldex.clr(reads=skin.UU.Ab.Sprame, conds=c(rep("Untreated", 12), rep("Abound",12)), denom="zero")
skin.UU.Ab.eff<-aldex.effect(skin.UU.Ab.clr)
skin.UU.Ab.ttest<-aldex.ttest(skin.UU.Ab.clr)
skin.UU.Ab.sig<-skin.UU.Ab.ttest[which(skin.UU.Ab.ttest$we.eBH<0.05),]
skin.UU.Ab.hybrid<-merge(skin.UU.Ab.eff, skin.UU.Ab.sig, by="row.names")
row.names(skin.UU.Ab.hybrid)=skin.UU.Ab.hybrid$Row.names
skin.UU.Ab.hybrid=skin.UU.Ab.hybrid[,2:12]
skin.UU.Ab.t<-merge(skin.UU.Ab.hybrid, as.data.frame(tax_table(skindata)), by="row.names")
skin.UU.Ab.testing<-merge(skin.UU.Ab.ttest, skin.UU.Ab.t, by="row.names")
skin.UU.Ab.h<-skin.UU.Ab.t[which(skin.UU.Ab.t$effect>0.5),]
skin.UU.Ab.l<-skin.UU.Ab.t[which(-0.5>skin.UU.Ab.t$effect),]
skin.UU.Ab.w<-rbind(skin.UU.Ab.h, skin.UU.Ab.l)
skin.UU.Ab.tests<-skin.UU.Ab.w[-which(skin.UU.Ab.w$Genus == "Unassigned"),] #omit unassigned taxa: not informative without other info?
#Hanseniaspora more prevalent in Abound

ggplot(skin.UU.Ab.tests, aes(x=Genus, y=effect))+
  geom_point(size=1)+  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust=0.5))  +
  scale_color_continuous()+
  labs(title="Abound Treatment (Skin)")


#Abound+NuFilm v Control Comparison Skin ###THIS ONE IS THE MOST UP TO DATE CODE~!!
skin.UU.AN.Sprame<-SprameSkin[, c(1:12, 25:36)]
skin.UU.AN.clr<-aldex.clr(reads=skin.UU.AN.Sprame, conds=c(rep("Untreated", 12), rep("Abound+NuFilm",12)), denom="zero")
skin.UU.AN.eff<-aldex.effect(skin.UU.AN.clr)
skin.UU.AN.ttest<-aldex.ttest(skin.UU.AN.clr)
skin.UU.AN.sig=skin.UU.AN.ttest[which(skin.UU.AN.ttest$we.eBH<0.05),]
skin.UU.AN.hybrid<-merge(skin.UU.AN.eff, skin.UU.AN.sig, by="row.names")
row.names(skin.UU.AN.hybrid)=skin.UU.AN.hybrid$Row.names
skin.UU.AN.hybrid=skin.UU.AN.hybrid[,2:12]
skin.UU.AN.t<-merge(skin.UU.AN.hybrid, as.data.frame(tax_table(skindata)), by="row.names")
skin.UU.AN.h<-skin.UU.AN.t[which(skin.UU.AN.t$effect>0.5),]
skin.UU.AN.l<-skin.UU.AN.t[which(-0.5>skin.UU.AN.t$effect),]
skin.UU.AN.w<-rbind(skin.UU.AN.h, skin.UU.AN.l)
skin.UU.AN.tests<-skin.UU.AN.w[-which(skin.UU.AN.w$Genus == "Unassigned"),] #omit unassigned taxa: not informative without other info?
#Hanseniaspora more prevalent in Abound+NuFilm

ggplot(skin.UU.AN.tests, aes(x=Genus, y=effect))+
  geom_point(size=1)+  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust=0.5))  +
  scale_color_continuous()+
labs(title="Abound+NuFilm Treatment (Skin)")

#Abound v Abound+ NuFilm Comparison Skin
skin.Ab.AN.Sprame<-SprameSkin[, c(13:24, 25:36)]
skin.Ab.AN.clr<-aldex.clr(reads=skin.Ab.AN.Sprame, conds=c(rep("Abound", 12), rep("Abound+NuFilm",12)), denom="zero")
skin.Ab.AN.eff<-aldex.effect(skin.Ab.AN.clr)
skin.Ab.AN.ttest<-aldex.ttest(skin.Ab.AN.clr)
skin.Ab.AN.sig=skin.Ab.AN.ttest[which(skin.Ab.AN.ttest$we.eBH<0.05),]

# NOTHING SIGNIFICANT BETWEEN THEM
#skin.Ab.AN.hybrid<-merge(skin.Ab.AN.eff, skin.Ab.AN.sig, by="row.names")
#row.names(skin.Ab.AN.hybrid)=skin.Ab.AN.hybrid$Row.names
#skin.Ab.AN.hybrid=skin.Ab.AN.hybrid[,2:12]
#skin.Ab.AN.t<-merge(skin.Ab.AN.hybrid, as.data.frame(tax_table(skindata)), by="row.names")

#skin.Ab.AN.h<-skin.Ab.AN.t[which(skin.Ab.AN.t$effect>0.5),]
#skin.Ab.AN.l<-skin.Ab.AN.t[which(-0.5>skin.Ab.AN.t$effect),]
#skin.Ab.AN.w<-rbind(skin.Ab.AN.h, skin.Ab.AN.l)
#skin.Ab.AN.tests<-skin.Ab.AN.w[-which(skin.Ab.AN.w$Genus == "Unassigned"),] #omit unassigned taxa: not informative without other info?


#ggplot(skin.Ab.AN.tests, aes(x=Genus, y=effect))+
#  geom_point(size=2)+  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust=0.5))  +
 # scale_color_continuous()+
  #labs(title="Abound Comparison (Skin)")

#Double Nickel Comparison Skin
skin.UU.DN.Sprame<-SprameSkin[, c(1:12,37:48)]
skin.UU.DN.clr<-aldex.clr(reads=skin.UU.DN.Sprame, conds=c(rep("Untreated", 12), rep("Double Nickel",12)), denom="zero")
skin.UU.DN.eff<-aldex.effect(skin.UU.DN.clr)
skin.UU.DN.ttest<-aldex.ttest(skin.UU.DN.clr)
skin.UU.DN.sig=skin.UU.DN.ttest[which(skin.UU.DN.ttest$we.eBH<0.05),]
skin.UU.DN.hybrid<-merge(skin.UU.DN.eff, skin.UU.DN.sig, by="row.names")
row.names(skin.UU.DN.hybrid)=skin.UU.DN.hybrid$Row.names
skin.UU.DN.hybrid=skin.UU.DN.hybrid[,2:12]
skin.UU.DN.t<-merge(skin.UU.DN.hybrid, as.data.frame(tax_table(skindata)), by="row.names")
skin.UU.DN.h<-skin.UU.DN.t[which(skin.UU.DN.t$effect>0.5),]
skin.UU.DN.l<-skin.UU.DN.t[which(-0.5>skin.UU.DN.t$effect),]
skin.UU.DN.w<-rbind(skin.UU.DN.h, skin.UU.DN.l)
skin.UU.DN.tests<-skin.UU.DN.w[-which(skin.UU.DN.w$Genus == "Unassigned"),] #omit unassigned taxa: not informative without other info?
##HANSENIASPORA MORE PREVALENT IN DN???

ggplot(skin.UU.DN.tests, aes(x=Genus, y=effect))+
  geom_point(size=1)+  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust=0.5))  +
  scale_fill_manual()+
  labs(title="Double Nickel (Skin)")


#skin.Ab.JA.Sprame<-SprameSkin[,c(13:24, 49:59)]
#skin.Ab.JA.clr<-aldex.clr(reads=skin.Ab.JA.Sprame, conds=c(rep("Abound", 12), rep("Conventional", 11)), denom="zero")
#skin.Ab.JA.eff<-aldex.effect(skin.Ab.JA.clr)
#skin.Ab.JA.ttest<-aldex.ttest(skin.Ab.JA.clr)
#skin.Ab.JA.sig=skin.Ab.JA.ttest[which(skin.Ab.JA.ttest$we.eBH<0.05),]
#Convention Comparison Skin
skin.UU.JA.Sprame<-SprameSkin[, c(1:12,49:59)]
skin.UU.JA.clr<-aldex.clr(reads=skin.UU.JA.Sprame, conds=c(rep("Untreated", 12), rep("Conventional",11)), denom="zero")
skin.UU.JA.eff<-aldex.effect(skin.UU.JA.clr)
skin.UU.JA.ttest<-aldex.ttest(skin.UU.JA.clr)
skin.UU.JA.t<-merge(skin.UU.JA.eff, as.data.frame(tax_table(skindata)), by="row.names")
skin.UU.JA.h<-skin.UU.JA.t[which(skin.UU.JA.t$effect>0.5),]
skin.UU.JA.l<-skin.UU.JA.t[which(-0.5>skin.UU.JA.t$effect),]
skin.UU.JA.w<-rbind(skin.UU.JA.h, skin.UU.JA.l)
skin.UU.JA.tests<-skin.UU.JA.w[-which(skin.UU.JA.w$Genus == "Unassigned"),] #omit unassigned taxa: not informative without other info?
#Hanseniaspora more prevalent in Conventional 

ggplot(skin.UU.JA.tests, aes(x=Genus, y=effect))+
  geom_point(size=1)+  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust=0.5))  +
  scale_fill_manual()+
  labs(title="Conventional Treatment (Skin)")

#Lifegard Comparison Skin
skin.UU.LG.Sprame<-SprameSkin[, c(1:12,60:71)]
skin.UU.LG.clr<-aldex.clr(reads=skin.UU.LG.Sprame, conds=c(rep("Untreated", 12), rep("Lifegard",12)), denom="zero")
skin.UU.LG.eff<-aldex.effect(skin.UU.LG.clr)
skin.UU.LG.ttest<-aldex.ttest(skin.UU.LG.clr)
skin.UU.LG.t<-merge(skin.UU.LG.eff, as.data.frame(tax_table(skindata)), by="row.names")
skin.UU.LG.h<-skin.UU.LG.t[which(skin.UU.LG.t$effect>1),]
skin.UU.LG.l<-skin.UU.LG.t[which(-1>skin.UU.LG.t$effect),]
skin.UU.LG.w<-rbind(skin.UU.LG.h, skin.UU.LG.l)
skin.UU.LG.tests<-skin.UU.LG.w[-which(skin.UU.LG.w$Genus == "Unassigned"),] #omit unassigned taxa: not informative without other info?
 #Hanseniaspora interestingly not more abundant in lifegard, unlike all other treatments...

ggplot(skin.UU.LG.tests, aes(x=Genus, y=effect))+
  geom_point(size=1)+  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust=0.5))  +
  scale_fill_manual()+
  labs(title="Lifegard (Skin)")


#Abound v Control Comparison pulp
pulp.UU.Ab.Sprame<-SpramePulp[, c(1:24)]
pulp.UU.Ab.clr<-aldex.clr(reads=pulp.UU.Ab.Sprame, conds=c(rep("Untreated", 12), rep("Abound",12)), denom="zero")
pulp.UU.Ab.eff<-aldex.effect(pulp.UU.Ab.clr)
pulp.UU.Ab.t<-merge(pulp.UU.Ab.eff, as.data.frame(tax_table(pulpdata)), by="row.names")
pulp.UU.Ab.h<-pulp.UU.Ab.t[which(pulp.UU.Ab.t$effect>0.5),]
pulp.UU.Ab.l<-pulp.UU.Ab.t[which(-0.5>pulp.UU.Ab.t$effect),]
pulp.UU.Ab.w<-rbind(pulp.UU.Ab.h, pulp.UU.Ab.l)
pulp.UU.Ab.tests<-pulp.UU.Ab.w[-which(pulp.UU.Ab.w$Genus == "Unassigned"),] #omit unassigned taxa: not informative without other info?


UU.Ab.DA.pulp<-ggplot(pulp.UU.Ab.tests, aes(x=Genus, y=effect))+
  geom_point(size=2)+  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust=0.5))  +
  scale_color_continuous()+
  labs(title="Abound Treatment (pulp)")


#Abound+NuFilm v Control Comparison pulp
pulp.UU.AN.Sprame<-SpramePulp[, c(1:12, 25:36)]
pulp.UU.AN.clr<-aldex.clr(reads=pulp.UU.AN.Sprame, conds=c(rep("Untreated", 12), rep("Abound+NuFilm",12)), denom="zero")
pulp.UU.AN.eff<-aldex.effect(pulp.UU.AN.clr)
pulp.UU.AN.t<-merge(pulp.UU.AN.eff, as.data.frame(tax_table(pulpdata)), by="row.names")
pulp.UU.AN.h<-pulp.UU.AN.t[which(pulp.UU.AN.t$effect>0.5),]
pulp.UU.AN.l<-pulp.UU.AN.t[which(-0.5>pulp.UU.AN.t$effect),]
pulp.UU.AN.w<-rbind(pulp.UU.AN.h, pulp.UU.AN.l)
##pulp.UU.AN.tests<-pulp.UU.AN.w[-which(pulp.UU.AN.w$Genus == "Unassigned"),] #omit unassigned taxa: not needed here


UU.AN.DA.pulp<-ggplot(pulp.UU.AN.w, aes(x=Genus, y=effect))+
  geom_point(size=2)+  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust=0.5))  +
  scale_color_continuous()+
  labs(title="Abound+NuFilm Treatment (Pulp)")

ggarrange(UU.Ab.DA.pulp, UU.AN.DA.pulp)

#Double Nickel Comparison Pulp
Pulp.UU.DN.Sprame<-SpramePulp[, c(1:12,37:46)]
Pulp.UU.DN.clr<-aldex.clr(reads=Pulp.UU.DN.Sprame, conds=c(rep("Untreated", 12), rep("Double Nickel",10)), denom="zero")
Pulp.UU.DN.eff<-aldex.effect(Pulp.UU.DN.clr)
Pulp.UU.DN.t<-merge(Pulp.UU.DN.eff, as.data.frame(tax_table(pulpdata)), by="row.names")
Pulp.UU.DN.h<-Pulp.UU.DN.t[which(Pulp.UU.DN.t$effect>0.5),]
Pulp.UU.DN.l<-Pulp.UU.DN.t[which(-0.5>Pulp.UU.DN.t$effect),]
Pulp.UU.DN.w<-rbind(Pulp.UU.DN.h, Pulp.UU.DN.l)
Pulp.UU.DN.tests<-Pulp.UU.DN.w[-which(Pulp.UU.DN.w$Genus == "Unassigned"),] #omit unassigned taxa: not informative without other info?


UU.DN.DA.pulp<- ggplot(Pulp.UU.DN.tests, aes(x=Genus, y=effect))+
  geom_point(aes(fill=cut(effect, c(-Inf, 0, Inf))),size=1)+  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust=0.5))  +
  scale_fill_manual(name="effect",
                    values= c("(-Inf,0]"="blue",
                              "(0, Inf]"="red"))+
  labs(title="Double Nickel (Pulp)")

#Conventional Comparison Pulp
Conventional.Sprame.Pulp<-SpramePulp[,c(1:12, 47:58)]
conventional.effects.clr<-aldex.clr(reads=Conventional.Sprame.Pulp, conds= c(rep("Untreated", 12), rep("Conventional",12)), denom="zero")
cef<-aldex.effect(conventional.effects.clr)
cef=cbind(as(cef, "data.frame"), as(tax_table(pulpdata)[rownames(cef),], "matrix"))


UU.JA<-merge(cef, as.data.frame(tax_table(pulpdata)), by="row.names")
UU.JA.h<-UU.JA[which(UU.JA$effect>0.5),]
UU.JA.l<-UU.JA[which(-0.5>UU.JA$effect),]
UU.JA.w<-rbind(UU.JA.h, UU.JA.l)

UU.JA.tests<-UU.JA.w[-which(UU.JA.w$Genus == "Unassigned"),] #omit unassigned taxa: not informative without other info?



UU.JA.DA.pulp<-ggplot(UU.JA.tests, aes(x=Genus, y=effect))+
  geom_point(size=3)+  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust=0.5))  +
  scale_fill_discrete()+
  labs(title="Conventional Treatment (Pulp)")

#Lifegardl Comparison Pulp
Pulp.UU.LG.Sprame<-SpramePulp[, c(1:12,59:70)]
Pulp.UU.LG.clr<-aldex.clr(reads=Pulp.UU.LG.Sprame, conds=c(rep("Untreated", 12), rep("Lifegard",12)), denom="zero")
Pulp.UU.LG.eff<-aldex.effect(Pulp.UU.LG.clr)
Pulp.UU.LG.t<-merge(Pulp.UU.LG.eff, as.data.frame(tax_table(pulpdata)), by="row.names")
Pulp.UU.LG.h<-Pulp.UU.LG.t[which(Pulp.UU.LG.t$effect>0.5),]
Pulp.UU.LG.l<-Pulp.UU.LG.t[which(-0.5>Pulp.UU.LG.t$effect),]
Pulp.UU.LG.w<-rbind(Pulp.UU.LG.h, Pulp.UU.LG.l)
Pulp.UU.LG.tests<-Pulp.UU.LG.w[-which(Pulp.UU.LG.w$Genus == "Unassigned"),] #omit unassigned taxa: not informative without other info?


UU.LG.DA.pulp<-ggplot(Pulp.UU.LG.tests, aes(x=Genus, y=effect))+
  geom_point(aes(fill=cut(effect, c(-Inf, 0, Inf))),size=2)+  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust=0.5))  +
  scale_fill_manual(name="effect",
                    values= c("(-Inf,0]"="blue",
                              "(0, Inf]"="red"))+
  labs(title="Lifegard (Pulp)")

#Abound Comparisons Pulp
Abound.Sprame.Pulp<-SpramePulp[,c(13:36)]
abn.pulp.clr<-aldex.clr(reads=Abound.Sprame.Pulp, conds=c(rep("Abound", 12), rep("Abound+NuFilm", 12)), denom="zero")
abn.pulp.ef<-aldex.effect(abn.pulp.clr)
abn.pulp.w<-merge(abn.pulp.ef, as.data.frame(tax_table(spraydata)), by="row.names")
abn.pulp.w.high<-abn.pulp.w[which(abn.pulp.w$effect>0.5),]
abn.pulp.w.low<-abn.pulp.w[which(-0.5>abn.pulp.w$effect),]
abn.pulp.whole<-rbind(abn.pulp.w.high, abn.pulp.w.low)
#abn.pulp.tests<-abn.pulp.whole[-which(abn.pulp.whole$Genus == "Unassigned"),] #omit unassigned taxa: not needed


ggplot(abn.pulp.whole, aes(x=Genus, y=effect))+
  geom_point(size=3)+  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust=0.5))  +
  scale_fill_discrete()+
  labs(title="Abound Comparison (Pulp)")

#NOTES: This is indicating that Abound with NuFilm has LESS Botrytis, MORE Alternaria/Epicoccum/etc.

#Abound Comparisons Skin
Abound.Sprame.skin<-SprameSkin[,c(13:36)]
abn.skin.clr<-aldex.clr(reads=Abound.Sprame.skin, conds=c(rep("Abound", 12), rep("Abound+NuFilm", 12)), denom="zero")
abn.skin.ef<-aldex.effect(abn.skin.clr)
abn.skin.ttest<-aldex.ttest(abn.skin.clr)
abn.skin.w<-merge(abn.skin.ef, as.data.frame(tax_table(spraydata)), by="row.names")
abn.skin.w.high<-abn.skin.w[which(abn.skin.w$effect>.5),]
abn.skin.w.low<-abn.skin.w[which(-.5>abn.skin.w$effect),]
abn.skin.whole<-rbind(abn.skin.w.high, abn.skin.w.low)
abn.skin.tests<-abn.skin.whole[-which(abn.skin.whole$Genus == "Unassigned"),] #omit unassigned taxa: not informative without other info?

ggplot(abn.skin.tests, aes(x=Genus, y=effect))+
  geom_point(size=3)+  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust=0.5))  +
  scale_fill_discrete()+
  labs(title="Abound Comparison (Skin)")


#Skin other pairs...
Skin.UU<-SprameSkin[,c(1:12)] #12
Skin.Ab<-SprameSkin[,c(13:24)] #12
Skin.AN<-SprameSkin[,c(25:36)] #12
Skin.DN<-SprameSkin[,c(37:48)] #12
Skin.JA<-SprameSkin[,c(49:59)] #11
Skin.LG<-SprameSkin[,c(60:71)] #12 
#For each combination, you can replace everything in the form of Pulp.X.Y with X=Trt1 Y=Trt2.
#Ab=Abound, AN= Abound+NuFilmP, DN=DOuble Nickle, JA=Conventional, LG=Lifegard. Make note to alter the number in the conditions statement if doing a comparison involving Conventional

Skin.Ab.LG.clr<-aldex.clr(reads=cbind(Skin.Ab, Skin.LG), conds=c(rep("Abound", 12), rep("LG", 12)), denom="zero")
Skin.Ab.LG.ttest<-aldex.ttest(Skin.Ab.LG.clr)
Skin.Ab.LG.sig=Skin.Ab.LG.ttest[which(Skin.Ab.LG.ttest$wi.eBH<0.05),] #unbalacned data set points towards using Wilcoxon over means
Skin.Ab.LG.eff<-aldex.effect(Skin.Ab.LG.clr)
Skin.Ab.LG.hybrid<-merge(Skin.Ab.LG.eff, Skin.Ab.LG.sig, by="row.names")
row.names(Skin.Ab.LG.hybrid)=Skin.Ab.LG.hybrid$Row.names
Skin.Ab.LG.hybrid=Skin.Ab.LG.hybrid[,2:12]
Skin.Ab.LG.t<-merge(Skin.Ab.LG.hybrid, as.data.frame(tax_table(skindata)), by="row.names")
#Successful skin comparisons include:
#Abound vs Lifegard (Hanseniaspora (OTU_249)) More prevalent in Abound (CHECK IF UVARUM)
#Conventional vs Lifegard (Hanseniaspora (OTU_249)) More prevalent in Conventional (!!!!!)
#Failed skin comparisons include:
#Abound vs Abound+NuFilm
#Abound vs Double Nickel
#Abound vs Conventional
#Abound+NuFilm vs Double Nickel
#Abound+NuFilm vs Conventional
#Abound+NuFilm vs Lifegard
#Double Nickel vs Conventional
#Double Nickel vs Lifegard (OTU_249 CLOSE though)
Skin.JA.LG.clr<-aldex.clr(reads=cbind(Skin.JA, Skin.LG), conds=c(rep("Conventional", 11), rep("Lifegard", 12)), denom="zero")
Skin.JA.LG.ttest<-aldex.ttest(Skin.JA.LG.clr)
Skin.JA.LG.sig=Skin.JA.LG.ttest[which(Skin.JA.LG.ttest$wi.eBH<0.05),] #unbalacned data set points towards using Wilcoxon over means
Skin.JA.LG.eff<-aldex.effect(Skin.JA.LG.clr)
Skin.JA.LG.hybrid<-merge(Skin.JA.LG.eff, Skin.JA.LG.sig, by="row.names")
row.names(Skin.JA.LG.hybrid)=Skin.JA.LG.hybrid$Row.names
Skin.JA.LG.hybrid=Skin.JA.LG.hybrid[,2:12]
Skin.JA.LG.t<-merge(Skin.JA.LG.hybrid, as.data.frame(tax_table(skindata)), by="row.names")





#Pulp other pairs...
Pulp.UU<-SpramePulp[,c(1:12)] #12
Pulp.Ab<-SpramePulp[,c(13:24)] #12
Pulp.AN<-SpramePulp[,c(25:36)] #12
Pulp.DN<-SpramePulp[,c(37:46)] #10
Pulp.JA<-SpramePulp[,c(47:58)] #12
Pulp.LG<-SpramePulp[,c(59:70)] #12 
#For each combination, you can replace everything in the form of Pulp.X.Y with X=Trt1 Y=Trt2.
#Ab=Abound, AN= Abound+NuFilmP, DN=DOuble Nickle, JA=Conventional, LG=Lifegard. Make note to alter the number in the conditions statement if doing a comparison involving double nickel

#Aboudn comparison redo with this version
Pulp.Ab.AN.clr<-aldex.clr(reads=cbind(Pulp.Ab, Pulp.AN), conds=c(rep("Abound", 12), rep("Abound+", 12)), denom="zero")
Pulp.Ab.AN.ttest<-aldex.ttest(Pulp.Ab.AN.clr)
Pulp.Ab.AN.sig=Pulp.Ab.AN.ttest[which(Pulp.Ab.AN.ttest$wi.eBH<0.05),] #unbalacned data set points towards using Wilcoxon over means
Pulp.Ab.AN.eff<-aldex.effect(Pulp.Ab.AN.clr)
Pulp.Ab.AN.hybrid<-merge(Pulp.Ab.AN.eff, Pulp.Ab.AN.sig, by="row.names")
row.names(Pulp.Ab.AN.hybrid)=Pulp.Ab.AN.hybrid$Row.names
Pulp.Ab.AN.hybrid=Pulp.Ab.AN.hybrid[,2:12]
Pulp.Ab.AN.t<-merge(Pulp.Ab.AN.hybrid, as.data.frame(tax_table(pulpdata)), by="row.names")
#The above is an example where there were differenes. Differences were also found for:
#Abound vs Lifegard (Botrytis (6, 746), Epicoccum dendrobii (2,288), Alternaria sp (17, 758, 506))
#Lifegard vs Double Nickel (OTU_6, Botrytis)
#Abound+NuFilm vs Conventional (OTU_6 and 746, both Botrytis)
#Lifegard vs Conventional (OTU_6 & 746, Botrytis, OTU_288, Epicoccum dendrobii, OTU_506 & 768, Alternaria sp)

###Potential antagonism between Epicoccum dendrobii/Alternaria sp and Botrytis


#Pulp Aboudn v Double Nickel (NO DIFFERENCES)
#Pulp.Ab.DN.clr<-aldex.clr(reads=cbind(Pulp.Ab, Pulp.DN), conds=c(rep("Abound", 12), rep("Double Nickel", 10)), denom="zero")
#Pulp.Ab.DN.ttest<-aldex.ttest(Pulp.Ab.DN.clr)
#Pulp.Ab.DN.sig=Pulp.Ab.DN.ttest[which(Pulp.Ab.DN.ttest$wi.eBH<0.05),]
#Pulp.Ab.DN.eff<-aldex.effect(Pulp.Ab.DN.clr)
#Pulp.Ab.DN.hybrid<-merge(Pulp.Ab.DN.eff, Pulp.Ab.DN.sig, by="row.names")
#row.names(Pulp.Ab.DN.hybrid)=Pulp.Ab.DN.hybrid$Row.names
#Pulp.Ab.DN.hybrid=Pulp.Ab.DN.hybrid[,2:12]
#Pulp.Ab.DN.t<-merge(Pulp.Ab.DN.hybrid, as.data.frame(tax_table(pulpdata)), by="row.names")

#The above example is one that had no difference. Others with no different taxa
#Abound vs Conventional
#Abound vs Double Nickel
#Abound+NuFilm vs Double Nickel
#Abound+NuFilm vs Lifegard
#Conventional vs Double Nickel
Pulp.AN.LG.clr<-aldex.clr(reads=cbind(Pulp.AN, Pulp.LG), conds=c(rep("Abound+NuFilm", 12), rep("Lifegard", 12)), denom="zero")
Pulp.AN.LG.ttest<-aldex.ttest(Pulp.AN.LG.clr)
Pulp.AN.LG.sig=Pulp.AN.LG.ttest[which(Pulp.AN.LG.ttest$wi.eBH<0.05),]
Pulp.AN.LG.eff<-aldex.effect(Pulp.AN.LG.clr)
Pulp.AN.LG.hybrid<-merge(Pulp.AN.LG.eff, Pulp.AN.LG.sig, by="row.names")
row.names(Pulp.AN.LG.hybrid)=Pulp.AN.LG.hybrid$Row.names
Pulp.AN.LG.hybrid=Pulp.AN.LG.hybrid[,2:12]
Pulp.AN.LG.t<-merge(Pulp.AN.LG.hybrid, as.data.frame(tax_table(pulpdata)), by="row.names")
