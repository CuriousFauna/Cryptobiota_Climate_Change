#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
#  Community Structure - Phyla and Families 
#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

library(dplyr)
library(plyr)
library(reshape2)
library(RVAideMemoire) 
library(ggplot2)


#><><><><><><><><><><
#  Figure 3A
#><><><><><><><><><><

PhyAnnotate<-read.csv("Annotations_Working.csv")
PhyAnnotate1<-as.data.table(PhyAnnotate[,c("ID","FinalPhylum")])
FamAnnotate1<-as.data.table(PhyAnnotate[,c("ID","FinalFamily")])

data<-read.csv("ARMS_Data_Working.csv")
meta<-read.csv("MetaData_ARMS.csv")

data2<-merge(meta,data, by = "Sample")

TreatProp<-as.data.table(data2)
TreatProp2<-TreatProp[,-c("Groups", "HeaderTank","Temperature","CO2","ARMS_Unit","Sample")]
TreatProp3<-melt(TreatProp2, id = c("Treatment"), variable = "ID")
TreatProp4<-merge(PhyAnnotate1,TreatProp3, by = "ID" )
TreatProp5<-ddply(TreatProp4, .(FinalPhylum, Treatment), summarize, Sum=sum(value))
TreatProp5$Phylum<-ifelse (TreatProp5$FinalPhylum == "CHORDATA" | 
                             TreatProp5$FinalPhylum == "PLATYHELMINTHES" | 
                             TreatProp5$FinalPhylum == "NEMERTEA" | 
                             TreatProp5$FinalPhylum == "BRYOZOA" | 
                             TreatProp5$FinalPhylum == "OCHROPHYTA" | 
                             TreatProp5$FinalPhylum == "GASTROTRICHA", 
                           "Other", TreatProp5$FinalPhylum )
TreatProp6<-ddply(TreatProp5, .(Treatment, Phylum), summarize, Sum = sum(Sum))
seqTotsTreat<-aggregate(Sum ~ Treatment, data = TreatProp6, FUN = sum)
colnames(seqTotsTreat)[2]<-"Total"
TreatProp7<-merge(TreatProp6, seqTotsTreat, by = "Treatment")
TreatProp7$RelAbun<-TreatProp7$Sum/TreatProp7$Total
#Other was null so removing it from graphs
TreatProp7<-subset(TreatProp7, Phylum != "Other")
PhyOR<-c("ANNELIDA","ARTHROPODA","CNIDARIA","ECHINODERMATA","MOLLUSCA",
         "PORIFERA","RHODOPHYTA","Other")
TreatProp7$Phylum<-factor(TreatProp7$Phylum, levels = PhyOR) 

ggplot(TreatProp7,aes(x=Treatment, y = RelAbun, fill=Phylum))+
  geom_bar(stat='identity', color = "black", size = .25)+
  scale_fill_manual(values = c( "ANNELIDA" = "#98D4C0", "ARTHROPODA" = "#9E6A28", "CNIDARIA" = "#BCB450",
                                "ECHINODERMATA" ="#A882BB","MOLLUSCA" = "#F9BFCB","PORIFERA" ="#C4E8F5",
                                "RHODOPHYTA" = "#F3EEBA", "Other" = "snow2"))+
  ylab("Proportional of Reads")+
  scale_x_discrete(drop = FALSE)+
  xlab("")+
  scale_y_continuous(expand = c(0, 0), limits = c(0,1) )+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.text=element_text(size=12),
        panel.border = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size=14),
        axis.line.x = element_line(color="black", size = .5),
        axis.line.y = element_line(color="black", size = .5))


#><><><><><><><><><><
#  Figure 3B
#><><><><><><><><><><

### Proportion By Phylum
PropPhy<-data2
#Randomly subsample of 5 units per treatment to account for uneven sampling
set.seed(42)
PropPhy2 <- PropPhy %>% 
  group_by(Treatment) %>% 
  sample_n(5) 
PropPhy3<-as.data.table(PropPhy2)
PropPhy4<-PropPhy3[,-c("Groups", "HeaderTank","Temperature","CO2","ARMS_Unit","Sample")]
PropPhy5<-melt(PropPhy4, id = c("Treatment"), variable = "ID", value.name = "Sum")
PropPhy6<-merge(PhyAnnotate1,PropPhy5, by = "ID" )
PropPhy7<-ddply(PropPhy6, .(FinalPhylum, Treatment), summarize, Sum=sum(Sum))
seqTots<-aggregate(Sum ~ FinalPhylum, data = PropPhy7, FUN = sum)
colnames(seqTots)[2]<-"Total"
PropPhy8<-merge(PropPhy7,seqTots, by = "FinalPhylum")
PropPhy9<-subset(PropPhy8, Sum > 0)
PropPhy9$RelAbun<-(PropPhy9$Sum/PropPhy9$Total)*100
PropPhy10<-PropPhy9 %>% filter(FinalPhylum %in% c("PORIFERA", "MOLLUSCA","ARTHROPODA","ANNELIDA","RHODOPHYTA","ECHINODERMATA", "CNIDARIA") )

Control<-subset(PropPhy10, Treatment == "Control")
namControl<-c("FinalPhylum","Control","SumControl","TotalPhySum","Rel")
names(Control)<-namControl

PropPhy11<-merge(Control, PropPhy10, by = "FinalPhylum") 
PropPhy11$Proportion<-PropPhy11$Sum/PropPhy11$SumControl

        #><><><><><><><><><><<><><><><><><<><
        #  DataFrame PropPhy 11 = Table S11
        #><><><><><><><><><><<><><><><><><<><

TreatOrder<-c("Control", "Acidified","Heated","AcidifiedHeated")
PropPhy11$Treatment<-factor(PropPhy11$Treatment, levels = TreatOrder)

ggplot(PropPhy11, aes(x=Treatment, y = Proportion, group = FinalPhylum))+
  geom_hline(yintercept=1, linetype="longdash",size = 1,color = "gray40")+
  geom_line(aes(color = FinalPhylum, fill = FinalPhylum), size = .75) +
  scale_colour_manual(values = c( "ANNELIDA" = "#98D4C0", "ARTHROPODA" = "#9E6A28", "CNIDARIA" = "#BCB450",
                                  "ECHINODERMATA" ="#A882BB","MOLLUSCA" = "#F9BFCB","PORIFERA" ="#C4E8F5",
                                  "RHODOPHYTA" = "#F3EEBA")) +
  geom_point(aes(color = FinalPhylum, fill = FinalPhylum), pch = 21, size = 8, col = "Black", lwd=3)+
  scale_fill_manual(values = c( "ANNELIDA" = "#98D4C0", "ARTHROPODA" = "#9E6A28", "CNIDARIA" = "#BCB450",
                                "ECHINODERMATA" ="#A882BB","MOLLUSCA" = "#F9BFCB","PORIFERA" ="#C4E8F5",
                                "RHODOPHYTA" = "#F3EEBA")) +
  ylab("Proportion of Reads by Phylum relative to Control Treatment")+
  xlab("")+
  ylim(0,4)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.text=element_text(size=12),
        axis.title.y = element_text(size = 14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1))


#><><><><><><><><><><
#  Figure 3C
#><><><><><><><><><><

# Treatment Bar Graphs 
# - starting off with random ARMS unit dataframe from Figure 3B
PropPhyBar2<-PropPhy3[,-c("Groups", "HeaderTank","Temperature","CO2","ARMS_Unit")]
PropPhyBar3<-melt(PropPhyBar2, id = c("Treatment","Sample"), variable = "ID")
PropPhyBar4<-merge(PhyAnnotate,PropPhyBar3, by = "ID" )
PropPhyBar5<-ddply(PropPhyBar4, .(FinalPhylum, Treatment, Sample), summarize, Sum=sum(value))
# Had to seperate out Rhodophyta because there was an Acid Sample without any but could not loose that sample for the test given the anova
PropPhyBar5a<-subset(PropPhyBar5, Treatment == "Acidified" & FinalPhylum == "RHODOPHYTA")
PropPhyBar5b<-subset(PropPhyBar5, Treatment != "Acidified" | FinalPhylum != "RHODOPHYTA")
PropPhyBar6<-subset(PropPhyBar5b, Sum > 0)
PropPhyBar6a<-rbind(PropPhyBar6, PropPhyBar5a)
ReadPhy<-aggregate(Sum ~ FinalPhylum, data = PropPhyBar6a, FUN = sum) # total Phyla sequences used in Table S3
colnames(ReadPhy)[2]<-"PhylaReads"
### Need to either Sum by Phyla or Sum by Treatment 
PropPhyBar7<-merge(PropPhyBar6a,ReadPhy, by = "FinalPhylum")
PropPhyBar7$RelAbun<-(PropPhyBar7$Sum/PropPhyBar7$PhylaReads)*100
PropPhyBar7$Treatment<-NULL
PropPhyBar8<-merge(meta, PropPhyBar7, by = "Sample")
PropPhyBar9<-PropPhyBar8 %>% filter(FinalPhylum %in% c("PORIFERA", "MOLLUSCA","ARTHROPODA","ANNELIDA","RHODOPHYTA","ECHINODERMATA", "CNIDARIA") )


PropPhyBar9$Treatment<-factor(PropPhyBar9$Treatment, levels = TreatOrder)
Treat_RelAbun <- ddply(PropPhyBar9, c("Treatment","FinalPhylum"), summarise,
                       N    = length(RelAbun),
                       mean = mean(RelAbun),
                       sd   = sd(RelAbun),
                       se   = sd / sqrt(N))

Treat_RelAbun$Treatment<-factor(Treat_RelAbun$Treatment, levels = TreatOrder)

cols <- c("Control" = "#56B4E9", "Acidified" = "#F0E442",
          "Heated" = "#E69F00","AcidifiedHeated" = "#999999")

TreatFacet<-ggplot(PropPhyBar9, aes(x=factor(Treatment), y=RelAbun)) + 
  geom_boxplot(position="identity", color = "black", width = 0.9) +
  xlab("")+
  ylab("Mean Relative Abundance - Reads")+
  facet_grid(.~FinalPhylum, scale = "free", space = "free")+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 25))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size = 14)) 

TreatFacet + aes(fill = Treatment)+
  scale_fill_manual(values = cols)

# Took pdf into Illustrator to create an more legible Figure 3c 


#><><><><><><><><><><<><><><><<><><><><><><><><<><><><><<<><><><><
#  Permutational ANOVAs on Phyla Relative Abundance - Table S13
#><><><><><><><><><><<><><><><<><><><><><><><><<><><><><<<><><><><

PropPhyBar9$Temperature<-as.factor(PropPhyBar9$Temperature)
PropPhyBar9$CO2<-as.factor(PropPhyBar9$CO2)
PropPhyBar9$Treatment<-as.factor(PropPhyBar9$Treatment)

## Annelida
ReadAnn<-subset(PropPhyBar9, FinalPhylum == "ANNELIDA")
perm.anova(RelAbun ~ Temperature * CO2, nperm=9999, data = ReadAnn)
perm.bartlett.test(RelAbun ~ Temperature * CO2, nperm=999, data = ReadAnn)
## Arthropoda
ReadArt<-subset(PropPhyBar9, FinalPhylum == "ARTHROPODA")
perm.anova(RelAbun ~ Temperature * CO2, nperm=9999, data = ReadArt)
perm.bartlett.test(RelAbun ~ Temperature * CO2, nperm=999, data = ReadArt)
## Cnidaria
ReadCni<-subset(PropPhyBar9, FinalPhylum == "CNIDARIA")
perm.anova(RelAbun ~ Temperature * CO2, nperm=999, data = ReadCni)
perm.bartlett.test(RelAbun ~ Temperature * CO2, nperm=999, data = ReadCni)
# Echinodermata
ReadEch<-subset(PropPhyBar9, FinalPhylum == "ECHINODERMATA")
perm.anova(RelAbun ~ Temperature * CO2, nperm=999, data = ReadEch)
perm.bartlett.test(RelAbun ~ Temperature * CO2, nperm=999, data = ReadEch)
# Mollusca
ReadMol<-subset(PropPhyBar9, FinalPhylum == "MOLLUSCA")
perm.anova(RelAbun ~ CO2*Temperature, nperm=999, data = ReadMol)
perm.bartlett.test(RelAbun ~ Temperature * CO2, nperm=999, data = ReadMol)
# Porifera
ReadPor<-subset(PropPhyBar9, FinalPhylum == "PORIFERA")
perm.anova(RelAbun ~ Temperature * CO2, nperm=999, data = ReadPor)
anova(aovp(RelAbun~Temperature*CO2, data =ReadPor, perm=""))
# Rhodophyta
ReadRho<-subset(PropPhyBar9, FinalPhylum == "RHODOPHYTA")
perm.anova(RelAbun ~ Temperature * CO2, nperm=999, data = ReadRho)
perm.bartlett.test(RelAbun ~ Temperature * CO2, nperm=999, data = ReadRho)


#><><><><><><><><><><
#  Figure 3D
#><><><><><><><><><><

Family <- PropPhy3 # this is the random subset of units from above
Family1<-as.data.table(Family)
Family2<-Family1[,-c("Groups", "HeaderTank","Temperature","CO2","ARMS_Unit","Sample")]
Family3<-melt(Family2, id = c("Treatment"), variable = "ID")
Family4<-merge(FamAnnotate1,Family3, by = "ID" )
Family4$PA<-ifelse(Family4$value > 0, 1, 0)
Family5<-dcast(Family4, ID + FinalFamily~ Treatment, value.var = "value", fun.aggregate = sum)
Family6<-ddply(Family5, .(FinalFamily), summarize,
               Control=sum(Control), Acidified=sum(Acidified), Heated=sum(Heated), AcidifiedHeated=sum(AcidifiedHeated))
Family6$Sum<-rowSums(Family6[,2:ncol(Family6)])
Family6$FinalFamily<-ifelse(Family6$FinalFamily == "", "Unk",Family6$FinalFamily)
Family6$RelCont<-round(Family6$Control/Family6$Sum, digits = 2)
Family6$RelAcid<-round(Family6$Acidified/Family6$Sum, digits = 2)
Family6$RelHeat<-round(Family6$Heated/Family6$Sum, digits = 2)
Family6$RelAcidHeat<-round(Family6$AcidifiedHeated/Family6$Sum, digits = 3)
Family6$MOTURelAbun<-round(Family6$Sum/sum(Family6$Sum), digits = 5)

# Getting known Families
Family7<-subset(Family6, FinalFamily != "Unk")

# Figuring out ~4-5% of Family annotated Seqeucnes
sum(Family7$Sum)*0.045 # ~10000
# Total Seq
Family8<-subset(Family7, Sum > 10000)
Family8<-Family8[,1:5]
Family8 <- droplevels(Family8[!Family8$FinalFamily == 'Plakinidae',])

Family9<-melt(Family8, id=c("FinalFamily"),variable = "Treatment" )

Amb<-subset(Family9, Treatment == "Control")
Amb2<-Amb[,c("FinalFamily","value")]
colnames(Amb2)[2]<-"ControlTotal"
Amb2<-subset(Amb2, ControlTotal != 0)

Family10<-merge(Family9, Amb2, by = c("FinalFamily"))

Family10$Prop<-Family10$value/Family10$ControlTotal
Family11<-subset(Family10, Treatment != "Control" & ControlTotal != 0)
# % of reads from these families
sum(Family10$value)/sum(Family1[,9:ncol(Family1)]) #53% of sequences

FamilyOrder<-c("Amphinomidae", "Cirratulidae","Gammaridae","Amphiuridae",
               "Ophiactidae","Hipponicidae","Vermetidae","Suberitidae")

Family11$FinalFamily<-factor(Family11$FinalFamily, levels = FamilyOrder)

ggplot(Family11, aes(x = Prop, y = FinalFamily, shape = Treatment, fill = Treatment)) +
  # geom_point(aes(shape=Treatment, colour = Treatment, fill = Treatment), size = 5)+  
  geom_point(color = 'black', size = 5) +
  scale_shape_manual(values = c(21, 24, 22))+
  theme_bw()+
  ylab("Order")+
  xlim(0,7.5)+
  geom_vline(xintercept = 1, linetype = "dashed")+
  scale_fill_manual(values = c("Acidified" = "#F0E442",
                               "Heated" = "#E69F00","AcidifiedHeated" = "#999999"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(legend.position="none")

#><><><><><><><><><><<><><><><<><><><><><><><><<><><><><<<><><><><
#  Permutational ANOVAs on Family Relative Abundance - Table S14
#><><><><><><><><><><<><><><><<><><><><><><><><<><><><><<<><><><><

FamStat<-Family1[,-c("Groups", "HeaderTank","Temperature","CO2","ARMS_Unit")]
FamStat2<-melt(FamStat, id = c("Treatment","Sample"), variable = "ID")
FamStat3<-merge(FamAnnotate1,FamStat2, by = "ID" )
FamStat4<- ddply(FamStat3, .(FinalFamily, Treatment, Sample),summarize, Sum=sum(value))
FamStat5<-FamStat4 %>% filter(FinalFamily %in% c("Amphinomidae","Amphiuridae","Cirratulidae","Gammaridae","Hipponicidae",
                                                 "Ophiactidae","Suberitidae","Vermetidae" ))
seqFam<-aggregate(Sum ~ FinalFamily, data = FamStat5, FUN = sum)
colnames(seqFam)[2]<-"Total"
FamStat6<-merge(FamStat5,seqFam, by = "FinalFamily")
FamStat6$RelAbun<-(FamStat6$Sum/FamStat6$Total)*100
FamStat6$Treatment<-NULL
FamStat7<-merge(FamStat6, meta, by = "Sample")

## Amphimonidae
Amph<-subset(FamStat7, FinalFamily == "Amphinomidae")
perm.anova(RelAbun~Temperature*CO2, nperm = 999, data = Amph)
perm.bartlett.test(RelAbun ~ Temperature * CO2, nperm=999, data = Amph)
## Amphiuridae
AmphIUR<-subset(FamStat7, FinalFamily == "Amphiuridae")
perm.anova(RelAbun~Temperature*CO2, nperm = 999, data = AmphIUR)
perm.bartlett.test(RelAbun ~ Temperature * CO2, nperm=999, data = AmphIUR)
## Cirratulidae
Cirr<-subset(FamStat7, FinalFamily == "Cirratulidae")
perm.anova(RelAbun~Temperature*CO2, nperm = 999, data = Cirr)
anova(aovp(RelAbun~Temperature*CO2, data =Cirr, perm=""))
perm.bartlett.test(RelAbun ~ Temperature * CO2, nperm=999, data = Cirr)
## Gammaridae
Gamma<-subset(FamStat7, FinalFamily == "Gammaridae")
perm.anova(RelAbun~Temperature*CO2, nperm = 999, data = Gamma)
perm.bartlett.test(RelAbun ~ Temperature * CO2, nperm=999, data = Gamma)
anova(aovp(RelAbun~Temperature*CO2, data =Gamma, perm=""))
## Hipponicidae
Hipp<-subset(FamStat7, FinalFamily == "Hipponicidae")
perm.anova(RelAbun~Temperature*CO2, nperm = 999, data = Hipp)
perm.bartlett.test(RelAbun ~ Temperature * CO2, nperm=999, data = Hipp)
## Ophiactidae
Ophi<-subset(FamStat7, FinalFamily == "Ophiactidae")
perm.anova(RelAbun~Temperature*CO2, nperm = 999, data = Ophi)
perm.bartlett.test(RelAbun ~ Temperature * CO2, nperm=999, data = Ophi)
## Suberitidae
Subb<-subset(FamStat7, FinalFamily == "Suberitidae")
perm.anova(RelAbun~Temperature*CO2, nperm = 999, data = Subb)
perm.bartlett.test(RelAbun ~ Temperature * CO2, nperm=999, data = Subb)
## Vermetidae
Verm<-subset(FamStat7, FinalFamily == "Vermetidae")
perm.anova(RelAbun~Temperature*CO2, nperm = 999, data = Verm)
perm.bartlett.test(RelAbun ~ Temperature * CO2, nperm=999, data = Verm)


#><><><><><><><><><><<><><><><<><><><><><><><><<><><><><
#  Figure S6 - MOTU data
#><><><><><><><><><><<><><><><<><><><><><><><><<><><><><

TreatMotu<-as.data.table(data2)
TreatMotu2<-TreatMotu[,-c("Groups", "HeaderTank","Temperature","CO2","ARMS_Unit","Sample")]
TreatMotu3<-melt(TreatMotu2, id = c("Treatment"), variable = "ID")
TreatMotu4<-merge(PhyAnnotate1,TreatMotu3, by = "ID" )
TreatMotu4a<-ddply(TreatMotu4, .(ID, Treatment, FinalPhylum), summarize, Sum=sum(value))
TreatMotu4b<-subset(TreatMotu4a, Sum > 0)
TreatMotu4b$PA<-rep(1,nrow(TreatMotu4b))
TreatMotu5<-ddply(TreatMotu4b, .(FinalPhylum, Treatment), summarize, Sum=sum(PA))
TreatMotu5<-as.data.table(TreatMotu5)
TreatMotu5$Phylum<-TreatMotu5$FinalPhylum
# Making these as Other cause MOTUs so few
TreatMotu5$Phylum<-ifelse(TreatMotu5$FinalPhylum == "CHORDATA" | 
                             TreatMotu5$FinalPhylum == "PLATYHELMINTHES" | 
                             TreatMotu5$FinalPhylum == "NEMERTEA" | 
                             TreatMotu5$FinalPhylum == "BRYOZOA" | 
                             TreatMotu5$FinalPhylum == "OCHROPHYTA" | 
                             TreatMotu5$FinalPhylum == "GASTROTRICHA", 
                           "Other", TreatMotu5$Phylum )

TreatMotu6<-ddply(TreatMotu5, .(Treatment, Phylum), summarize, Sum = sum(Sum))
TreatAll<-ddply(TreatMotu4, .(ID, Treatment), summarize, Sum=sum(value))
TreatAll2<-subset(TreatAll, Sum > 0)
TreatAll2$PA<-rep(1,nrow(TreatAll2))
TreatAll3<-aggregate(PA~Treatment, data=TreatAll2, FUN=sum)
TreatMotu7<-merge(TreatMotu6, TreatAll3, by = "Treatment")
TreatMotu7$RelAbun<-TreatMotu7$Sum/TreatMotu7$PA

PhyOR<-c("ANNELIDA","ARTHROPODA","CNIDARIA","ECHINODERMATA","MOLLUSCA",
         "PORIFERA","RHODOPHYTA","Other")
TreatMotu7$Phylum<-factor(TreatMotu7$Phylum, levels = PhyOR) 

ggplot(TreatMotu7,aes(x=Treatment, y = RelAbun, fill=Phylum))+
  geom_bar(stat='identity', color = "black", size = .25)+
  scale_fill_manual(values = c( "ANNELIDA" = "#98D4C0", "ARTHROPODA" = "#9E6A28", "CNIDARIA" = "#BCB450",
                                "ECHINODERMATA" ="#A882BB","MOLLUSCA" = "#F9BFCB","PORIFERA" ="#C4E8F5",
                                "RHODOPHYTA" = "#F3EEBA", "Other" = "snow1")) +
  ylab("Proportion of MOTUs")+
  xlab("")+
  scale_x_discrete(drop = FALSE)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text=element_text(size=12),
        axis.title.y = element_text(size = 14))


#><><><><><><><><><><<><>
#  Figure S6b - MOTU data
#><><><><><><><><><><<><>

#Taking the random ARMS subset dataframe
MotuPhy1<-as.data.table(PropPhy3)
MotuPhy2<-MotuPhy1[,-c("Groups", "HeaderTank","Temperature","CO2","ARMS_Unit","Sample")]
MotuPhy3<-melt(MotuPhy2, id = c("Treatment"), variable = "ID")
MotuPhy4<-merge(PhyAnnotate1,MotuPhy3, by = "ID" )
MotuPhy4$PA<-ifelse(MotuPhy4$value > 0, 1, 0)
MotuPhy5<-dcast(MotuPhy4, ID ~ Treatment, value.var = "PA", fun.aggregate = sum)
# Getting total MOTU per treatment
RelAmbientMOTUTotal<-specnumber(MotuPhy5[,2:ncol(MotuPhy5)], MARGIN = 2)
RelAmbientMOTUTotal
MotuPhy6<-as.data.frame(ifelse(MotuPhy5[2:5] > 0, 1, 0))
MotuPhy6$ID<-MotuPhy5$ID
MotuPhy7<-merge(PhyAnnotate1, MotuPhy6, by = "ID")
MotuPhy8<-melt(MotuPhy7,id=c("ID", "FinalPhylum"), variable = "Treatment")
MotuPhy9<-ddply(MotuPhy8, .(Treatment, FinalPhylum), summarize, Sum=sum(value))

PhyMOTUCount<-count(PhyAnnotate1$FinalPhylum)
colnames(PhyMOTUCount)[1]<-"FinalPhylum"

MotuPhy10<-merge(MotuPhy9,PhyMOTUCount, by = "FinalPhylum")
MotuPhy10$RelAbun<-(MotuPhy10$Sum/MotuPhy10$freq)*100
MotuPhy11<-MotuPhy10 %>% filter(FinalPhylum %in% c("PORIFERA", "MOLLUSCA","ARTHROPODA","ANNELIDA","RHODOPHYTA","ECHINODERMATA", "CNIDARIA") )

ControlMotu<-subset(MotuPhy11, Treatment == "Control")
namContMotu<-c("FinalPhylum","Control","SumControl","TotalPhySum","RelControl")
names(ControlMotu)<-namContMotu

MotuPhy12<-merge(ControlMotu, MotuPhy11, by = "FinalPhylum")
MotuPhy12$MotuProportion<-MotuPhy12$Sum/MotuPhy12$SumControl

TreatOrder<-c("Control", "Acidified","Heated","AcidifiedHeated")
MotuPhy12$Treatment<-factor(MotuPhy12$Treatment, levels = TreatOrder) 


            #><><><><><><><><><><<><><><><><><<><
            #  DataFrame MotuPhy12 11 = Table S12
            #><><><><><><><><><><<><><><><><><<><

ggplot(MotuPhy12, aes(x=Treatment, y = MotuProportion, group = FinalPhylum, color=FinalPhylum))+
  geom_hline(yintercept=1, linetype="longdash",size = 1,color = "gray40")+
  geom_line(aes(color = FinalPhylum, fill = FinalPhylum), size = .75) +
  scale_colour_manual(values = c( "ANNELIDA" = "#98D4C0", "ARTHROPODA" = "#9E6A28", "CNIDARIA" = "#BCB450",
                                  "ECHINODERMATA" ="#A882BB","MOLLUSCA" = "#F9BFCB","PORIFERA" ="#C4E8F5",
                                  "RHODOPHYTA" = "#F3EEBA")) +
  geom_point(aes(color = FinalPhylum, fill = FinalPhylum), pch = 21, size = 8, col = "Black", lwd=3)+
  scale_fill_manual(values = c( "ANNELIDA" = "#98D4C0", "ARTHROPODA" = "#9E6A28", "CNIDARIA" = "#BCB450",
                                "ECHINODERMATA" ="#A882BB","MOLLUSCA" = "#F9BFCB","PORIFERA" ="#C4E8F5",
                                "RHODOPHYTA" = "#F3EEBA")) +
  ylab("Proportion of MOTUs realtive to Control treatment")+
  xlab("")+
  ylim(0,2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text=element_text(size=12))+
  theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.title.y = element_text(size = 14))

