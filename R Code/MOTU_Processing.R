
#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
#  OTU Table Processing 
#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

# This script takes the OTU table output from JAMP processing - https://github.com/VascoElbrecht/JAMP
# for further processing.

library(data.table) # read and manipulate large data
library(dplyr) # data manipulation
library(vegan) # diversity estimates
library(seqRFLP) # create fasta file
library(EcolUtils) # Permutational Rarefraction

MOTU<-fread("3_Raw_MOTU_Table.csv")

# creating dataframe of ID tied to sequence for MASCSE step 
SequenceCheck<-MOTU[,c("ID","sequ")]

#removing columns not needed
MOTU$sort<-NULL
MOTU$sequ<-NULL

# Abundance filtering to reduce the number of false positives due to PCR and sequencing errors 
MOTU_samples<-as.data.frame(MOTU[,2:ncol(MOTU)])
sum(MOTU_samples)
rownames(MOTU_samples)<-MOTU$ID
total_MOTU<-colSums(MOTU_samples)
MOTU_rel <- as.data.frame(MOTU_samples)
for (i in 1:ncol(MOTU_samples))  MOTU_rel[,i] <- MOTU_samples[,i]/total_MOTU[i] 
MOTU_rel_2<-as.data.frame(lapply(MOTU_rel, function(x){replace(x,x <= 0.0001,0)}))
MOTU_rel_2$ID<-MOTU$ID
MOTU_rel_3 <- MOTU_rel_2 %>% dplyr::select(ID, everything())
MOTU_rel_3$MOTUReads<-rowSums(MOTU_rel_3[,2:ncol(MOTU_rel_3)])
MOTU_rel_4<-subset(MOTU_rel_3, MOTUReads != 0)
MOTU_rel_4$MOTUReads<-NULL

# Converting dataframe back to the raw counts of the OTUS that were not removed
MOTU2<-MOTU_rel_4[,2:ncol(MOTU_rel_4)]
MOTU3<-MOTU2
for (i in 1:ncol(MOTU2))  MOTU3[,i] <- MOTU2[,i]*total_MOTU[i] 
MOTU3$ID<-MOTU_rel_4$ID
postRelSum<-MOTU3
postRelSum$ID<-NULL
sum(postRelSum)
MOTU4 <- MOTU3 %>% dplyr::select(ID, everything())

#><><><><><><><><><>
#  Annotation File 
#><><><><><><><><><>
      
annotate<-fread("Tota_MOTU_Annotatations.csv") 
annotate1<-annotate[,c("ID","Reads", "FinalKingdom","FinalPhylum","FinalClass","FinalOrder","FinalFamily","FinalGenus", "FinalSpecies","Calcify", "TaxaFrequency")]

# merging the annotated file with the sequence file 
annotate2<-merge(annotate1,MOTU4,by = "ID") 
# Taking only Metazoans and MacroAlgae
annotate3<-subset(annotate2,FinalKingdom == "Metazoa" | FinalKingdom == "Plantae")
# Removing MOTUs Unclassified to Phylum
annotate4<-subset(annotate3, FinalPhylum != "Unclassified") # 279 remaining - 63%

# Obtaining sequences for IDs to run through MASCE which looks for pseudogenes
MASCE<-as.data.frame(annotate4[,c("ID")])
# Merging dataframe back to SequenceCheck inorder to acquire the actual sequences for each ID
MACSE2<-merge(MASCE, SequenceCheck, by = "ID")
# creating fasta file to run through MACSE
dataframe2fas(MACSE2, file = "OA_MACSE.fasta")

# bringing back MACSE output based on ID; these are the MOTUs not identified as pseudogenes
outputMACSE<-read.table("MACSE_output.txt", col.names = T)
colnames(outputMACSE)[1]<-"ID"

# Checking the number of identified pseudogenes 
Questionable<-anti_join(annotate4,outputMACSE, by = "ID") # there were 4

# Removing those IDs identified as pseudogenes from dataset
annotate5<-merge(outputMACSE, annotate4, by = "ID")
# Working Annotation File to Match Sequences
write.csv(annotate5, "Annotations_Working.csv",row.names = F)

# Subsampling to correct for sequence size
subsample<-as.data.frame(t(annotate5[,12:ncol(annotate5)]))
colnames(subsample)<-annotate5$ID
subsample2<-rrarefy.perm(subsample, n=100, round.out = T)
subsample3<-as.data.table(subsample2)
# Removing any zero OTUs columns if any
subsample4<-subsample3[,colSums(subsample3 != 0) > 0, with = F]
subsample4$Sample<-rownames(subsample2)
subsample5 <- subsample4 %>% dplyr::select(Sample, everything())
# Now having cleaned working data frame of samples and sequences
write.csv(subsample5, "ARMS_Data_Working.csv", row.names = F)


      #><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
      #  Table S1 and S2 
      #><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
      
      # Raw Sequences and MOTUs
      PreSeq<-as.data.frame(colSums(MOTU_samples))
      colnames(PreSeq)[1]<-"RawSampleSequences"
      PreSeq$RawMOTU<-specnumber(MOTU_samples, MARGIN = 2)
      # Post Filtration Sequences and MOTUs
      PostRelSumSeq<-as.data.frame(colSums(MOTU4[,2:ncol(MOTU4)]))
      colnames(PostRelSumSeq)[1]<-"PostFiltrationSequences"
      PostRelSumSeq$PostFiltrationMOTU<-specnumber(MOTU4[,2:ncol(MOTU4)], MARGIN = 2)
      MetaAlg<-as.data.frame(colSums(annotate3[,12:ncol(annotate3)]))
      colnames(MetaAlg)[1]<-"MetaAlgSeq"
      MetaAlg$MOTUMetaAlg<-specnumber((annotate3[,12:ncol(annotate3)]), MARGIN = 2)
      # Classified to Phylum Sequences and MOTUs
      Classify<-as.data.frame(colSums(annotate5[,12:ncol(annotate5)]))
      colnames(Classify)[1]<-"ClassifyPhylumSequences"
      Classify$ClassifyMOTU<-specnumber((annotate5[,12:ncol(annotate5)]), MARGIN = 2)
      #Final Sequences and MOTU richness post rarefraction
      rarefiedSeq<-as.data.frame(rowSums(subsample5[,2:ncol(subsample5)]))
      colnames(rarefiedSeq)[1]<-"SequencesPostRarefied"
      rarefiedSeq$MOTU_Richness<-specnumber((subsample5[,2:ncol(subsample5)]))
     
      # Binding all of this together into one dataframe to export to create Table S1 and S2
      TabS1<-cbind(PreSeq,PostRelSumSeq, MetaAlg,Classify,rarefiedSeq)
      TabS1$Sample<-rownames(TabS1)
      fwrite(TabS1, "TableS1_S2.csv", row.names = F)
      
     
      #><><><><><><><><><><><><><><><><><><><><><><><><><><><>><><><><><><><>
      #  Sequence Reads and MOTUs associated with taxonomic groups - Table S3
      #><><><><><><><><><><><><><><><><><><><><><><><><><><><>><><><><><><><>
      #Higher Taxa Group
      SumAll<-annotate2
      SumAll$Sum<-SumAll$Sum<-rowSums(SumAll[,12:ncol(SumAll)])
      KingdomSum<-aggregate(Sum ~ FinalKingdom, data = SumAll, FUN = sum)
      TotalKingdomSum<-sum(KingdomSum$Sum)
      KingdomSum$SequencePercent<-(KingdomSum$Sum/TotalKingdomSum)*100
      # Metazoan and Macroalgae, classified and unclassified
      MetaPlan<-KingdomSum[KingdomSum$FinalKingdom == "Metazoa" | KingdomSum$FinalKingdom == "Plantae",]
      sum(MetaPlan$SequencePercent) # 77.4%
      # Getting sum of PhylaSeq Post rarefraction
      PhySeqRare<-subsample5
      PhySeqRare1<-as.data.frame(t(PhySeqRare[,2:ncol(PhySeqRare)]))
      colnames(PhySeqRare1)<-PhySeqRare$Sample
      PhySeqRare1$ID<- rownames(PhySeqRare1)
      PhySeqRare2<-merge(annotate1,PhySeqRare1, by = "ID")
      PhySeqRare2$Sum<-rowSums(PhySeqRare2[,12:ncol(PhySeqRare2)])
      PhylaSequences<-aggregate(Sum ~ FinalPhylum, data = PhySeqRare2, FUN = sum)
      PhylaSequences$tally<- count(PhySeqRare2, FinalPhylum)
      
      #><><><><><><><><><><><><><><><><><><><><><><><><><><><>><><><><>
      #  Figure S4 - Species Annotations
      #><><><><><><><><><><><><><><><><><><><><><><><><><><><>><><><><>      
      
            # Table S4 was produced based on the MOTUs identified to species 
            # from the annotate5 dataframe ("ARMS_Data_Working.csv")

      #><><><><><><><><><><><><><><><><><><><><><><><><><><><>><><><><>
      #  Figure S3 - MOTU annotations across classification levels
      #><><><><><><><><><><><><><><><><><><><><><><><><><><><>><><><><>

      freqTaxa<- annotate5
      freqTaxa1<-count(freqTaxa, TaxaFrequency)
      sum(freqTaxa1$n)

      TaxaOrder<-c("Species", "Genus","Family","Order","Class","Phylum")
      freqTaxa1$Taxa<-TaxaOrder
      freqTaxa1$Taxa <- factor(freqTaxa1$Taxa, levels = TaxaOrder)

      ggplot(freqTaxa1, aes(x = Taxa, y = n,fill=Taxa))+ 
         geom_bar(aes(fill=Taxa),stat="identity", colour = "black",lwd=.5) +
        scale_fill_manual(values = c(Species ="#020202",
                                      Genus ="#383434",
                                      Family ="#606060",
                                      Order = "#999898",
                                      Class = "#C1BFBF",
                                      Phylum = "#E5E5E5"))+
        xlab("")+
        ylab("Number of MOTUs")+
        scale_y_continuous(expand = c(0, 0), limits = c(0, 280), n.breaks = 15)+
        theme_bw()+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              legend.position = "none",
              panel.background = element_blank(),
              axis.title.y = element_text(size = 14)) 

        ## Took this Graph into Illustrator to make FigureS3 by copying bars and stacking 
        ## them accordingly.

