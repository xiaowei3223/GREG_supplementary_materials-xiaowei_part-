#------
#title: "Data profiling in GREG(bigbin)"
#author: xiaowei
#date: November 19, 2019
#Usedatabase:bigbin
#OutputCSV:
#       labels
#       rel_details
#------

#*******************************
#************Outline************
#*******************************
#1. Number of nodes of each type (TF, chr1, etc)
    #1.1 Count all nodes
    #1.2 Count nodes of each type labels 
    #1.3 plot numbers of labels nodes 

#2. Number of relationships of each type (bind, etc)
    #2.1 Count all relationships
    #2.2 Count each type relationships 
    #2.3 Count relationships between two kinds of labels
    #2.4 count and plot numbers of relationship type
    #2.5 plot counts of relationship TF Bind chromosome
    #2.6 plot counts of relationship LncRNA Bind chromosome
    #2.7 plot counts of relationship chromosome Inclusion chromosome Range


##########################################################################
#connect neo4j with RNeo4j
library(RNeo4j)
graph = startGraph("http://localhost:7474/db/data/", username="neo4j", password="xiaowei")


##########################################################################
#1. Number of nodes of each type (TF, chr1, etc)
##########################################################################

#========================================================================
#1.1 Count all nodes
#========================================================================
totalnodes <- cypher(graph, query ="MATCH (n) RETURN count(n) as AllNodes" )
totalnodes
#========================================================================
#1.2 Count nodes of each type labels 
#========================================================================
query = "MATCH (n) RETURN distinct labels(n) as label, count(labels(n)) as count"
labels <- cypher(graph,query) 
labels
write.csv(labels, file = "labels.csv", sep=",",quote = F, col.names = TRUE, row.names = FALSE)
#========================================================================
#1.3 plot numbers of labels nodes 
#========================================================================
labels1 = labels
labels1$label <- factor(labels1$label，levels = c(labels1$label))
par(mfrow=c(1,2))

p = function(data){
  library(ggplot2)
  ggplot() + 
    geom_bar(aes(x = data$label, y = data$count),
             stat="identity", fill = "lightblue")+ 
    labs(x="node type",y = "count")+ 
    theme(axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90 ),
          plot.title = element_text(face ="bold" ,hjust = 0.5))
}
par(mfrow=c(1,2))
p1 <- p(labels1[c(1:25,50),]) + ggtitle("Numbers of Nodes with Chromosome Labels ")
p2 <- p(labels1[26:49,]) + ggtitle("Numbers of Nodes with Chromosome Range Labels")

cowplot::plot_grid(p1, p2, labels = c("A", "B"), ncol = 1)

##########################################################################
#2. Number of relationships of each type (bind, etc)
##########################################################################
#========================================================================
#2.1 Count all relationships
#========================================================================
query = "MATCH ()-[r]->() RETURN count(r) as namber_of_all_relationships;"
RNeo4j::cypher(graph,query)

#========================================================================
#2.2 Count each type relationships 
#========================================================================
query = "MATCH ()-[r]->() RETURN type(r) as relationshipType, count(*) as count;"
relationshipType <- RNeo4j::cypher(graph,query)
relationshipType

#========================================================================
#2.3 Count relationships between two kinds of labels
#========================================================================
query <- "MATCH (A)-[r]->(B) 
return labels(A) AS A, 
        type(r) AS relationshipType, 
        labels(B) AS B ,
        count(r) AS count"
rel_details <- cypher(graph, query)
write.csv(rel_details, file = "Count_relationships_details.csv", sep=",",quote = F, col.names = TRUE, row.names = FALSE)
#========================================================================
#2.4 count and plot numbers of relationship type
#========================================================================

BindRel <- rel_details[which(rel_details$relationshipType == "Bind"), ] 

#-----TF Binding relationship-----
BindRel_TF <- BindRel[which(BindRel$A == "TF"),]
BindRel_TFSUM <- sum(BindRel_TF$count)

#-------LncRNA Bind relationship---------------
BindRel_LncRNA <- BindRel[which(BindRel$A == "LncRNA"),]
BindRel_LncRNASUM <- sum(BindRel_LncRNA$count)

#-----Interaction relationship --------------
InteractionRel <- rel_details[rel_details$relationshipType =="Interaction",]

#---------TF - TF interaction relationships-------------
TFInteractionRel <- InteractionRel[InteractionRel$A == "TF",]
TFInteractionRelSUM <- sum(TFInteractionRel$count)

#--------chr_Range - chr_Range raltionships-------------
chr_RangeRel <- InteractionRel[InteractionRel$A != "TF",]
chr_RangeRelSUM <- sum(chr_RangeRel$count)

#--------numbers of relationships type------------------------
relationship_Type = c('TF binding ', 'lncRNA binding',
                      'TF-TF interactions', 'DNA-DNA interactions')
SumCount = c(BindRel_TFSUM,BindRel_LncRNASUM ,TFInteractionRelSUM, chr_RangeRelSUM)
RelSUM = data.frame(relationship_Type, SumCount)

#-------plot numbers of relationships type--------------------

ggplot() + 
  geom_bar(aes(x = RelSUM$relationship_Type, y = RelSUM$SumCount),
           stat="identity", fill = "lightblue", width = 0.4,position=position_dodge(0.01))+ 
  labs(x="relationship type",y = "Sum count")+
  ggtitle("Numbers of relationships type") + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90 ),
        plot.title = element_text(face ="bold" ,hjust = 0.5))+
  geom_text(aes(x = RelSUM$relationship_Type, y = RelSUM$SumCount, label=RelSUM$SumCount), size = 2.5,hjust= 0.5, vjust = -0.5 )

#========================================================================
#2.5 plot counts of relationship TF Bind chromosome
#========================================================================
#-------TF Bind chromosome---------------
library(ggplot2)
chr_level = c(paste0("chr",1:22),"chrX", "chrY")
BindRel <- rel_details[which(rel_details$relationshipType == "Bind"), ] 
BindRel_TF <- BindRel[which(BindRel$A == "TF"),]
BindRel_TF$B <- factor(BindRel_TF$B, levels = chr_level)

p3 <- ggplot() + 
  geom_bar(aes(x = BindRel_TF$B, y = BindRel_TF$count),
           stat="identity", fill = "lightblue")+ 
  labs(x="label",y = "count")+
  ggtitle("TF Bind Chromosome") + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90 ),
        plot.title = element_text(face ="bold" ,hjust = 0.5))

p3


#========================================================================
#2.6 plot counts of relationship LncRNA Bind chromosome
#========================================================================
#-------LncRNA Bind chromosome---------------
BindRel_LncRNA <- BindRel[which(BindRel$A == "LncRNA"),]
BindRel_LncRNA$B <- factor(BindRel_LncRNA$B, levels = chr_level)

p4 <- ggplot() + 
  geom_bar(aes(x = BindRel_LncRNA$B, y = BindRel_LncRNA$count),
           stat="identity", fill = "lightblue")+ 
  labs(x="label",y = "count")+
  ggtitle("LncRNA Bind Chromosome") + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90 ),
        plot.title = element_text(face ="bold" ,hjust = 0.5))

p4

#========================================================================
#2.7 plot counts of relationship chromosome Inclusion chromosome Range
#========================================================================
#--------chr inclusion chr_Range----------------------------
InclusionRel <- rel_details[rel_details$relationshipType =="Inclusion",]

InclusionRel$rel <- apply(InclusionRel[1:3], 1， FUN = function(X){ y = paste0(X[1],"-",X[2],"->",X[3]) })
chrRange_level <- c(paste0("chr",1:22,"-Inclusion->chr",1:22,"_Range"),
                    "chrX-Inclusion->chrX_Range",
                    "chrY-Inclusion->chrY_Range")
InclusionRel$rel <- factor(InclusionRel$rel, levels = chrRange_level)

ggplot() + 
  geom_bar(aes(x = InclusionRel$rel, y = InclusionRel$count),
           stat="identity", fill = "lightblue")+ 
  labs(x="Relationship",y = "count")+
  ggtitle("Chromosome Inclusion Chromosome Range") + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90 ),
        plot.title = element_text(face ="bold" ,hjust = 0.5))+
  coord_flip()+ #将条形图坐标轴转换过来变成横向条形图
  geom_text(aes(x = InclusionRel$rel, y = InclusionRel$count, label=InclusionRel$count), size = 2.5,hjust=-0.5 )


