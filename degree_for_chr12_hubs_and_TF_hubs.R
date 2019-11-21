#------
#title: "degree for chr12 hubs and TF hubs"
#author: xiaowei
#date: November 20, 2019
#Usedatabase:bigbin
#OutputCSV:
#       chr12_hubs
#       TF_hubs
#------
#*******************************
#************Outline************
#*******************************

#1. List of hubs for chr12 (a hub is a node with a high degree)
#1.1 Get hubs for chr12 (a hub is degree > 50)
#1.2 plot top 20 of degree of chr12 hubs 
#1.3 plot frequency and distribution of degree for chr12 hubs 

#2. List of hubs for TF (a hub is a node with a high degree)
#2.1 get hubs for TF (a hub is degree > 50)
#2.2 plot top 20 of degree of TF hubs 
#2.3 plot frequency and distribution of degree for TF hubs 


##########################################################################
#connect neo4j with RNeo4j
library(RNeo4j)
graph = startGraph("http://localhost:7474/db/data/", username="neo4j", password="xiaowei")


##########################################################################
#1. List of hubs for chr12 (a hub is a node with a high degree)
##########################################################################
#========================================================================
#1.1 Get hubs for chr12 (a hub is degree > 50)
#========================================================================
#Use the whole database graph to make degree algorithm, and return chr12 nodes
query = "
CALL algo.degree.stream(
                        'MATCH (allnodes) return allnodes',
                        'Interaction | Inclution | Bind',
                        {direction: 'both'})
YIELD nodeId, score
with nodeId, score
MATCH (chr12:chr12)
where id(chr12) = nodeId and score > 50
RETURN chr12.Name AS Name, score as degree
ORDER By score DESC
"
chr12_hubs <- cypher(graph, query)

head(chr12_hubs)
dim(chr12_hubs)
write.table(chr12_hubs, file = "chr12_hubs.csv",sep=",",quote = F, col.names = TRUE, row.names = FALSE)
#========================================================================
#1.2 plot top 20 of degree of chr12 hubs 
#========================================================================
chr12_distribution <- as.data.frame(table(chr12_hubs$degree))
library(ggplot2)
chr12_hubs_20 <- chr12_hubs[1:20,]
chr12_hubs_20$Name <- factor(chr12_hubs_20$Name, levels = chr12_hubs_20$Name)
ggplot() + 
  geom_bar(aes(x = chr12_hubs_20$Name, y = chr12_hubs_20$degree),
           stat="identity", fill = "lightblue")+ 
  labs(x="Relationship",y = "count")+
  ggtitle("Chromosome Inclusion Chromosome Range") +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90 ),
        plot.title = element_text(face ="bold" ,hjust = 0.5))

#========================================================================
#1.3 plot frequency and distribution of degree for chr12 hubs 
#========================================================================
ggplot(chr12_hubs, aes(x=degree))+ 
  geom_histogram(binwidth=30,colour="black", fill = "lightblue")+ #distribution plot
  geom_vline(aes(xintercept=mean(degree, na.rm=T)), color="red", linetype="dashed", size=1)+   #mean
  #ggtitle("Ferquency of hubs' degree in chromosome 12") +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90 ),
        plot.title = element_text(face ="bold" ,hjust = 0.5))

ggplot(chr12_hubs, aes(x=degree)) + 
  geom_histogram(aes(y=..density..),binwidth=30,colour="black", fill = "lightblue")+ #distribution plot
  geom_density(alpha=.2, color="red")+  #density line,Smoothed density estimates
  geom_vline(aes(xintercept=mean(degree, na.rm=T)), color="red", linetype="dashed", size=1)+  #mean
  #ggtitle("Distribution of hubs' degree in chromosome 12")+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90 ),
        plot.title = element_text(face ="bold" ,hjust = 0.5))

##########################################################################
#2. List of hubs for TF (a hub is a node with a high degree)
##########################################################################

#========================================================================
#2.1 get hubs for TF (a hub is degree > 50)
#========================================================================

#Use the whole database graph to make degree algorithm, and return TF nodes
query = "
CALL algo.degree.stream(
                        'MATCH (allnodes) return allnodes',
                        'Interaction | Inclution | Bind',
                        {direction: 'both'})
YIELD nodeId, score
with nodeId, score
MATCH (TF:TF)
where id(TF) = nodeId and score > 50
RETURN TF.Name AS Name, score as degree
ORDER By score DESC
"
TF_hubs <- cypher(graph, query)
write.table(TF_hubs, file = "TF_hubs.csv",sep=",",quote = F, col.names = TRUE, row.names = FALSE)

#========================================================================
#2.2 plot top 20 of degree of TF hubs 
#========================================================================

TF_hubs_20 <- TF_hubs[1:20,]
TF_hubs_20$Name <- factor(TF_hubs_20$Name, levels = TF_hubs_20$Name)
ggplot() + 
  geom_bar(aes(x = TF_hubs_20$Name, y = TF_hubs_20$degree),
           stat="identity", fill = "lightblue")+ 
  labs(x="node",y = "degree")+
  ggtitle("Top 20 of degree of TF nodes")+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90 ),
        plot.title = element_text(face ="bold" ,hjust = 0.5))

#========================================================================
#2.3 plot frequency and distribution of degree for TF hubs 
#========================================================================

TF_distribution <- as.data.frame(table(TF_hubs$degree))
library(ggplot2)
#ggplot(distribution, aes(x=Freq))+ geom_histogram(binwidth=.5)

ggplot(TF_hubs, aes(x=degree))+ 
  geom_histogram(binwidth=300,colour="black", fill = "lightblue")+ #distribution plot
  geom_vline(aes(xintercept=mean(degree, na.rm=T)), color="red", linetype="dashed", size=1)+   #mean
  ggtitle("Ferquency of degree in TF") +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90 ),
        plot.title = element_text(face ="bold" ,hjust = 0.5))


ggplot(TF_hubs, aes(x=degree)) + 
  geom_histogram(aes(y=..density..),binwidth=300,colour="black", fill = "lightblue")+ #distribution plot
  geom_density(alpha=.2, color="red")+  #density line,Smoothed density estimates
  geom_vline(aes(xintercept=mean(degree, na.rm=T)), color="red", linetype="dashed", size=1)+  #mean
  ggtitle("Distribution of hubs' degree in TF")+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90 ),
        plot.title = element_text(face ="bold" ,hjust = 0.5))

