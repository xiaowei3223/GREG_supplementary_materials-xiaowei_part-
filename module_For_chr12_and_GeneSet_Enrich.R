#------
#title: "make module for chr12 and Gene set enrich"
#author: xiaowei
#date: November 19, 2019
#Usedatabase:bigbin
#Output:
#       dataframe included gene names,nodeid and node labels: AllgeneName, AllgeneName_TF_LncRNA,  LncRNAName,  TFName
#       module for chr12: mod1, mod2
#       relationships between nodes in the same module: mod2_Relationships
#       target nodes of target genes: TargetNodes,  TargetNodes2 
#       Gene set enrich GO and KEGG results: mod_enrich
#function:
#       getAllGeneName: to get all TF, LncRNA, genes in chromosome
#       geneMatchIndex: To map nodes of gene names, TF names or lncRNA names
#       gene_match_mod: To get modules which included nodes of target genes/TF/LncRNA
#       Module_Relationship: To get relationships between nodes in the same module
#------


#*******************************
#************Outline************
#*******************************

#1  Make modules for chr12
#   1.1How many module for chr12 by algo.louvain.stream
#2  Find out nodes of target genes/TF/LncRNA
    ##2.1 Step 1: Build a function to get all TF, LncRNA, genes in chromosome
    ##2.2 Step 2: Get all TF, LncRNA, genes in chromosome 
    ##2.3 Step 3: Build a function `geneMatchIndex` to find out nodes of target genes/TF/LncRNA
    ##2.4 Step 4: Find out nodes of target genes/TF/LncRNA
    ##2.5 Step 5: To map gene names, TF names, lncRNA in each modules

#3  To get modules which included nodes of target genes/TF/LncRNA
    ##3.1 Step 1: Make a function to judege modules which included nodes of target genes/TF/LncRNA
    ##3.2 Step 2: Get the modules included target genes/TF/LncRNA

#4  Get all relationships in modules which included nodes of target genes/TF/LncRNA
    ##4.1 Step 1: build a function to get relationships during nodes in each community(module)
    ##4.2 Step 2: get relationship of mod2

#5  Gene Set enrich
    ##5.1 enrich GO
    ##5.2 enrich KEGG

##########################################################################
#1  Make modules for chr12
##########################################################################

#connect neo4j with RNeo4j
library(RNeo4j)
graph = startGraph("http://localhost:7474/db/data/", username="neo4j", password="xiaowei")

# cypher command to get modules for chr12
query = "
        // to make louvain modules for chr12
        CALL algo.louvain.stream(
          'MATCH (n:chr12)--(m) with id(n) as chr12,id(m) AS m   MATCH (al) where id(al) = m or  id(al) = chr12 with distinct al As al return id(al) as id',
          'MATCH (n:chr12)-[r]-(m) RETURN id(n) AS source, id(m) AS target',
          {graph:'cypher',direction: 'both'})
        YIELD nodeId, community
        
        MATCH (n)
        where id(n) = nodeId
        
        with collect (id(n)) AS AllnodeId, 'module'+community AS community, collect(n.Name) AS AllNodeNames
        RETURN AllnodeId, community, AllNodeNames, size(AllnodeId) AS size_of_allnode
        "

#use function `cypher` in {RNeo4j} to run cypher command and return the results as table(data.frame) 
mod1 = cypher(graph, query)
class(mod1)
names(mod1)
#AllnodeId just is columns of nodes ids in each module.
#community is columns of modules names.
#AllNodeNames is columns for all nodes those have Name properties.
#size_of_allnode is calculated by AllnodeId.
#========================================================================
#1.1 How many module for chr12 by algo.louvain.stream
#========================================================================
nrow(mod1)

##########################################################################
#2  Find out nodes of target genes/TF/LncRNA
##########################################################################

#========================================================================
##2.1 Step 1: Build a function to get all TF, LncRNA, genes in chromosome
#========================================================================
#function format:getAllGeneName(node_labels)
# node_labels: one of "genesInchr", "TF", "LncRNA", "all", "chr1", "chr2", "chr3" ... "chr22", "chrX" and "chrY"

#values:
# a data.frame which included geneName, NodeID, NodeLabel and NodeNames
library(RNeo4j)
graph = startGraph("http://localhost:7474/db/data/", username="neo4j", password="xiaowei")
getAllGeneName = function(node_labels){
  
  #get all nodes of gene names in chromosome
  if (node_labels == 'genesInChr'){
  query = "
      MATCH (n)
      where exists(n.Details)
      with split(n.Details, ';') AS Details,id(n) AS NodeID, labels(n) AS NodeLabel, n.Name AS NodeNames
      UNWIND Details AS details
      with details, NodeID, NodeLabel, NodeNames
      where details contains 'gene_name'
      with details as genes, NodeID, NodeLabel, NodeNames
      return substring(genes,11) as geneName, NodeID, NodeLabel, NodeNames
  "
  allGenes <- RNeo4j::cypher(graph, query) 
  }
  
  #get all nodes of TF names 
  else if (node_labels == 'TF'){
    query = "MATCH (n:TF) RETURN n.Name AS geneName, id(n) AS NodeID, labels(n) AS NodeLabel, n.Name AS NodeNames"
    allGenes <- RNeo4j::cypher(graph, query) }
  
  #get all nodes of LncRNA names
  else if (node_labels == 'LncRNA'){ 
    query = "MATCH (n:LncRNA)  RETURN n.Name AS geneName, id(n) AS NodeID, labels(n) AS NodeLabel, n.Name AS NodeNames"
    allGenes <- RNeo4j::cypher(graph, query) }
  
  #get all nodes that contains TF, lncRNA and gene names
  else if (node_labels == 'all'){
    query = "
      MATCH (n)
      where exists(n.Details)
      with split(n.Details, ';') AS Details,id(n) AS NodeID, labels(n) AS NodeLabel, n.Name AS NodeNames
      UNWIND Details AS details
      with details, NodeID, NodeLabel, NodeNames
      where details contains 'gene_name'
      with details as genes, NodeID, NodeLabel, NodeNames
      return substring(genes,11) as geneName, NodeID, NodeLabel, NodeNames
  "
    GenesInchr <- RNeo4j::cypher(graph, query) 
    
    query = "MATCH (n:TF) RETURN n.Name AS geneName, id(n) AS NodeID, labels(n) AS NodeLabel, n.Name AS NodeNames"
    TF <- RNeo4j::cypher(graph, query)
    
    query = "MATCH (n:LncRNA)  RETURN n.Name AS geneName, id(n) AS NodeID, labels(n) AS NodeLabel, n.Name AS NodeNames"
    LncRNA <- RNeo4j::cypher(graph, query)
    
    allGenes <- rbind( TF, LncRNA, GenesInchr)
  }
  
  # get one chromosome nodes of genes names
  else{
    query = paste0("MATCH (n:", 
                   node_labels,
                   ")
      where exists(n.Details)
      with split(n.Details, ';') AS Details,id(n) AS NodeID, labels(n) AS NodeLabel, n.Name AS NodeNames
      UNWIND Details AS details
      with details, NodeID, NodeLabel, NodeNames
      where details contains 'gene_name'
      with details as genes, NodeID, NodeLabel, NodeNames
      return substring(genes,11) as geneName, NodeID, NodeLabel, NodeNames
  ")
    allGenes <- RNeo4j::cypher(graph, query ,chr_label = node_labels) 
  }
  allGenes = unique( allGenes )
  
  #return results
  return(allGenes)
}


#========================================================================
##2.2 Step 2: Get all TF, LncRNA, genes in chromosome 
#========================================================================
AllgeneName_TF_LncRNA <- getAllGeneName('all')  #get all nodes that contains TF, lncRNA and gene names
TFName <- getAllGeneName('TF')                  #get all nodes of TF names 
LncRNAName <- getAllGeneName('LncRNA')          #get all nodes of LncRNA names
AllgeneName <- getAllGeneName('genesInChr')     #get all nodes of gene names in chromosome
chr1 <- getAllGeneName('chr1')                  #get one chromosome nodes of genes names

#========================================================================
##2.3 Step 3: Build a function `geneMatchIndex` to find out nodes of target genes/TF/LncRNA
#========================================================================
#The format of `geneMatchIndex` is:  geneMatchIndex(gene, allGeneNode, FullMatch = 1, by = "geneName")

#gene： a vector, what genes you want to match.
#allGeneNode: allGeneNode is data you chose which one as your match. It must is the results from getAllGeneName function.
#FullMatch : Default is 1. please set argument FullMatch = 0 ,if you want use fuzzy match for those gene names contain strings of your input.
#by : Default is "geneName". one type of "geneName", "NodeID", "NodeLabel" and "NodeNames", it should be as your genes types.

#values:
#A list contains `TargetNodes` and `NoGeneIndex`.
#`TargetNodes`` is a frame data of those matched in allGeneNode. `NoGeneIndex`` is a vector included those didnot find out in allGeneNode.


geneMatchIndex = function(gene, allGeneNode, FullMatch = 1, by = "geneName"){
  geneIndex = c() 
  NoGeneIndex =  c()
  Ag = allGeneNode[,by]
  if (FullMatch == 1){
    #----------which-----------------
    for (i in 1:length(gene) ) {
      geneIndex1 = which(Ag == gene[i] )
      if (length(geneIndex1) == 0){NoGeneIndex = c(gene[i], NoGeneIndex)}
      else{geneIndex = c(geneIndex,geneIndex1)}
    }
  }
  
  else{
    #----------grep---------------
    for (i in 1:length(gene) ) {
      geneIndex1 = grep(gene[i],Ag, fixed = TRUE)
      if (length(geneIndex1) == 0){NoGeneIndex = c(gene[i], NoGeneIndex)}
      else{geneIndex = c(geneIndex,geneIndex1)}
    }
  }
  
  TargetNodes <- allGeneNode[geneIndex, ]
  geneInGREG = list(TargetNodes, NoGeneIndex)
  names(geneInGREG) = c("TargetNodes", "NoGeneIndex")
  #返回匹配上的TargetNodes和没有找到基因NoGeneIndex
  return (geneInGREG )
}



#========================================================================
##2.4 Step 4: Find out nodes of target genes/TF/LncRNA
#========================================================================
gene <- c('NANOG', 'CTCF', 'EP300', 'POLR2A', 'YY1', 'RAD21', 'SMC3', 'STAG1', 'MED1', 'MED12')
#Genes must is in modules
MustIn = c("CTCF", "EP300","POLR2A","YY1")
#Genes may not in modules, but at least one gene in modules
OrIn1 = c("RAD21","SMC3","STAG1")
OrIn2 = c("MED1","MED12")

#get node which gene names is
geneIndex1 = geneMatchIndex('NANOG',AllgeneName)
geneIndex2 = geneMatchIndex(gene[-1],TFName)
TargetNodes <- rbind(geneIndex1$TargetNodes,geneIndex2$TargetNodes)
TargetNodes 


#========================================================================
##2.5 Step 5: To map gene names, TF names, lncRNA in each modules
#========================================================================
mod2 <- mod1

#mod1_match_gene to save all match results.
mod1_match_gene <- list() #用来存放每个模块nodeid匹配后对应的gene列表，里面包括
for (i in 1:nrow(mod2)){
  Mod_nodeId_list = unlist(mod1$AllnodeId[i])
  Mod_nodes_list = geneMatchIndex(Mod_nodeId_list, AllgeneName_TF_LncRNA, FullMatch = 1, by = "NodeID")
  Mod_gene_nodes_list = unique(Mod_nodes_list$TargetNodes) #得到每个模块nodeID匹配后基因的列表,加上一个unique来避免重复
  
  mod1_match_gene[[i]] <-  Mod_gene_nodes_list #保存匹配后的结果
  
  #Get all gene , TF and LncRNA names in each modules
  mod2$allGeneTFLncRNA[[i]] <- unique(Mod_gene_nodes_list$geneName) #得到每个模块对应的基因数量，加上一个unique来避免重复
  mod2$size_of_allGeneTFLncRNA[i] <- length(mod2$allGeneTFLncRNA[[i]])  #统计每个匹配后基因的数量
  
  #Get all chr12 names in each modules
  mod2$allchr12Name[[i]] <- unique(Mod_gene_nodes_list$NodeNames[which(Mod_gene_nodes_list$NodeLabel == 'chr12')]) #获取每个模块里chr12的节点名称,加上一个unique来避免重复
  mod2$size_of_allchr12Name[i] <- length(mod2$allchr12Name[[i]]) #统计每个模块匹配后chr12节点的数量
  
  #Get all gene names at chr12 in each modules
  mod2$allGeneAtchr12[[i]] <- unique(Mod_gene_nodes_list$geneName[which(Mod_gene_nodes_list$NodeLabel == 'chr12')])
  mod2$size_of_allGeneAtchr12[i] <- length(mod2$allGeneAtchr12[[i]])
  
  #Get all TF names in each modules
  mod2$AllTF[[i]] <- unique(Mod_gene_nodes_list$NodeNames[which(Mod_gene_nodes_list$NodeLabel == 'TF')]) ##获取每个模块里TF的节点名称,加上一个unique来避免重复
  mod2$size_of_AllTF[i] <- length(mod2$AllTF[[i]]) #统计每个模块匹配后TF节点的数量
}
names(mod1_match_gene) <- mod2$community

##########################################################################
#3  To get modules which included nodes of target genes/TF/LncRNA
##########################################################################

#========================================================================
##3.1 Step 1: Make a function to judege modules which included nodes of target genes/TF/LncRNA
#========================================================================
#gene_match_mod function format is: gene_match_mod(TargetNodes,mod,MustIn,OrIn1,OrIn2)

#TargetNodes: TargetNodes is a data frame of that results you used geneMatchIndex function.
#mod: just like mod2 in this tutorial object. A data frame contains a column `AllnodeId` for list of all nodes ids in each modules, 
#     and `community` column for each modules names.
#MustIn: A vector of all of gene names must be in module as you design.
#OrIn1: default as null, A vector of all of gene names may be in module as you design, at least one gene in module.
#OrIn2: default as null, same of OrIn1

#Values:
#Results is a list included `TargetNodes_mod`, `mod_targetGene`, `mod_want`.
#`TargetNodes_mod`: a data frame contains columns of "geneName", "NodeID", "NodeLabel" and "module". In "module" column tells you which 
#                   module this gene belong to.
#`mod_targetGene`: a data frame contains columns of "module", "targetGene" and "NoTargetGene"; "targetGene" columns tells you those genes
#                  are in this modules and "NoTargetGene" columns tells you those genes are not in this modules.
#`mod_want`: This is select from mod if there had modules you ask.



gene_match_mod <- function(TargetNodes,mod,MustIn,OrIn1=NULL ,OrIn2=NULL){
  TargetNodes_mod <- TargetNodes
  
  # find out genes in modules
  for (i in 1:length(TargetNodes$NodeID)){
    for (j in 1:nrow(mod)){
      if( TargetNodes$NodeID[i] %in% mod$AllnodeId[[j]]){ 
        TargetNodes_mod$module[i] = mod$community[j]  #module column in TargetNodes_mod will record which module included the gene
      }
    }
  }
  
  #show you what target genes in this module and not in this module
  #Those modules do not have one target gene node will not in the data.frame. 
  mod_targetGene <- aggregate(TargetNodes_mod$geneName, by = TargetNodes_mod['module'], FUN = paste0) 
  names(mod_targetGene)[2] <- 'targetGene'
  
  module_of_having_targetGene = unique(TargetNodes_mod$module)
  
  
  mod_want_index = c() #save those modules as you asked
  for (i in 1:length(module_of_having_targetGene)){
    mg = TargetNodes_mod[which( TargetNodes_mod$module == module_of_having_targetGene[i]),]$geneName #Those target genes in this modules
    Nmg = setdiff(gene,mg)# Those target genes not in this modules
    mod_targetGene$NoTargetGene[i] = paste0(Nmg,collapse = ', ')
    
    #target genes you asked in the same modules but actually not in this modules
    Nmg_MustIn = setdiff(MustIn,mg) #Nmg_MustIn > 0，Not you want 
    case1 = length(Nmg_MustIn) > 0 
    #target genes you asked maybe in the same modules but at least one in this modules
    if (length(OrIn1) > 0){
      mg_OrIn1 = intersect(OrIn1,mg) 
      case2 = length(mg_OrIn1) > 0 # mg_OrIn1 > 0, yes, that as you defined
    }else{case2 = TRUE }
    
    if (length(OrIn2) > 0){
      mg_OrIn2 = intersect(OrIn2,mg)# mg_OrIn2 > 0, yes, that as you defined
      case3 = length(mg_OrIn2) > 0 
    }else{case3 = TRUE }
    
    # If really as you defined, that get this modules
    if ( (!case1) & case2 &case3 ) { 
      mod_index1 = which(mod$community == module_of_having_targetGene[i])
      mod_want_index <- c(mod_index1 ,mod_want_index)
    }
  }
  
  mod_want <- mod[mod_want_index, ] 
  
  results <- list(TargetNodes_mod, mod_targetGene, mod_want)
  names(results) <- c('TargetNodes_mod', 'mod_targetGene', 'mod_want')
  
  return(results)
}

#========================================================================
##3.2 Step 2: Get the modules included target genes/TF/LncRNA
#========================================================================

#--------
#Example 1
#Select the modules that include NANOG AND CTCF AND EP300 AND POLR2A AND YY1 AND cohesin (RAD21, SMC3, OR STAG1) AND mediator (MED1 or MED12)
#--------

gene <- c('NANOG', 'CTCF', 'EP300', 'POLR2A', 'YY1', 'RAD21', 'SMC3', 'STAG1', 'MED1', 'MED12')
#Genes must is in modules
MustIn = c("CTCF", "EP300","POLR2A","YY1")
#Genes may not in modules, but at least one gene in modules
OrIn1 = c("RAD21","SMC3","STAG1")
OrIn2 = c("MED1","MED12")

#get nodes which gene names is
geneIndex1 = geneMatchIndex('NANOG',AllgeneName)
geneIndex2 = geneMatchIndex(gene[-1],TFName)
TargetNodes <- rbind(geneIndex1$TargetNodes,geneIndex2$TargetNodes)
TargetNodes 

#To judge which module included target genes
Select1 <-  gene_match_mod(TargetNodes,mod2,MustIn,OrIn1,OrIn2)
#check which modules target genes are in 
Select1$TargetNodes_mod
#check what target genes are in the module and not in the module
Select1$mod_targetGene


#--------
#Example 2
#Select the modules that include CTCF AND EP300 AND POLR2A AND YY1 AND cohesin (RAD21, SMC3, OR STAG1) AND mediator (MED1 or MED12 
#--------
gene2 <- c('CTCF', 'EP300', 'POLR2A', 'YY1', 'RAD21', 'SMC3', 'STAG1', 'MED1', 'MED12')
#Genes must is in modules
MustIn = c("CTCF", "EP300","POLR2A","YY1")
#Genes may not in modules, but at least one gene in modules
OrIn1 = c("RAD21","SMC3","STAG1")
OrIn2 = c("MED1","MED12")

#get nodes which gene names is
gene2Index = geneMatchIndex(gene[-1],TFName)
TargetNodes2 <- gene2Index$TargetNodes
TargetNodes2 

#To judge which module included target genes
Select2 <-  gene_match_mod(TargetNodes2,mod2,MustIn,OrIn1,OrIn2)
#check which modules target genes are in 
Select2$TargetNodes_mod
#check what target genes are in the module and not in the module
Select2$mod_targetGene
#check the module as you defined
Select2$mod_want

##########################################################################
#4  Get all relationships in modules which included nodes of target genes/TF/LncRNA
##########################################################################

#========================================================================
##4.1 Step 1: build a function to get relationships during nodes in each community(module)
#========================================================================
#This function format is: Module_Relationship(Module)
#Module: A data frame of contains the first column included list of all nodes ids in each modules and has a column 
#        names "community" as modules names.

#Values:
#Results is a list includes relationships data frames of during nodes in all modules.
#Each relationships data frame contains two nodes' id, labels, properties and relationship id, type, and properties.


#----------------建一个函数来获取每一个模块节点之间的关系--------------------------------------
#
#The argument is result we get  from step2
#module就是我们前面获取的结果
Module_Relationship = function(Module){
  
  #Module_Relationship as an object for save our results
  Module_Relationship = list()
  
  #To get relationships during nodes in the same modules
  #The results include nodes' id ,label and properties, and relationships' type, id and properties.
  for( i in 1:nrow(Module) ){
    #for the vector about all nodes id in the same modules
    Module_NodeID = unlist(Module[[1]][[i]])
    query = "
            // the first we should get all node in the same modules, result is a list
            MATCH (n)
            where id(n) in {Module_NodeID}
            with collect( id(n) ) AS Module_Node_ID
        
            // get relationships 
            MATCH (A)-[r]->(B)
            where id(A) in Module_Node_ID and id(B) in Module_Node_ID
        
            //at here, properties() function for all properties at an node or a relationship    
            RETURN id(A) AS fromNodeID, 
                  labels(A) AS fromNodeLabel, 
              
                  id(B) AS toNodeID,
                  labels(B) AS toNodeLabel,
              
                  id(r) AS Relationship_ID,
                  type(r) AS Relationship,
          
                  properties(A) AS AProperties,
                  properties(B) AS BProperties,
                  properties(r) AS RProperties
        "
    #run the query for results   
    #Module_Relationship[[i]]= cypher(graph, query, Module_NodeID = Module_NodeID )
    ABRelationship= cypher(graph, query, Module_NodeID = Module_NodeID )
    
    #Now, we split all properties of nodes using unnest_nodes() function in {neo4r}
    #A
    Ap = ABRelationship["AProperties"]
    names(Ap)[1] = "properties"
    App =  neo4r::unnest_nodes(Ap, what = "properties")
    names(App) = paste0("fromNode_", names(App) )
    
    #B
    Bp = ABRelationship["BProperties"]
    names(Bp)[1] = "properties"
    Bpp =  neo4r::unnest_nodes(Bp, what = "properties")
    names(Bpp) = paste0("toNode_", names(Bpp) )
    
    
    #split properties of relationships,we need make the same formats. If not, there show some wrong and we can not fixed.
    Rp = ABRelationship[,c("Relationship_ID","Relationship","fromNodeID","toNodeID","RProperties")]
    names(Rp) = c("id", "type", "stratNode", "endNode", "properties")
    Rpp =  neo4r::unnest_relationships(Rp)
    names(Rpp)[1:4] = c("Relationship_ID","Relationship","fromNodeID","toNodeID")
    names(Rpp)[c(-(1:4))] = paste0("R_", names(Rpp)[c(-(1:4))] )
    
    #union all data
    alldataRelationship = cbind(Rpp,App,Bpp, ABRelationship[,c("fromNodeLabel","toNodeLabel")])
    ColNames = c("fromNodeID","fromNodeLabel","toNodeID","toNodeLabel","Relationship_ID","Relationship",
                 names(App),
                 names(Bpp),
                 names(Rpp)[c(-(1:4))])
    Module_Relationship[[i]] = alldataRelationship[,ColNames]
  }
  
  names(Module_Relationship) = Module$community
  
  return(Module_Relationship)
}

#========================================================================
##4.2 Step 2: get relationship of mod2
#========================================================================
#Because we don't get modules as we defined, so here make an example for showing get relationship
mod2_Relationships =Module_Relationship(mod2)
head(mod2_Relationships[[1]])


##########################################################################
#5  Gene Set enrich
##########################################################################
# for each module, do gene set enrich analysis
mod_enrich <- data.frame(mod2[,'community'])
#BiocManager::install(c('org.Hs.eg.db','clusterProfiler'))

#========================================================================
##5.1 enrich GO
#========================================================================

library(clusterProfiler)
for (i in 1:nrow(mod2)){
  module_gene <- unlist(mod2$allGeneAtchr12[[i]])
  module_enrichGO <- enrichGO(module_gene, 'org.Hs.eg.db', keyType = 'SYMBOL', pvalueCutoff = 0.05)
  mod_enrich$enrichGO[i] <- paste0(module_enrichGO@result$ID, collapse = "; ")
  
}


#========================================================================
##5.2 enrich KEGG
#========================================================================
library(org.Hs.eg.db)
library(clusterProfiler)
for (i in 1:nrow(mod2)){
  module_gene <- unlist(mod2$allGeneAtchr12[[i]])
  geneIDselect <-select(org.Hs.eg.db, #.db是这个芯片数据对应的注释包
                        keys= module_gene,
                        columns=c("ENTREZID","ENSEMBL","GENENAME"), #clolumns参数是你要转换的ID类型是什么，这里选择三个。
                        keytype="SYMBOL" )#函数里面的keytype与keys参数是对应的，keys是你输入的那些数据，keytype是指这些数据是属于什么类型的数据。 
  module_enrichKEGG = enrichKEGG(geneIDselect$ENTREZID, organism = 'hsa', keyType = "kegg", pvalueCutoff = 0.05)
  mod_enrich$enrichKEGG[i] <- paste0(module_enrichKEGG@result$ID, collapse = "; ")
}

#=======All results of enrich GO / KEGG show you likes that:=============
head(mod_enrich)[1:2,]