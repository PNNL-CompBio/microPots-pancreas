##this script loads human data

library(dplyr)
source("synapseUtil.R")
syn<-synapseLogin()

#################################################Metadata processing

##read in metadata from synaspe
library(readxl)
isletMeta <- readxl::read_xlsx(syn$get('syn30070383')$path,sheet='Metadata table')%>%#"~/Desktop/2022/Paul/Hubmap/IsletMeta_New_correct.csv") %>% 
  mutate(X_Y = paste(`Islet Number`, "S", `X-coord`, `Y-coord`, sep = "_")) %>% 
  tibble::column_to_rownames("X_Y") %>% 
  mutate(`IsletNumber` = as.character(`Islet Number`))%>%
  dplyr::rename(IsletStatus='Islet status')%>%
  mutate(IsletOrNot=stringr::str_replace_all(IsletStatus,'Proximal|Distal','NonIslet'))

metadata <- isletMeta%>%
  select(IsletStatus,IsletOrNot,Plex,`Grid Number`)%>%
  tibble::rownames_to_column('Spot')


labels <- readxl::read_xlsx(syn$get('syn30070383')$path,sheet='Labels')%>%
  tidyr::pivot_longer(-Plex,names_to='OtherPlex',values_to='mwDist')

################################################Median centered data processing

crosstabs<-list(`0`='syn30389393',`10`='syn30385467',`1`='syn30385349',`2`='syn30385407',
                `3`='syn30385441',`4`='syn30385450',`7`='syn30385456')



##load in median-centered crosstabs into single file
crosstabList <- lapply(names(crosstabs), function(file_i) 
{
  islet_i <- file_i#sub("/.*", "", file_i)
  
  crosstab_i <- read.delim(syn$get(crosstabs[[file_i]])$path, row.names = 1)
  
  colnames(crosstab_i) <- paste0(islet_i, "_", colnames(crosstab_i))
  crosstab_i <- crosstab_i %>% tibble::rownames_to_column("feature")
  return(crosstab_i)
})



##create single matrix
crosstab <- crosstabList %>% purrr::reduce(full_join, by = "feature") %>% 
  tibble::column_to_rownames("feature") %>% 
  as.matrix()


normtab <- do.call(rbind,lapply(crosstabList,function(x){
  newx<-x%>%
    tidyr::separate(feature,into=c('id','upid','protein'),
                    sep='\\|')%>%
    dplyr::select(-c(id,upid))
  #ibble::column_to_rownames('protein')
  newx<-newx%>%
    tidyr::pivot_longer(2:ncol(newx),names_to='Spot', values_to='logRatio')%>%
    tidyr::separate(Spot,into=c('Image','S','Xcoord','Ycoord'),sep='_',remove=F)%>%
    dplyr::select(-S)
  newx
})
)%>%
  left_join(metadata)


################################################ ORiginal data processing

#there are also v2 crosstabs. these are NOT median centered originallly
v2_crosstabs<-list(`0`='syn31972408',`10`='syn31972599',`1`='syn31972430',`2`='syn31972462',
                   `3`='syn31972520',`4`='syn31972565',`7`='syn31972584')

##load in crosstabs into single file -without median
crosstabList2<- lapply(names(v2_crosstabs), function(file_i) 
{
  islet_i <- file_i#sub("/.*", "", file_i)
  
  crosstab_i <- read.delim(syn$get(v2_crosstabs[[file_i]])$path, row.names = 1)
  
  colnames(crosstab_i) <- paste0(islet_i, "_", colnames(crosstab_i))
  crosstab_i <- crosstab_i %>% tibble::rownames_to_column("feature")
  return(crosstab_i)
})




##create single matrix
crosstab2 <- crosstabList2 %>% purrr::reduce(full_join, by = "feature") %>% 
  tibble::column_to_rownames("feature") %>% 
  as.matrix()

fulltab <- do.call(rbind,lapply(crosstabList2,function(x){
  newx<-x%>%
    tidyr::separate(feature,into=c('id','upid','protein'),
                    sep='\\|')%>%
    dplyr::select(-c(id,upid))
  #ibble::column_to_rownames('protein')
  newx<-newx%>%
    tidyr::pivot_longer(2:ncol(newx),names_to='Spot', values_to='logRatio')%>%
    tidyr::separate(Spot,into=c('Image','S','Xcoord','Ycoord'),sep='_',remove=F)%>%
    dplyr::select(-S)
  newx
})
)%>%
  left_join(metadata)

##now save the fulltab to supp data table 1

write.table(fulltab,'suppTable1.csv',sep=',',quote=F,row.names=F)

isletMeta <- isletMeta[colnames(crosstab), ]


###############################utility functions
#' plotLeapR
#' plot leapR result
plotLeapR<-function(leapr.result,subCat=NA){
  
  if(!is.na(subCat))
    leapr.result<-leapr.result[grep(subCat,rownames(leapr.result)),]
  
  fisher=F
  
  if(length(grep('zscore',colnames(leapr.result)))==0){
    fisher=T
    leapr.result<-dplyr::rename(leapr.result,zscore='oddsratio')
  }
  
  pos<-leapr.result%>%
    subset(zscore>0)%>%
    dplyr::arrange(desc(zscore))%>%
    dplyr::mutate(Direction='Up')%>%
    subset(BH_pvalue>0)
  
  neg<-leapr.result%>%
    subset(zscore<0)%>%
    dplyr::arrange(zscore)%>%
    dplyr::mutate(Direction='Down')%>%
    subset(BH_pvalue>0)
  
  top20<-rbind(pos[1:10,],neg[1:10,])%>%
    subset(!is.na(zscore))%>%
    dplyr::arrange(zscore)%>%
    tibble::rownames_to_column('pathway')%>%
    dplyr::mutate(pathway=stringr::str_replace_all(pathway,'_',' '))
  
  if(!is.na(subCat))
    top20<-top20%>%
    dplyr::mutate(pathway=stringr::str_replace_all(pathway,subCat,''))
  
  top20$pathway<-factor(top20$pathway,levels=top20$pathway)
  if(fisher){
    top20<-top20%>%dplyr::rename(`Odds Ratio`='zscore')
  
    p<-top20%>%
      ggplot(aes(x=`Odds Ratio`,y=pathway,fill=Direction))+geom_bar(stat='identity')
  }
  else{
    p<-top20%>%
      ggplot(aes(x=zscore,y=pathway,fill=Direction))+geom_bar(stat='identity')
  }
  
  p

}







#' Created a function that fixed missing proteins by median expressio
#' and return matrix
correctMissingProteins<-function(fulltab){
  prot_avgs<-fulltab%>%
    dplyr::mutate(Grid=as.factor(`Grid Number`))%>%
    dplyr::group_by(protein)%>%
    dplyr::summarize(medExp=median(logRatio,na.rm=T))
  
  #now we deal with the larger crosstab matrix
  crosstab<-fulltab%>%
    dplyr::select(Spot,protein,logRatio)%>%
    tidyr::pivot_wider(names_from='Spot',values_from='logRatio')%>%
    tibble::column_to_rownames('protein')
  
  ##now we have to update the median expression again
  fixed.crosstab<-do.call('rbind',lapply(rownames(crosstab),function(x){
    avg<-subset(prot_avgs,protein==x)[,'medExp']
    protvals<-crosstab[x,]
    protvals[is.na(protvals)]<-unlist(avg)
    protvals
  }))%>%
    as.matrix()
  
  var0<-which(apply(fixed.crosstab,1,var)==0)
  if(length(var0)>0)
    fixed.crosstab<-fixed.crosstab[-var0,]
  return(fixed.crosstab)
}

###load go file
gosigs <- leapR::read_gene_sets('GO_Biological_Process_2021.txt')