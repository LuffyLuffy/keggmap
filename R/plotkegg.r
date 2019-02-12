#' A keggmap plot function
#'
#' This function allows you to plot kegg pathway based on xml file.
#' @param kegg pathway Which pathway you want to plot, plot_kegg('ec00010').
#' @keywords keggplot
#' @export
#' @examples
#' plot_kegg("ec00020")
plot_kegg<-function(pg){
  library(rvest)
  library(purrr)
  library(dplyr)
  library(stringr)
  #pg='ec00020'
  pg1<-system.file('extdata',paste(pg,'.xml',sep=''),package = 'keggmap')
  pg <- read_html(pg1)
  links <- html_nodes(pg, "reaction")
  df1<-bind_rows(lapply(xml_attrs(links), function(x) data.frame(as.list(x), stringsAsFactors=FALSE)))
  library(tidyverse)
  library(xml2)
  library(XML)

  x <- read_xml(pg1)
  doc <- xmlParse(pg1)
  # NUMBER OF reaction
  mtg_num <- length(xpathSApply(doc, "//reaction"))
  # DATAFRAME LIST OF EXPANDED MEETING ATTRS
  reaction_list <- lapply(seq(mtg_num), function(i) {
    races_num <- length(xpathSApply(doc, sprintf("//reaction[%s]/substrate", i)))
    #races_num=74
    data.frame(
      react_id = rep(xpathSApply(doc, sprintf("//reaction[%s]/@id", i)), races_num),
      name = rep(xpathSApply(doc, sprintf("//reaction[%s]/@name", i)), races_num),
      type= rep(xpathSApply(doc, sprintf("//reaction[%s]/@type", i)), races_num),
      substrate = rep(xpathSApply(doc, sprintf("//reaction[%s]/substrate/@name", i)), races_num),
      substrate_id = rep(xpathSApply(doc, sprintf("//reaction[%s]/substrate/@id", i)), races_num),
      product = rep(xpathSApply(doc, sprintf("//reaction[%s]/product/@name", i)), races_num),
      product_id = rep(xpathSApply(doc, sprintf("//reaction[%s]/product/@id", i)), races_num)
    )
  })
  final_df <- do.call(rbind, reaction_list)
  #final_df <- cbind(do.call(rbind, reaction_list),XML:::xmlAttrsToDataFrame(getNodeSet(doc, "//reaction/substrate")))
  graphs<-cbind(XML:::xmlAttrsToDataFrame(getNodeSet(doc, "//entry")),
                XML:::xmlAttrsToDataFrame(getNodeSet(doc, "//entry/graphics")))

  mtg_num <- length(xpathSApply(doc, "//relation"))
  # DATAFRAME LIST OF EXPANDED MEETING ATTRS
  relation_list <- lapply(seq(mtg_num), function(i) {
    races_num <- length(xpathSApply(doc, sprintf("//relation[%s]/subtype", i)))
    #races_num=74
    data.frame(
      subtype_id = rep(xpathSApply(doc, sprintf("//relation[%s]/subtype/@name", i)), races_num),
      subtype_value = rep(xpathSApply(doc, sprintf("//relation[%s]/subtype/@value", i)), races_num)
    )
  })

  final_relation <-cbind(do.call(rbind, relation_list),XML:::xmlAttrsToDataFrame(getNodeSet(doc, "//relation")))
  graphs$x<-as.numeric(as.character(graphs$x))
  graphs$y<-as.numeric(as.character(graphs$y))
  graphs$y<- - graphs$y
  graphs$id<-as.character(graphs$id)
  pos<-graphs[,c('id','x','y')]
  pro_sub<-data.frame(apply(final_df, MARGIN = 2, as.character),stringsAsFactors = F)
  df1<-merge(pro_sub,pos,by.x='substrate_id',by.y = 'id')
  df2<-merge(df1,pos,by.x='product_id',by.y='id')
  #library(ggplot2)
  #library(grid)
  graphs<-graphs[,-c(2,3)]

  #load('cpd.names.rda')
  temp<-merge(graphs,cpd.names,by.x='name',by.y='ACCESSION_NUMBER',all=T)
  temp1<-as.data.frame(graphs$name[grepl('^C\\d{4}',graphs$name)])
  temp1<-temp1[!duplicated(temp1[,1]),,drop=FALSE]
  rownames(temp1)<-temp1[,1]
  colnames(temp1)<-'cpd'
  temp2<-merge(temp1,cpd.names,by.x='cpd',by.y='ACCESSION_NUMBER',all=TRUE)
  temp2<-temp2[!duplicated(temp2$cpd),c(1,5)]
  temp2<-na.omit(temp2)
  rownames(temp2)<-temp2$cpd
  temp2$cpd<-as.character(temp2$cpd)
  temp2$NAME<-as.character(temp2$NAME)
  graphs$name<-as.character(graphs$name)


  for(i in 1:nrow(graphs)){
    if(graphs$name[i] %in% temp2$cpd){graphs$name[i] <- temp2[graphs$name[i],'NAME']}
  }

  #graphs<-graphs[!duplicated(graphs),]

  entry12<-data.frame(apply(final_relation, MARGIN = 2, as.character),stringsAsFactors = F)
  df3<-merge(entry12,pos,by.x='entry1',by.y = 'id')
  df4<-merge(df3,pos,by.x='entry2',by.y='id')
  df4$name<-df4$subtype_value
  df5<-merge(df4,pos,by.x='subtype_value',by.y='id')
  temp<-graphs$name
  graphs$name<-str_wrap(graphs$name,width=25)
  ggplot(graphs, aes(x, y, label = name))+
    geom_segment(data=df2,mapping = aes(x=x.x, y=y.x, xend=x.y, yend=y.y),
                 lineend = 'butt',linejoin = 'bevel',
                 arrow=arrow(type='closed',length = unit(0.1, "inches")),size=1,color='black')+
    geom_segment(data=df5,mapping = aes(x=x.x, y=y.x, xend=x, yend=y,linetype=type),
                 lineend = 'butt',linejoin = 'bevel',
                 arrow=arrow(type='closed',length = unit(0.1, "inches")),size=1,color='black')+
    geom_segment(data=df5,mapping = aes(x=x, y=y, xend=x.y, yend=y.y,linetype=type),
                 lineend = 'butt',linejoin = 'bevel',
                 arrow=arrow(type='closed',length = unit(0.1, "inches")),size=1,color='black')+
    geom_label(aes(x,y),data=graphs[graphs$type=='rectangle',])+
    geom_label(aes(x,y),color='darkgreen',data=graphs[graphs$type=='roundrectangle',])+
    geom_point(aes(x,y),shape=16,data=graphs[graphs$type=='circle',],
               color="blue",position = position_dodge(0.5), size=3)+
    geom_text(aes(x,y,label=name),data=graphs[graphs$type=='circle',],
              vjust=-0.8,hjust=0, color="blue",position = position_dodge(0.5), size=3)+
    xlim(min(graphs$x)*0.8,max(graphs$x)*1.1)+
    theme_void()
}

