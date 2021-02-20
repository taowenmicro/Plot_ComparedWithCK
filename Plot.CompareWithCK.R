

# j = "Phylum"


Plot.CompareWithCK <- function(
  ps = psdata,
  result = result,
  CK = "OE",
  j = "Genus",
  abun = 0.001){
  
  map = sample_data(ps)
  TF = c()
  for (i in 1: nrow(result)) {
    a <- result %>%
      select(ends_with("level")) %>%
      filter(row_number() == i) %>%
      as.matrix() %>%
      as.vector() %>%
      unique()
    a
    if (length(a) ==1) {
      if (a == "nosig") {
        TF[i] = 0
      }else {TF[i] = 1}
    } else {
      TF[i] = 1
    }
    
  }
  result$TF = TF
  
  
  
  result$mean <- result %>%
    # filter(TF == 1) %>%
    select(one_of(unique(map$Group))) %>%
    rowMeans() %>%
    as.vector()
  
  Sresult <- result %>%
    filter(TF == 1) %>%
    filter(mean > abun)
  
  
  
  Sresult$ID = row.names(Sresult)
  rank_names(ps)
  b <- colnames(
    result %>%
      select(ends_with("level"))
  )
  
  longda <- reshape2::melt(Sresult,
                           id.vars = c("ID",b),#需要保留不参与聚合的变量,
                           measure.vars = c(as.character(unique(map$Group))),#用于聚合的变量,
                           variable.name='treat',
                           value.name='abundance') %>%
    filter(treat != CK)
  level = c()
  for (i in 1:nrow(longda)) {
    level[i] <- longda[i,] %>%
      select(contains(as.character(longda$treat[i]))) %>%
      as.matrix() %>%
      as.vector()
  }
  longda$level = level
  
  ck <- reshape2::melt(Sresult,
                       id.vars = c("ID",b),#需要保留不参与聚合的变量,
                       measure.vars = c(as.character(unique(map$Group))),#用于聚合的变量,
                       variable.name='treat',
                       value.name='abundance') %>%
    filter(treat == CK)  %>%
    select("ID","abundance")
  
  colnames(ck)[2] = paste("CK",colnames(ck)[2],sep = "_")
  plotda <- longda %>% left_join(ck)  %>% 
    mutate(level2 = abundance - CK_abundance,.keep = "all") %>% arrange(ID)
  plotda$level2 <- plotda$level2 > 0
  plotda$abundance[plotda$level2 == F] = -plotda$abundance[plotda$level2 == F]
  
  
  
  Taxonomies_x = plyr::ddply(plotda,"ID", summarize,
                             label_sd = cumsum(abundance),
                             label_y = cumsum(abundance) - 0.5*abundance)
  
  plotdata <- cbind(plotda,Taxonomies_x[,-1]) 
  head(plotdata)
  
  
  plotdata$treat = factor(plotdata$treat,levels = as.character(unique(plotdata$treat)[4:1]))
  
  c = c()
  for (i in 1:nrow(plotdata)) {
    if (plotdata$level[i] %in% c("enriched","depleted") ) {
      c[i] = "*"
    }
    if (plotdata$level[i] == "nosig") {
      c[i] = ""
    }
  }
  plotdata$level3 = c
  plotdata$ID = factor(plotdata$ID,levels = unique(plotdata$ID)[length( unique(plotdata$ID)):1])
  
  p <- ggplot(plotdata) +　
    geom_bar(aes(y = ID,x = abundance,group = treat,fill = treat),stat = "identity",color  = "black",size = 0.5) +
    geom_vline(aes(xintercept=0), colour="black") + 
    geom_text(aes(y = ID,x = label_y,label = level3),color = "white") + 
    labs(title = "Control",y = "ASV of microbiome",
         x = "Abundance") + theme_bw()#  + scale_fill_manual(values = brewer.pal(9,"Set1"))
  p
  
  return(p)
}



tax_glom_wt <- function(ps = ps,ranks = "Phylum") {
  
  
  otu <- as.data.frame(t(vegan_otu(ps)))
  tax <- as.data.frame(vegan_tax(ps))
  
  # building group
  tax[[ranks]][is.na(tax[[ranks]])] = "Unknown"
  tax[[ranks]][tax[[ranks]] == ""] = "Unknown"
  split <- split(otu,tax[[ranks]])
  #calculate sum by group
  apply <- lapply(split,function(x)colSums(x[]))
  # result chack
  otucon <- do.call(rbind,apply)
  
  taxcon <- tax[1:match(ranks,colnames(tax))]
  taxcon <- taxcon[!duplicated(tax[[ranks]]),]
  #-tax name with NA wound be repeated with unknown
  taxcon[[ranks]][is.na(taxcon[[ranks]])] = "unknown"
  row.names(taxcon) <- taxcon[[ranks]]
  
  
  pscon <- phyloseq(
    otu_table( as.matrix(otucon),taxa_are_rows = TRUE),
    tax_table(as.matrix(taxcon)),
    sample_data(ps)
  )
  return(pscon)
}


