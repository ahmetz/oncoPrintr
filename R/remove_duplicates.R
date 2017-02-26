 remove_duplicates <- function(df){
    cat("Filtering duplicates\n")
    choose_alteration_type <- function(x){
      for (i in 1:length(x)){
        if(grepl("Frameshift", x[i])){
          return(x[grep("Frameshift", x[i])])
        }else if(grepl("InFrame", x[i])){
          return(x[grep("InFrame", x[i])])
        }else if(grepl("Missense", x[i])){
          return(x[grep("Missense", x[i])])
        }else if(grepl("Mutation", x[i])){
          return(x[grep("Mutation", x[i])])
        }else if(grepl("Nonsense", x[i])){
          return(x[grep("Nonsense", x[i])])
        }else if(grepl("Promoter", x[i])){
          return(x[grep("Promoter", x[i])])
        }else if(grepl("Splicing", x[i])){
          return(x[grep("Splicing", x[i])])
        }else if(grepl("stopgain", x[i])){
          return(x[grep("stopgain", x[i])])
        }else if(grepl("insertion", x[i])){
          return(x[grep("insertion", x[i])])
        }else if(grepl("deletion", x[i])){
          return(x[grep("deletion", x[i]) ])
        }else if(grepl("splicing", x[i])){
          return(x[grep("splicing", x[i]) ])
        }else if(grepl("nonsynonymous", x[i])){
          return(x[grep("nonsynonymous", x[i]) ])
        }else if(grepl("upstream", x[i])){
          return(x[grep("upstream", x[i]) ])
        }else if(grepl("Present", x[i])){
          return(x[grep("Present", x[i]) ])
        }else if(grepl("NotTested", x[i])){
          return(x[grep("NotTested", x[i]) ])
        }else if(grepl("del", x[i])){
          return(x[grep("del", x[i]) ])
        }else if(grepl("homo-del", x[i])){
          return(x[grep("homo-del", x[i]) ])
        }else if(grepl("CN-LOH", x[i])){
          return(x[grep("CN-LOH", x[i]) ])
        }else if(grepl("LOH", x[i])){
          return(x[grep("LOH", x[i]) ])
        }else if(grepl("None", x[i])){
          return(x[grep("None", x[i]) ])
        }else if(grepl("Fusion", x[i])){
          return(x[grep("Fusion", x[i]) ])
        }else if(grepl("Amplification", x[i])){
          return(x[grep("Amplification", x[i]) ])
        }else if(grepl("Deletion", x[i])){
          return(x[grep("Deletion", x[i]) ])
        }else if(grepl("Yes", x[i])){
          return(x[grep("Yes", x[i]) ])
        }else if(grepl("No", x[i])){
          return(x[grep("No", x[i]) ])
        }else if(grepl("Pathogenic", x[i])){
          return(x[grep("Pathogenic", x[i]) ])
        }else if(grepl("IntragenicDeletion", x[i])){
          return(x[grep("IntragenicDeletion", x[i]) ])
        }else{
          return(x)
        }
      }
    }
    df <- df %>% group_by(Gene, Sample) %>% summarise(VarClass = choose_alteration_type(VarClass)) %>% ungroup()
    #df %>% group_by(Gene, Sample) %>% unique() %>% ungroup()
    cat("Finished filtering duplicates\n")
    return(df)
  }
