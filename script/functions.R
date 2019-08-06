# Paul Wackers
# 06-08-2019
# Functions used for the comparison approach as published in 
# Schaap et. al. 2012, Schaap et. al. 2015, Heusinkveld et. al. 2018, Luijten et. al. 2019 (submitted)
# R version 3.5.1 (2018-07-02)
# paul.wackers@rivm.nl
# +31 30 274 7512

# Function to read data
read.data.fun <- function(directory){
  fls <- dir(directory, pattern = ".txt", full.names = TRUE)
  print(paste("Found ", length(fls), " files", sep=""))
  files.list <- NULL
  for(iii in 1:length(fls)){
    compound <- as.vector(unlist(lapply(strsplit(x = basename(fls[iii]), split = "\\."), '[[', 1)))
    print(paste("Read: ", compound, sep=""))
    files.list[[compound]] <- read.delim(fls[iii], row.names=1)
  }
  return(files.list)
}

# Function to assign scores over compounds
score.over.compounds.fun <- function(TC1, TC2, n) {
  {
    # TC1 are the T-statistics of compound 1
    # TC2 are the T-statistics of compound 2
    # n is number of top genes to use in analysis
  }
  
  TC1.ids <- names(TC1)
  TC2.ids <- names(TC2)
  
  abs.T <- abs(c(TC1, TC2)) # make T statistics absolute
  o <- order(abs.T, decreasing = TRUE) # order
  linear_over <-
    c(seq(n, 1, by = -1), rep(0, (length(abs.T) - n))) # assign scores to the first n genes from n to 1 , rest 0
  names(linear_over) <- names(abs.T)[o] # assign names
  
  # Match TC1 to original order
  linear.TC1 <- linear_over[grep("_TC1", names(linear_over))]
  linear.TC1 <- linear.TC1[match(TC1.ids, names(linear.TC1))]
  # Match TC2 to original order
  linear.TC2 <- linear_over[grep("_TC2", names(linear_over))]
  linear.TC2 <- linear.TC2[match(TC2.ids, names(linear.TC2))]
  # Combine TC1 and TC2
  df <- as.data.frame(cbind(linear.TC1 = linear.TC1,
                            linear.TC2 = linear.TC2))
}

# Extract top n genes
top.genes.fun <- function(tstat, score, annotation, geneid, n) {
  {
    # tstat is T statistics
    # score is normalized score from score.over.compounds.fun
    # annotation is gene annotation
    # geneid is geneid columnname in annotation
    # n is number of top genes to use in analysis
  }
  top <- NULL
  names(tstat) <- annotation[, geneid]
  idx <- order(abs(tstat), decreasing = T)[1:n]
  top <- as.data.frame(cbind(
    geneid = annotation[idx, geneid],
    gene.idx = idx,
    sign = sign(tstat[idx]),
    tstat = tstat[idx],
    score = score[idx]
  ))
  rownames(top) <- top$GeneID
  return(top)
}

# Determine overlapping genes between the compounds
overlap.genes.fun <- function(TC1.top, TC2.top, annotation, id) {
  {
    # TC1.top is result of top.genes.fun of TC1
    # TC2.top is result of top.genes.fun of TC2
    # annotation is gene annotation
    # id is geneid columnname
  }
  
  # Intersection of genes in the same direction
  overlap.down <-
    intersect(TC1.top$geneid[which(TC1.top$sign == -1)], TC2.top$geneid[which(TC2.top$sign == -1)])
  overlap.up <-
    intersect(TC1.top$geneid[which(TC1.top$sign == 1)], TC2.top$geneid[which(TC2.top$sign == 1)])
  genes <- c(overlap.down, overlap.up)
  
  # Add scores
  score.over.compounds <-
    (as.numeric(TC1.top[match(genes, TC1.top$geneid), "score"]) + as.numeric(TC2.top[match(genes, TC2.top$geneid), "score"]))
  score.tstat <-
    (abs(as.numeric(TC1.top[match(genes, TC1.top$geneid), "tstat"])) + abs(as.numeric(TC2.top[match(genes, TC2.top$geneid), "tstat"])))
  result <-
    as.data.frame(
      cbind(
        annotation[match(genes, annotation[, id]),],
        score.over.compounds = score.over.compounds,
        score.tstat = score.tstat,
        regulation = c(rep("-", length(overlap.down)), rep("+", length(overlap.up)))
      )
    )
  
  return(result)
}

# Function to match two compounds
match.compounds.fun <- function(TC1, TC2, n, annotation, geneid) {
  {
    # TC1 are the T-statistics of compound 1
    # TC2 are the T-statistics of compound 2
    # n is number of top genes to use in analysis
    # annotation is gene annotation
    # geneid is geneid columnname in annotation
  }
  
  result <- NULL
  for (iii in 1:length(TC1)) {
    for (jjj in 1:length(TC2)) {
      # save match names
      match.name <-
        paste(colnames(TC1)[iii], "-", colnames(TC2)[jjj], sep = "")
      
      # Extract analysis data
      TC1.analysis <- TC1[, iii]
      #names(TC1.analysis) <- paste(rownames(TC1), "_TC1", iii, sep="")
      names(TC1.analysis) <- paste(rownames(TC1), "_TC1", sep = "")
      TC2.analysis <- TC2[, jjj]
      #names(TC2.analysis) <- paste(rownames(TC2), "_TC2", jjj, sep="")
      names(TC2.analysis) <- paste(rownames(TC2), "_TC2", sep = "")
      
      # call function for scoring over compounds scores (B in flowdiagram)
      score.over.compounds <-
        score.over.compounds.fun(TC1.analysis, TC2.analysis, n)
      
      # call function to extract top genes (A in flowdiagram)
      TC1.top.genes <-
        top.genes.fun(
          tstat = TC1[, iii],
          score = score.over.compounds$linear.TC1,
          annotation = annotation,
          geneid = geneid,
          n = n
        )
      TC2.top.genes <-
        top.genes.fun(
          tstat = TC2[, jjj],
          score = score.over.compounds$linear.TC2,
          annotation = annotation,
          geneid = geneid,
          n = n
        )
      
      # call function to determine overlap (C in flowdiagram)
      overlap.genes <-
        overlap.genes.fun(TC1.top = TC1.top.genes,
                          TC2.top = TC2.top.genes,
                          annotation,
                          geneid)
      result[[match.name]] <- overlap.genes
    }
  }
  return(result)
}

# Main function to analyse multiple compounds
comparison.approach.fun <-
  function(tstats, n=50, annotation, geneid, incl = FALSE) {
    {
      # tstats is T statistic list, a list where per compound the T-statistics are stored
      # n is number of top genes to use in analysis
      # annotation is annotation
      # geneid is geneid columnname in annotation file
      # incl: should the compound match it self?
    }
    
    result <- NULL
    compound.list <- NULL
    # Loop over T-statistic list
    for (iii in 1:length(tstats)) {
      TC1.compound <- names(tstats)[iii]
      compounds <- names(tstats[[iii]])
      compound.list <- as.data.frame(
        rbind(compound.list, cbind(compound=rep(TC1.compound, length(compounds)), compound.overall=compounds)))
      print(paste("Testing: ", TC1.compound, sep = ""))
      
      if (incl == TRUE) {
        TC2.data.idx <- seq(1, length(tstats))
      }
      if (incl == FALSE) {
        TC2.data.idx <- seq(1, length(tstats))[-iii]
      }
      importance.overall <- c()
      
      # Save Parameters
      values <-
        c(TC1.compound, paste(names(tstats)[TC2.data.idx], collapse = ";"), n, geneid)
      params <- c("Tested compound", "Database", "Topgenes", "GeneID")
      param.df <- as.data.frame(cbind(params, values))
      colnames(param.df) <- c("Parameters", "Values")
      
      comparison.list <- NULL
      for (jjj in TC2.data.idx) {
        TC2.compound <- names(tstats)[jjj]
        # call match.compounds.fun to match two compounds
        comparison <- match.compounds.fun(
          TC1 = tstats[[iii]],
          TC2 = tstats[[jjj]],
          n = n,
          annotation = annotation,
          geneid = geneid
        )
        comparison.list <- c(comparison.list, comparison)
        result[["Result"]][[TC1.compound]][["Parameters"]] <- param.df
      }
      
      # save results in results object
      result[["Compounds"]] <- compound.list
      result[["Result"]][[TC1.compound]][["Genes"]] <-
        comparison.list #
    }
    return(result)
  }

# Extract result tables
extract.result.fun <- function(comparison.object) {
  {
    # comparison.object is obcject returned from comparison.approach.fun
  }
   
    print("Extract results")
    
    result <- NULL
    
    compounds <- comparison.object$Compounds$compound.overall
    n.compounds <- length(compounds)
    
    # Create master tables
    {
      score.df <- matrix(NA, nrow = n.compounds, ncol = n.compounds)
      rownames(score.df) <- compounds
      colnames(score.df) <- compounds
      
      tstat.df <- matrix(NA, nrow = n.compounds, ncol = n.compounds)
      rownames(tstat.df) <- compounds
      colnames(tstat.df) <- compounds
      
      hits.df <- matrix(NA, nrow = n.compounds, ncol = n.compounds)
      rownames(hits.df) <- compounds
      colnames(hits.df) <- compounds
    }
    
    # Fill master tables
    for (i in 1:length(comparison.object$Result)) {
      tmp.matches <- comparison.object$Result[[i]]$Genes
      
      for (iii in 1:length(tmp.matches)) {
        compound <- names(tmp.matches[iii])
        yc <- sapply(strsplit(compound, "-"), function(x)
          x[1])
        xc <- sapply(strsplit(compound, "-"), function(x)
          x[2])
        idxy <- match(yc, compounds)
        idxx <- match(xc, compounds)
        score.df[idxy, idxx] <-
          sum(tmp.matches[[iii]]$score.over.compounds)
        tstat.df[idxy, idxx] <-
          round(sum(tmp.matches[[iii]]$score.tstat), digits = 0)
        hits.df[idxy, idxx] <- nrow(tmp.matches[[iii]])
      }
    }
    
    # create long list master table
    longlist.master <- NULL
    for(iii in 1:n.compounds){
      if(iii != n.compounds){
        label.tmp <- paste(compounds[iii], "-", compounds[(iii+1):n.compounds], sep="")
        score.tmp <- score.df[(iii+1):n.compounds, iii]
        tstat.tmp <- tstat.df[(iii+1):n.compounds, iii]
        hits.tmp <- hits.df[(iii+1):n.compounds, iii]
        tmp <- as.data.frame(cbind(Comparison = label.tmp, 
                                   Score=score.tmp, 
                                   TstatisticScore=tstat.tmp,
                                   NumberOfHits=hits.tmp))
        longlist.master <- as.data.frame(rbind(longlist.master, tmp))
      }
    }
    # remove NA's (i.e. within compound comparisons)
    remove <- apply(longlist.master[,2:4],1,function(x){all(is.na(x) == TRUE)})
    remove <- c(1:nrow(longlist.master))[remove]
    if(length(remove) > 0){
      longlist.master <- longlist.master[-remove,]
    }
    
    # Store master tables in result object
    result[["Master"]][["Score"]] <- score.df
    result[["Master"]][["TstatScore"]] <- tstat.df
    result[["Master"]][["NumberOfHits"]] <- hits.df
    result[["Master"]][["All"]] <- longlist.master
    
    compounds.overall <- unique(comparison.object$Compounds$compound)
    
    # Create best match tables
    {
      score.max <-
        matrix(NA,
               nrow = length(compounds.overall),
               ncol = length(compounds.overall))
      rownames(score.max) <- compounds.overall
      colnames(score.max) <- compounds.overall
      
      match.max <-
        matrix(NA,
               nrow = length(compounds.overall),
               ncol = length(compounds.overall))
      rownames(match.max) <- compounds.overall
      colnames(match.max) <- compounds.overall
      
      
      tstat.max <-
        matrix(NA,
               nrow = length(compounds.overall),
               ncol = length(compounds.overall))
      rownames(tstat.max) <- compounds.overall
      colnames(tstat.max) <- compounds.overall
      
      hits.max <-
        matrix(NA,
               nrow = length(compounds.overall),
               ncol = length(compounds.overall))
      rownames(hits.max) <- compounds.overall
      colnames(hits.max) <- compounds.overall
      
    }
    
    # Fill best match tables
     for (ccc in 1:length(compounds.overall)) {
      TC1.compound <- compounds.overall[ccc]
      TC2.compounds <- compounds.overall[-ccc]
      TC1.idx <- which(TC1.compound == comparison.object$Compounds$compound)

      for (mmm in 1:length(TC2.compounds)) {
        TC2.compound <- TC2.compounds[mmm]
        TC2.idx <- which(TC2.compound == comparison.object$Compounds$compound)
        
        score.max.compound <-
          max(score.df[TC2.idx, TC1.idx]) # max score
        score.max.idx <-  which(score.df[TC2.idx, TC1.idx] == max(score.df[TC2.idx, TC1.idx]), arr.ind = TRUE)
        TC1.score.max <- colnames(score.df[TC2.idx, TC1.idx])[score.max.idx[1,2]]
        TC2.score.max <- rownames(score.df[TC2.idx, TC1.idx])[score.max.idx[1,1]]
        score.max.match <- paste(TC1.score.max, TC2.score.max, sep="-")
        
        tstat.max.compound <-
          max(tstat.df[TC2.idx, TC1.idx]) # max tstats 

        hits.max.compound <- max(hits.df[TC2.idx, TC1.idx]) # max hits

        score.max[TC2.compound, TC1.compound] <- score.max.compound
        tstat.max[TC2.compound, TC1.compound] <- tstat.max.compound
        hits.max[TC2.compound, TC1.compound] <- hits.max.compound
        match.max[TC2.compound, TC1.compound] <- score.max.match
      }
    }
    
    # create long best match
    longlist.bestmatch <- NULL
    n.compounds.overall <- length(compounds.overall)
    for(iii in 1:n.compounds.overall){
      if(iii != n.compounds.overall){
        label.max.tmp <- paste(compounds.overall[iii], "-", compounds.overall[(iii+1):n.compounds.overall], sep="")
        match.max.tmp <- match.max[(iii+1):n.compounds.overall, iii]
        score.max.tmp <- score.max[(iii+1):n.compounds.overall, iii]
        tstat.max.tmp <- tstat.max[(iii+1):n.compounds.overall, iii]
        hits.max.tmp <- hits.max[(iii+1):n.compounds.overall, iii]
        tmp <- as.data.frame(cbind(Comparison = label.max.tmp,
                                   MatchName = match.max.tmp,
                                   Score=score.max.tmp, 
                                   TstatisticScore=tstat.max.tmp,
                                   NumberOfHits=hits.max.tmp))
        longlist.bestmatch <- as.data.frame(rbind(longlist.bestmatch, tmp))
      }
    }
    # remove NA's 
    remove <- apply(longlist.bestmatch[,2:4],1,function(x){all(is.na(x) == TRUE)})
    remove <- c(1:nrow(longlist.bestmatch))[remove]
    if(length(remove) > 0){
      longlist.bestmatch <- longlist.bestmatch[-remove,]
    }
    rownames(longlist.bestmatch) <- seq(1, nrow(longlist.bestmatch), by = 1)
    
    # Store best match tables in result object
    result[["BestMatch"]][["Score"]] <- score.max
    result[["BestMatch"]][["TstatScore"]] <- tstat.max
    result[["BestMatch"]][["NumberOfHits"]] <- hits.max
    result[["BestMatch"]][["All"]] <- longlist.bestmatch
    
    # Ordered tables
    score.max.o <- NULL
    tstat.max.o <- NULL
    hits.max.o <- NULL
    
    for (ccc in 1:length(compounds.overall)) {
      # score
      score.order <-
        order(score.max[, ccc], decreasing = T, na.last = TRUE)
      score.matches.o <-
        paste(rownames(score.max)[score.order], ": ", score.max[score.order, ccc], sep =
                "")
      score.max.o <-
        as.data.frame(cbind(score.max.o, score.matches.o))
      
      # tstat
      tstat.order <-
        order(tstat.max[, ccc], decreasing = T, na.last = TRUE)
      tstat.order.o <-
        paste(rownames(tstat.max)[tstat.order], ": ", tstat.max[tstat.order, ccc], sep =
                "")
      tstat.max.o <-
        as.data.frame(cbind(tstat.max.o, tstat.order.o))
      
      
      # number of hits
      hits.order <-
        order(hits.max[, ccc], decreasing = T, na.last = TRUE)
      hits.order.o <-
        paste(rownames(hits.max)[hits.order], ": ", hits.max[hits.order, ccc], sep =
                "")
      hits.max.o <- as.data.frame(cbind(hits.max.o, hits.order.o))
      
    }
    colnames(score.max.o) <- compounds.overall
    colnames(tstat.max.o)  <- compounds.overall
    colnames(hits.max.o) <- compounds.overall
    
    # Order longlist best match; score & hits
    order.longlist <- order(as.numeric(longlist.bestmatch$Score), as.numeric(longlist.bestmatch$NumberOfHits), decreasing = T) 
    longlist.bestmatch.o <- longlist.bestmatch[order.longlist, ]
    rownames(longlist.bestmatch.o) <- seq(1, nrow(longlist.bestmatch.o), by = 1)
    
    # Store ordered tables in result object
    result[["BestMatchOrdered"]][["Score"]] <- score.max.o
    result[["BestMatchOrdered"]][["TstatScore"]] <- tstat.max.o
    result[["BestMatchOrdered"]][["NumberOfHits"]] <- hits.max.o
    result[["BestMatchOrdered"]][["All"]] <- longlist.bestmatch.o
    
    return(result)
}

# Run analysis
run.analysis.fun <- function(tstats, n=50, annotation, geneid, incl = FALSE){
  {
    # tstats is T statistic list, a list where per compound the T-statistics are stored
    # n is number of top genes to use in analysis
    # annotation is annotation
    # geneid is geneid columnname in annotation file
    # incl: should the compound match it self?
  }
  
  analysis.result <- NULL
  
  # call comparison approach function
  comparison.result <- comparison.approach.fun(tstats=tstats, 
                                               n=n, 
                                               annotation=annotation, 
                                               geneid=geneid)
  
  analysis.result[["comparison.result"]] <- comparison.result
  
  # extract results
  result.tables <- extract.result.fun(comparison.object = comparison.result)
  
  analysis.result[["tables"]] <- result.tables
  
  return(analysis.result)
}

# function to extract hits from analysis object
getGenes.fun <- function(analysis.object, match.name){
  {
    # analysis.object is object returned from run.analysis.fun
    # match.name is compound match
  }
  
  match.idx <- which(analysis.object$tables$BestMatch$All$MatchName == match.name)
  TC1  <- as.vector(unlist(lapply(strsplit(analysis.object$tables$BestMatch$All$Comparison[match.idx], split = "-"), '[[', 1)))
  genelist <- analysis.object$comparison.result$Result[[TC1]]$Genes[[match.name]]
  return(genelist)
}

# function to extract result tables from analysis object
getResult.fun <- function(analysis.object, table, metric){
  {
    # analysis.object is object returned from run.analysis.fun
    # table is table to show, options:  Master (All vs All table)
    #                                   BestMatch (Compound vs Compound best match)
    #                                   BestMatchOrdered (Compound vs Compound best match, ordered by Score)
    # metric is score to show, options: Score (linear over compounds score)
    #                                   TstatScore (T-statistic score)
    #                                   NumberOfHits (Number of hits between compounds)
    #                                   All (All scores combined in a list)
  }
  
  return(analysis.object$tables[[table]][[metric]])
}
