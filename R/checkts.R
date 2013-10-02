## arguments: ts1, ts2, the two time-series,
## datatype: (numerical, categorical)
## 

checkts <- function(ts1, ts2, datatype, thrshd){

  ## take as input the two sequences and normalize them by a pre-defined
  ## threshold constant.

  ln1 = length(ts1)
  ln2 = length(ts2)
  
  
  if (datatype == "numeric"){ ## just make sure that the vec are numeric
    ts1 = as.numeric(as.matrix( ts1) )
    ts2 = as.numeric(as.matrix( ts2) )
  }
  
  
  if(datatype == "categorical"){ ## recode sp into numerical vec.
    
    ts1 = as.character( as.matrix(ts1) )
    ts2 = as.character( as.matrix(ts2) )
    
    sprecoded = numerify(ts1, ts2) ## assign numerical indeces to categorical factors. 
    ts1 = sprecoded[[1]]; ts2 = sprecoded[[2]]
    
  }
  
  if ( ln1 != ln2 ){    ## timeseries can differ in length
    df = abs(ln1 - ln2)
    
    ## TODO if length of time-series is different there is a threshold
    ## find optimal solution to accepted threshold
    
    if (df <= thrshd){     ## threshold to adjust | remove pair
      
      if (ln2 > ln1){ ts2 = ts2[1:(length(ts2)-df)] } ## remove final time-points to series of equal length
      if (ln1 > ln2){ ts1 = ts1[1:(length(ts1)-df)] }
      
    }         
  }
  
  if ( length(ts1) == length(ts2) ){
    
    return (list(cbind(ts1, ts2), TRUE) )
    
  } else {
    
    return (list(df, FALSE) )## how many units are the sequences differing. 
            
  }
  
}


numerify <- function( ts1, ts2 ){

  ## transform categorical series into numerical series
  ## arguments: two categorical sequences 

  nwts1 = nwts2  = vector()
  
  objects = sort(unique(c(ts1, ts2 )))
  ids = seq(1, length(objects), 1)
  
  for ( o in 1:length(objects) ){
    nwts1[which(ts1 == objects[o])] = ids[o]
    nwts2 [which(ts2  == objects[o])] = ids[o]
  }
  
  return( list( nwts1, nwts2  ) )
  
}
