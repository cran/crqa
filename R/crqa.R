## crqa, adapted from a Matlab code developed at
## Cincinnati school (we need full ref. here)
## written in R by Moreno I. Coco, 2013, (moreno.cocoi@gmail.com)

## arguments to pass to crqa:
## ts1, ts2: times series of integers indicating the states
## delay = nr. of lags
## embed = the embedding dimension, i.e., the lag intervals
## rescale = the normalization for the distance; 
##           if 1 (Mean Distance); if 2 (Max Distance)
## radius = the maximum distance to be calculated (set it very 
##           small, if the series are categorical in nature
## normalize = rescale factor for source data; 
##           if 1 (Unit interval); if 2 (z-score) 
## mindiagline = set a minimum diagonal line length
## mindiagline = set a minimum vertical line length

## delay = 1; embed = 1; rescale = 1; radius = 0.00001;
## normalize = 0; mindiagline = 2; minvertline = 2;
## whiteline = TRUE - flag to compute or not white vertical lines
##                    in the recurrence plot. Note, white lines are not
##                    yet used to derive any particular measure
## recpt = TRUE - flag to indicate whether the input ts1 is already
##                a recurrence plot

packageName <- 'crqa'

crqa <- function(ts1, ts2, delay, embed, rescale,
                 radius, normalize, mindiagline, minvertline,
                 whiteline, recpt){
    
    v11 = v21 = NULL ## stupid initializations to please CRAN
    
   ## require("fields") ## to compute the Euclidean distance matrix
   ## require("Matrix")  ## to manipulate sparse matrices 


    ## check if the input is a recurrence plot 
    if (recpt == FALSE){

    
        ts1 = as.vector(as.matrix(ts1)) ## make sure data is a vector
        ts2 = as.vector(as.matrix(ts2))
    
        if (is.matrix(ts1)){ stop("Your data must consist of a single column of data.")}  
        if (is.matrix(ts2)){ stop("Your data must consist of a single column of data.")}
    
    
    ##chop of sequences if they are off different lengths
    
        if (length(ts1) != length(ts2)){
            shortest = min(c(length(ts1), length(ts2)) );
            ts1 = ts1[1:shortest];
            ts2 = ts2[1:shortest];
        }
    
    ##rescale the data if really necessary
        
        if (normalize > 0){
            switch (normalize,
                    {1 
                     ts1 = (ts1 - min(ts1));
                     ts1 = ts1 / max(ts1);
                     ts2 = (ts2 - min(ts2));
                     ts2 = ts2 / max(ts2);},
                    
                    {2                      
                     ts1 = (ts1 - mean(ts1))/sd(ts1) # zscore(ts1, na.rm = TRUE, robust = FALSE); ## using R.basics
                     ts2 = (ts2 - mean(ts2))/sd(ts2) # zscore(ts2, na.rm = TRUE, robust = FALSE)
                 }
                    )
        }
        
        ## start to compute recurrence
        ## do it twice for the two series
        
        for (loop in 1:embed){
            vectorstart = (loop-1) * delay + 1;
            vectorend = length(ts1) - ( (embed-loop)* delay);
            assign(paste("v1", loop, sep =""), ts1[vectorstart:vectorend]);
        }
        
        for (loop in 1:embed){
            vectorstart = (loop-1) * delay + 1;
            vectorend = length(ts2) - ( (embed-loop)* delay);
            assign(paste("v2", loop, sep =""), ts2[vectorstart:vectorend]);
        }
        
        ## Create matrix from vectors to use for distance matrix calcs

        dimts1 = dimts2 = vector() ## initialize dims of embedding 
        
        for (loop in 1:embed){
            if (loop == 1){ dimts1 = v11 }
            else{
                eval(
                    parse(
                        text = paste("dimts1 = cbind(dimts1,",
                            paste( "v1", loop, sep = ""),
                            ", deparse.level = 0)", sep = "")
                        )
                    )
            }
        }
    
        for (loop in 1:embed){
            if (loop == 1){ dimts2 = v21 }
            else{
                eval(
                    parse(
                        text = paste("dimts2 = cbind(dimts2,",
                            paste( "v2", loop, sep = ""),
                            ", deparse.level = 0)", sep = "")
                        )
                    )
            }
        }
        
        
        ## Compute euclidean distance matrix
        vlength = length(v11) ## just to have the length of matrix saved
        dm = rdist(dimts1,dimts2);
        
        ## Find indeces of the distance matrix that fall
        ## within prescribed radius.
        
        switch(rescale,
               {1  ## Create a distance matrix that is re-scaled
                   ## to the mean distance
                
                rescaledist = mean( dm )    
                dmrescale=(dm/rescaledist)*100},
               
               {2  ## Create a distance matrix that is re-scaled
                   ## to the max distance
                
                rescaledist = max(dm);
                dmrescale = (dm/rescaledist)*100}
           )
        
        ## Compute recurrence matrix
        
        ind = which(dmrescale <= radius, arr.ind=T);
        r = ind[,1]; c = ind[,2]

    } else { ## take as input an RP directly

        ind = which(ts1 > 0, arr.ind = T)
        ## just a trick to reduce the number of lines
        ## of the code
        r = ind[,1]; c = ind[,2]
        
    }
    
    if (length(r) != 0 & length(c) != 0){ ##avoid cases with no recurrence
        S = sparseMatrix(r, c) ## this is the recurrent plot
  
        spdiagonalize = spdiags(S) ##  spdiags should have decent speed 
        B = spdiagonalize$B
        
        ##calculate percentage recurrence by taking all non-zeros
        
        numrecurs = length(which(B == TRUE));
        percentrecurs = (numrecurs/((vlength^2)))*100;
        
####################################################################
####################################################################
        
        ## Computing the line counts

        ## This section finds the index of the zeros in the matrix B,
        ## which contains the diagonals of one triangle of the
        ## recurrence matrix (the identity line excluded).

        ## The find command indexes the matrix sequentially
        ## from 1 to the total number of elements.
        ## The element numbers for a 2X2 matrix would be [1 3; 2 4].
        ## You get a hit for every zero. If you take the difference
        ## of the resulting vector, minus 1, it yields the length of an
        ## interceding vector of ones, a line. Here is an e.g.
        ## using a row vector rather than a col. vector, since it types
        ## easier: B=[0 1 1 1 0], a line of length 3.
        ## find( B == 0 ) yields [1 5], diff( [1 5] ) -1 = 3,
        ## the line length.
        ## So this solution finds line lengths in the interior of
        ## the B matrix, BUT fails if a line butts up against either
        ## edge of the B matrix, e.g. say  B = [0 1 1 1 1],
        ## which( B == 0) returns a 1, and you miss the line of length 4.
        ## A solution is to "bracket" B with a row of zeros at each
        ## top and bottom.

        ## Bracket B with zeros
        if (is.vector(B)) {
            false = rep(FALSE, length(B)) ##cases where B is a vector
            B = rbind(false, B, false, deparse.level = 0)
        } else  {
            false = rep(FALSE, ncol(B))
            B = as.matrix(B)
            ## need to transform the sparseMat into normal to bracket it
            B = rbind(false, B, false, deparse.level = 0)
        }
        
        ## Get list of line lengths, sorted from largest to smallest
        diaglines = sort( diff(which(B == FALSE) ) -1, decreasing = TRUE)
        
        ## Delete line counts less than the minimum diagonal.
        diaglines = diaglines[-which(diaglines < mindiagline)]
        ## diaglines(diaglines>200)=[]; # Can define a maximum line length too.

        ## exlude the rare cases where there are no diaglines
        
        if(length(diaglines) != 0){
            
            numdiaglines = length(diaglines) ## extract the length of diag
            maxline = max(diaglines)
            meanline = mean(diaglines)
            
            tabled = as.data.frame(table(diaglines))
            
            total = sum(tabled$Freq)       
            p = tabled$Freq/total
            
            ##remove zero probability..it should not be necessary
            del = which(p == 0 )
            if (length(del) > 0) {
                p = p[-del]
            }
            
            ## entropy log2, and relative entropy divided by max
            entropy = -1 * sum(p*log2(p))    
            relEntropy = entropy/(-1*log2(1/nrow(tabled)))

            ## entropy/max entropy: comparable across contexts and conditions.
        
            pdeter = (sum(diaglines)/numrecurs)*100
            ## percent determinism: the predictability of the dynamical system 
        
            ## calculate laminarity and trapping time
            restt = tt(S, minvertline, whiteline)
            lam = restt$lam; TT = restt$TT
                    
        } else {
            
            numdiaglines = 0; maxline = 0; pdeter = NA;
            entropy = NA; relEntropy = NA; meanline = 0
            lam = 0; TT = 0; RP = NA
        }
        
        results = list(rec = percentrecurs, det = pdeter, 
            nrline = numdiaglines, maxline = maxline, 
            meanline = meanline, entropy = entropy, 
            relEntropy = relEntropy,
            lam = lam, TT = TT, RP = S)
        
    } else { print (paste ("No recurrence found") )
             results = list(rec = 0, det = NA, nrline = 0, maxline = 0, 
                 meanline = 0, entropy = NA, relEntropy = NA,
                 lam = NA, TT = NA, RP = NA)}  
    
    return (results)
    
}
