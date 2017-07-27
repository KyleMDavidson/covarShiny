#ag=1,ac=2,at=3,aa=4
#ga=5,gc=6, gg=7, gt=8
#ca=9, cc=10, cg=11, ct=12
#ta=13, tc=14, tg=15, tt=16
#t or u decided by dC (diffChar parameter) passed to multiSeq, getExpDints

multiSeq=function(seqs,bp.count,pos1,pos2, dC = 'RNA'){
  length = length(seqs)
  if (missing(pos1)){
    pos1=1
  }
  if (missing(pos2)){
    pos2=length(seqs[[1]])
  }
  if (missing(bp.count)){
    bp.count = array(0, c(pos2,pos2,16))
  }

  chi.categories = array(0, c(pos2, pos2, 16))
  for (i in 1:length(seqs)){
    prebpcount= bpCount(seqs[[i]],seqs[[i]],pos1,pos2, dC)
    bp.count = bp.count + prebpcount[[1]]
    print(dim(prebpcount[[2]]))
    chi.categories = chi.categories + prebpcount[[2]]
  }
  chi.sum.categories = matrix(0,pos2,pos2)
  for (i in pos1:pos2){
    for (j in pos1:pos2){
      chi.sum.categories[i,j] = sum(chi.categories[i,j,]!=0)
      if (i==j){
        chi.sum.categories[i,j] = 0
        }
      }
    }
    
    bpcount.chi = list("bpcount"=bp.count, "chicat"=chi.sum.categories)
}


## windowedEnergySums2=function(bp.array, window){
##   folds.energy = vector('numeric', nrow(bp.array)-window)
##   for (i in 1:(nrow(bp.array)-window)){
##     folds.energy[i] = sum(bp.array[i:(i+window),i:(i+window)])
##   }
##   folds.energy
## }


## #energy.sums should be a list of energy sums
## compareWindowedEnergySums=function(sequences, window){
##   energy.sums = lapply(sequences, windowedEnergySums2, window = window)
##   intersect.values = Reduce(intersect, energy.sums)
##   }


#used to actually count dinucleotide pairs
bpCount=function(x,y,pos1,pos2, dC='RNA'){
  if (dC ==RNA){
    diffChar = 'u'
  }else{
    diffChar = 't'
  }
  if (missing(pos1)){
    pos1=1
  }
  if (missing(pos2)){
    pos2=length(x)
  }
 if (pos2-pos1 < 0){
      stop('negative range to sum energies across. please choose a positive window')
      }
  bp.count = array(0,c(pos2, pos2, 16))
  chi.Categories = array(0,c(pos2,pos2,16))
  for (i in pos1:pos2){
    for(j in i:pos2){
      if(x[i]=='a'){
        if(y[j]=='g'){
          chi.Categories[i,j,1] = 1
          bp.count[i,j,1]=bp.count[i,j,1]+1
        }else if(y[j]=='c'){
          chi.Categories[i,j,2] = 1
          bp.count[i,j,2]=bp.count[i,j,2]+1
        }else if(y[j]=='a'){
          chi.Categories[i,j,4] = 1
          bp.count[i,j,4]=bp.count[i,j,4]+1
        }else if (y[j]==diffChar){
          bp.count[i,j,3]=bp.count[i,j,3]+1
          chi.Categories[i,j,3] = 1
        }
      }else if(x[i]=='g'){
        if(y[j]=='a'){
          chi.Categories[i,j,5] = 1
          bp.count[i,j,5]=bp.count[i,j,5]+1
        }else if (y[j]=='c'){
          chi.Categories[i,j,6] = 1
          bp.count[i,j,6]=bp.count[i,j,6]+1
        }else if (y[j]=='g'){
          chi.Categories[i,j,7] = 1
          bp.count[i,j,7]=bp.count[i,j,7]+1
        }else if (y[j]==diffChar){
          bp.count[i,j,8]=bp.count[i,j,8]+1
          chi.Categories[i,j,8] = 1
        }
      }else if(x[i]=='c'){
        if(y[j]=='a'){
          chi.Categories[i,j,9] = 1
          bp.count[i,j,9]=bp.count[i,j,9]+1
        }else if (y[j]=='c'){
        chi.Categories[i,j,10] = 1
          bp.count[i,j,10]=bp.count[i,j,10]+1
        }else if (y[j]=='g'){
          chi.Categories[i,j,11] = 1
          bp.count[i,j,11]=bp.count[i,j,11]+1
        }else if (y[j]==diffChar){
          bp.count[i,j,12]=bp.count[i,j,12]+1
          chi.Categories[i,j,12] = 1
        }
      }else if(x[i]==diffChar){
        if(y[j]=='a'){
          chi.Categories[i,j,13] = 1
          bp.count[i,j,13]=bp.count[i,j,13]+1
        }else if (y[j]=='c'){
          chi.Categories[i,j,14] = 1
          bp.count[i,j,14]=bp.count[i,j,14]+1
        }else if (y[j]=='g'){
          chi.Categories[i,j,15] = 1
          bp.count[i,j,15]=bp.count[i,j,15]+1
        }else if (y[j]==diffChar){
          bp.count[i,j,16]=bp.count[i,j,16]+1
          chi.Categories[i,j,16] = 1
        }
      }
    }
  }

list("bpcount"=bp.count, "categories"=chi.Categories)
}


binder=function(a,b){
  tval = FALSE
  testv=c(a,b)
  if ('c'%in% testv && 'g' %in% testv && a!=b){
    tval = TRUE
  }else if ('g'%in% testv && a!=b&& 'u' %in% testv){
    tval=TRUE
  }else if ('g' %in% testv  && a!=b && 'c'%in% testv ){
    tval = TRUE
  }else if('g' %in% testv  && a!=b && 'u' %in% testv){
    tval=TRUE
  }else if('t' %in% testv && a!=b  && 'g'%in% testv ){
    tval = TRUE
  }else if('t' %in% testv && a!=b  && 'a' %in% testv){
    tval = TRUE
  }else if('a' %in% testv && a!=b && 't' %in% testv){
    tval = TRUE
  }
  tval
}


#returns the expected nucleotide frequency as derived from pairs of valid dinucleotides 
#between anywhere i and j in clustal output. Excludes any i or j part of i,j inclusive of NA chars like '-' or '.'
#---------additionally returns counts of dint pair participating nt at i,j
#same function as above but with t-->u
getExpNtFreq=function(seqs,pos1,pos2, dC = 'RNA'){
  if (dC == 'RNA'){
    diffChar = 'u'
  }else if (dC == 'DNA'){
    diffChar = 't'
  }
  if (missing(pos1)){
    pos1 = 1
  }
  if (missing(pos2)){
    pos2 = length(seqs[[1]])
  }
  if (dC == 'RNA'){
    alphabet = c('a','g','c','u')
  }else if (dC == 'DNA'){
    alphabet = c('a','g','c','t')
  }

  validcount=matrix(0,pos2,pos2)
  exp.ntfreq=array(0,c(pos2,pos2,8))
  for (i in pos1:pos2){
    for (j in i:pos2){
      count=0
      for (k in 1:length(seqs)){
        a= seqs[[k]][i]
        b= seqs[[k]][j]
        truthVectori= alphabet %in% a
        truthVectorj= alphabet %in% b
        if (TRUE %in% truthVectori){
          if (TRUE %in% truthVectorj){
            count=count+1
            if(a=='a'){
              exp.ntfreq[i,j,1]=exp.ntfreq[i,j,1]+1
            }else if(a=='g'){
              exp.ntfreq[i,j,2]=exp.ntfreq[i,j,2]+1
            }else if(a=='c'){
              exp.ntfreq[i,j,3]=exp.ntfreq[i,j,3]+1
            }else if(a==diffChar){
              exp.ntfreq[i,j,4]=exp.ntfreq[i,j,4]+1
            }
            if(b=='a'){
              exp.ntfreq[i,j,5]=exp.ntfreq[i,j,5]+1
            }else if(b=='g'){
              exp.ntfreq[i,j,6]=exp.ntfreq[i,j,6]+1
            }else if(b=='c'){
              exp.ntfreq[i,j,7]=exp.ntfreq[i,j,7]+1
            }else if(b==diffChar){
              exp.ntfreq[i,j,8]=exp.ntfreq[i,j,8]+1
            }
          }
        }
      }#end seqs for a given i,j
      validcount[i,j]=count
      if (count>0){
        exp.ntfreq[i,j,]=exp.ntfreq[i,j,]/count
      }
    }
  }
  #returns the expected nucleotide frequency as derived from pairs of valid dinucleotides 
  #between anywhere i and j in clustal output. Excludes any i or j part of i,j inclusive of NA chars like '-' or '.'
  #additionally returns counts of dint pair participating at i,j
  list('ExpNt'=exp.ntfreq,'validcount'=validcount)
}
#this uses individual sequence's dinucleotide output to calculate expectations for dinucleotides.


#same as above, but with getExpNtFreq--> getExpNtFreq
getExpDints=function(seqs,pos1,pos2, dC = 'RNA'){
  if (missing(pos1)){
    pos1 = 1
  }
  if (missing(pos2)){
    pos2=length(seqs[[1]])
    
  }
  exp.nt.freq.and.validcount=getExpNtFreq(seqs,pos1,pos2,dC)
  exp.ntfreq=exp.nt.freq.and.validcount[[1]]
  validcount=exp.nt.freq.and.validcount[[2]]
  length=length(seqs[[1]])
  exp.dints=array(0,c(pos2,pos2,16))
  for (i in 1:pos2){
    for (j in i:pos2){
      exp.ntfreq.ij=exp.ntfreq[i,j,]
      exp.dints[i,j,1]=exp.ntfreq.ij[1] * exp.ntfreq.ij[6]
      exp.dints[i,j,2]=exp.ntfreq.ij[1] * exp.ntfreq.ij[7]
      exp.dints[i,j,3]=exp.ntfreq.ij[1] * exp.ntfreq.ij[8]
      exp.dints[i,j,4]=exp.ntfreq.ij[1] * exp.ntfreq.ij[5]
      exp.dints[i,j,5]=exp.ntfreq.ij[2] * exp.ntfreq.ij[5]
      exp.dints[i,j,6]=exp.ntfreq.ij[2] * exp.ntfreq.ij[7]
      exp.dints[i,j,7]=exp.ntfreq.ij[2] * exp.ntfreq.ij[6]
      exp.dints[i,j,8]=exp.ntfreq.ij[2] * exp.ntfreq.ij[8]
      exp.dints[i,j,9]=exp.ntfreq.ij[3] * exp.ntfreq.ij[5]
      exp.dints[i,j,10]=exp.ntfreq.ij[3] * exp.ntfreq.ij[7]
      exp.dints[i,j,11]=exp.ntfreq.ij[3] * exp.ntfreq.ij[6]
      exp.dints[i,j,12]=exp.ntfreq.ij[3] * exp.ntfreq.ij[8]
      exp.dints[i,j,13]=exp.ntfreq.ij[4] * exp.ntfreq.ij[5]
      exp.dints[i,j,14]=exp.ntfreq.ij[4] * exp.ntfreq.ij[7]
      exp.dints[i,j,15]=exp.ntfreq.ij[4] * exp.ntfreq.ij[6]
      exp.dints[i,j,16]=exp.ntfreq.ij[4] * exp.ntfreq.ij[8]
      exp.dints[i,j,]=exp.dints[i,j,] * validcount[i,j]
    }
  }
  
  exp.dints
}


#this deploys the packages methods in order to get
# for each i,j positions in sequence alignments 
# the sum(from 1-16) of (observed # of dinucleotides - expected # of dinucleotides)^2 / expected # dinucleotides
#--------returns a two dimensional matrix with chi square values.
#this takes  either 1. specification of filename by file.choose OR a manual entry passed as args
# this does not name rows and columncs for csv/excel output
covarRNA=function(filepath,filename,pos1,pos2, dC = 'RNA'){
  if (missing(filename)){
  f = file.choose()
  seqs = coveIn(f)
  
  }else{
    seqs = coveIn(filepath,filename)
  }
  if (missing(pos1)){
    pos1 = 1
  }
  if (missing(pos2)){
    pos2 = length(seqs[[1]])
  }
  
  exp.dints = getExpDints(seqs,pos1,pos2, dC)
  bp.count = array(0, c(pos2,pos2,16))
  dints.chi = multiSeq(seqs,bp.count,pos1,pos2, dC)
  obs.dints = dints.chi[[1]]
  chi.categories = dints.chi[[2]]


  bp.chi = array(0,c(pos2,pos2,16))
  bp.chi = obs.dints - exp.dints
  bp.chi = bp.chi * bp.chi
  bp.chi = bp.chi / exp.dints
  bp.chi[is.na(bp.chi)] = 0
  
  covar.grid = array(0,c(pos2,pos2))
  
  
  for(i in pos1:pos2){
    for(j in i:pos2){
        x = sum(bp.chi[i,j,])
        covar.grid[i,j] = pchisq(x, chi.categories[i,j])
    }
  }
  for (i in pos1:pos2){
    covar.grid[i,i] = 0
  }

  covar.grid
}


#meant for use with ui that determines file(s) to be passed as arguments
#covarShiny package calls it from interface
covarRNAui=function(f,pos1,pos2,dC='RNA'){
  seqs= coveIn(f)
  length = length(seqs[[1]])
  if (missing(pos1)){
    pos1=1
  }
  if (missing(pos2)){
    pos2=length
  }
  if (is.na(pos1)){
    pos1=1
  }
  if (is.na(pos2)){
    pos2=length
  }
  if (pos2 > length){
    pos2=length
  }
  
  exp.dints = getExpDints(seqs,pos1,pos2,dC)
  bp.count = array(0, c(pos2,pos2,16))
  dints.chi = multiSeq(seqs,bp.count,pos1,pos2,dC)
  ## obs.dints holds observed dinucleotide pair counts in 3 d table
  obs.dints = dints.chi[[1]]
  ##chi.categories holds number of types of dinucleotide pairs found at i,j in 2 d table
  chi.categories = dints.chi[[2]]
  
  
  bp.chi = array(0,c(pos2,pos2,16))
  bp.chi = obs.dints - exp.dints
  bp.chi = bp.chi * bp.chi
  bp.chi = bp.chi / exp.dints
  bp.chi[is.na(bp.chi)] = 0
  
  covar.grid = array(0, c(pos2, pos2))

  
  
  
  for(i in pos1:pos2){
    for(j in i:pos2){
      covar.grid[i,j] = sum(bp.chi[i,j,])
    }
  }
  for (i in pos1:pos2){
    covar.grid[i,i] = 0
  }

    pval.grid = matrix(0, pos2,pos2)
  for (i in 1:16){
    j = which(chi.categories==i)
    trans = chi.categories[j]
    transcovar = covar.grid[j]
    pvals = vector("numeric", length=length(trans))
    for ( h in 1:length(trans)){
      pvals[h] = 1-pchisq(transcovar[h],trans[h])
      }
    pval.grid[j] = pvals
    }

  g = covar.grid[pos1:pos2,pos1:pos2]
  row.names(g) = pos1:pos2
  colnames(g) = pos1:pos2
  list("covar"=g, "pvals"=pval.grid)
 }

  ## degrees of freedom are n-1 where n is types of dinucleotide pairs occurring at i,j.




#filepatht="../../Documents/Syllabi/Grabowskilab/C1_alignments/"
#filenamet= "tRNA.fasta.txt"
