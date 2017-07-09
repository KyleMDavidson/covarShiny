#returns DNA sequence i,j pairs with the following paradigm for the third dimension table specifying type of bp:
#ag=1,ac=2,at=3,aa=4
#ga=5,gc=6, gg=7, gt=8
#ca=9, cc=10, cg=11, ct=12
#ta=13, tc=14, tg=15, tt=16
multiSeq=function(seqs,bp.count,pos1,pos2){
  length = length(seqs)
  if (missing(pos1)){
    pos1=1
  }
  if (missing(pos2)){
    pos2=length(seqs[[1]])
  }

  for (i in 1:length(seqs)){
    bp.count = bpCount(seqs[[i]],seqs[[i]],bp.count,pos1,pos2)
  }
  bp.count
}

#returns for RNA (t->u)
multiSeqRNA=function(seqs,bp.count,pos1,pos2){
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
  chi.categories = array(0, c(pos2, pos2, 8))
  for (i in 1:length(seqs)){
    prebpcount= bpCountRNA(seqs[[i]],seqs[[i]],pos1,pos2)
    bp.count = bp.count + prebpcount[[1]]
    chi.categories = chi.categories + prebpcount[[2]]
  }
  chi.product.categories = matrix(0,pos2,pos2)
  for (i in pos1:pos2){
    for (j in pos1:pos2){
      chitemp = which(chi.categories[i,j,1:4] != 0)
      l1 = length(chitemp)
      chitemp = which(chi.categories[i,j,5:8] != 0)
      l2 = length(chitemp)
      n = l1*l2 -1
      chi.product.categories[i,j] = n
      }
    }
  bpcount.chi = list("bpcount"=bp.count, "chicat"=chi.product.categories)
}
#calculate the bp.array with bpCount beforehand

windowedEnergySums1=function(bp.array, window){
  for (i in 1:(nrow(bp.array)-window)){
    b <- sum(bp.array[i:(i+window),i:(i+window)])
    a <- as.character(i)
    assign(a, b)
   
  }
    names = list(as.character(1:(nrow(bp.array)-window)))
    print(names)
    energy.sums = lapply(names, get)
  
}


windowedEnergySums2=function(bp.array, window){
  folds.energy = vector('numeric', nrow(bp.array)-window)
  for (i in 1:(nrow(bp.array)-window)){
    folds.energy[i] = sum(bp.array[i:(i+window),i:(i+window)])
  }
  folds.energy
}


#energy.sums should be a list of energy sums
compareWindowedEnergySums=function(sequences, window){
  energy.sums = lapply(sequences, windowedEnergySums2, window = window)
  intersect.values = Reduce(intersect, energy.sums)
  }


bpCount=function(x,y,bp.count,pos1,pos2){
  if (missing(pos1)){
    pos1=1
  }
  if (missing(pos2)){
    pos2=length(x)
  }
  if (pos2-pos1 < 0){
      stop('negative range to sum energies across. please choose a positive window')
      }
  bp.count = array(0,c(pos2-pos1, pos2-pos1, 16))
  for (i in pos1:pos2){
    for(j in i:pos2){
      if(x[i]=='a'){
        if(y[j]=='g'){
          bp.count[i,j,1]=bp.count[i,j,1]+1
        }
        if(y[j]=='c'){
          bp.count[i,j,2]=bp.count[i,j,2]+1
        }
        if(y[j]=='a'){
          bp.count[i,j,4]=bp.count[i,j,4]+1
        }
        if(y[j]=='t'){
          bp.count[i,j,3]=bp.count[i,j,3]+1
        }
      }else if(x[i]=='g'){
        if(y[j]=='a'){
          bp.count[i,j,5]=bp.count[i,j,5]+1
        }
        if(y[j]=='c'){
          bp.count[i,j,6]=bp.count[i,j,6]+1
        }
        if(y[j]=='g'){
          bp.count[i,j,7]=bp.count[i,j,7]+1
        }
        if(y[j]=='t'){
          bp.count[i,j,8]=bp.count[i,j,8]+1
        }
      }else if(x[i]=='c'){
        if(y[j]=='a'){
          bp.count[i,j,9]=bp.count[i,j,9]+1
        }
        if(y[j]=='c'){
          bp.count[i,j,10]=bp.count[i,j,10]+1
        }
        if(y[j]=='g'){
          bp.count[i,j,11]=bp.count[i,j,11]+1
        }
        if(y[j]=='t'){
          bp.count[i,j,12]=bp.count[i,j,12]+1
        }
      }else if(x[i]=='t'){
        if(y[j]=='a'){
          bp.count[i,j,13]=bp.count[i,j,13]+1
        }
        if(y[j]=='c'){
          bp.count[i,j,14]=bp.count[i,j,14]+1
        }
        if(y[j]=='g'){
          bp.count[i,j,15]=bp.count[i,j,15]+1
        }
        if(y[j]=='t'){
          bp.count[i,j,16]=bp.count[i,j,16]+1
        }
      }
    }
  }
  bp.count
}

#used to actually count dinucleotide pairs
bpCountRNA=function(x,y,pos1,pos2){
  if (missing(pos1)){
    pos1=1
  }
  if (missing(pos2)){
    pos2=length(x)
  }
 if (pos2-pos1 < 0){
      stop('negative range to sum energies across. please choose a positive window')
      }
  bp.count = array(0,c(pos2-pos1+1, pos2-pos1+1, 16))
  chi.Categories = array(0,c(pos2, pos2,8))
  chi.sumCategories = matrix(0,pos2,pos2)
  for (i in pos1:pos2){
    for(j in i:pos2){
      if(x[i]=='a'){
        chi.Categories[i,j,1] = 1
        if(y[j]=='g'){
          chi.Categories[i,j,6] = 1
          bp.count[i,j,1]=bp.count[i,j,1]+1
        }else if(y[j]=='c'){
          chi.Categories[i,j,7] = 1
          bp.count[i,j,2]=bp.count[i,j,2]+1
        }else if(y[j]=='a'){
          chi.Categories[i,j,5] = 1
          bp.count[i,j,4]=bp.count[i,j,4]+1
        }else if (y[j]=='u'){
          bp.count[i,j,3]=bp.count[i,j,3]+1
          chi.Categories[i,j,8] = 1
        }
      }else if(x[i]=='g'){
          chi.Categories[i,j,2] = 1
        if(y[j]=='a'){
          chi.Categories[i,j,5] = 1
          bp.count[i,j,5]=bp.count[i,j,5]+1
        }else if (y[j]=='c'){
        chi.Categories[i,j,7] = 1
          bp.count[i,j,6]=bp.count[i,j,6]+1
        }else if (y[j]=='g'){
        chi.Categories[i,j,6] = 1
          bp.count[i,j,7]=bp.count[i,j,7]+1
        }else if (y[j]=='u'){
          bp.count[i,j,8]=bp.count[i,j,8]+1
          chi.Categories[i,j,8] = 1
        }
      }else if(x[i]=='c'){
        chi.Categories[i,j,3] = 1
        if(y[j]=='a'){
          chi.Categories[i,j,5] = 1
          bp.count[i,j,9]=bp.count[i,j,9]+1
        }else if (y[j]=='c'){
        chi.Categories[7] = 1
          bp.count[i,j,10]=bp.count[i,j,10]+1
        }else if (y[j]=='g'){
          chi.Categories[i,j,6] = 1
          bp.count[i,j,11]=bp.count[i,j,11]+1
        }else if (y[j]=='u'){
          bp.count[i,j,12]=bp.count[i,j,12]+1
          chi.Categories[i,j,8] = 1
        }
      }else if(x[i]=='u'){
          chi.Categories[i,j,4] = 1
        if(y[j]=='a'){
          bp.count[i,j,13]=bp.count[i,j,13]+1
        }else if (y[j]=='c'){
        chi.Categories[i,j,7] = 1
          bp.count[i,j,14]=bp.count[i,j,14]+1
        }else if (y[j]=='g'){
          chi.Categories[i,j,6] = 1
          bp.count[i,j,15]=bp.count[i,j,15]+1
        }else if (y[j]=='u'){
          bp.count[i,j,16]=bp.count[i,j,16]+1
          chi.Categories[i,j,8] = 1
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


getBasepairs=function(x,y,pos1,pos2){
  if (missing(pos1)){
  pos1 = 1
  } 
  if (missing(pos2)){
  pos2 = length(x)
  }
  bindCounts = matrix(0,nrow=length(x),ncol=length(y))
  for (i in pos1:pos2){
    for(j in i:pos2){
      if (binder(x[i],y[j])){
        bindCounts[i,j] = bindCounts[i,j]+1
      }
    }
  }
  bindCounts
}


#returns the expected nucleotide frequency as derived from pairs of valid dinucleotides 
#between anywhere i and j in clustal output. Excludes any i or j part of i,j inclusive of NA chars like '-' or '.'
#---------additionally returns counts of dint pair participating nt at i,j
getExpNtFreq=function(seqs,pos1,pos2){
  if (missing(pos1)){
    pos1 = 1
  }
  if (missing(pos2)){
    pos2 = length(seqs[[1]])
  }
  alphabet = c('a','g','c','t')
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
              }else if(a=='t'){
                exp.ntfreq[i,j,4]=exp.ntfreq[i,j,4]+1
              }
              if(b=='a'){
                exp.ntfreq[i,j,5]=exp.ntfreq[i,j,5]+1
              }else if(b=='g'){
                exp.ntfreq[i,j,6]=exp.ntfreq[i,j,6]+1
              }else if(b=='c'){
                exp.ntfreq[i,j,7]=exp.ntfreq[i,j,7]+1
              }else if(b=='t'){
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
 
  list('ExpNt'=exp.ntfreq,'validcount'=validcount)
}
#same function as above but with t-->u
getExpNtFreqRNA=function(seqs,pos1,pos2){
  if (missing(pos1)){
    pos1 = 1
  }
  if (missing(pos2)){
    pos2 = length(seqs[[1]])
  }
  alphabet = c('a','g','c','u')
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
            }else if(a=='u'){
              exp.ntfreq[i,j,4]=exp.ntfreq[i,j,4]+1
            }
            if(b=='a'){
              exp.ntfreq[i,j,5]=exp.ntfreq[i,j,5]+1
            }else if(b=='g'){
              exp.ntfreq[i,j,6]=exp.ntfreq[i,j,6]+1
            }else if(b=='c'){
              exp.ntfreq[i,j,7]=exp.ntfreq[i,j,7]+1
            }else if(b=='u'){
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
getExpDints=function(seqs,pos1,pos2){
  if (missing(pos1)){
    pos1 = 1
  }
  if (missing(pos2)){
    pos2=length(seqs[[1]])
  }
exp.nt.freq.and.validcount=getExpNtFreq(seqs,pos1,pos2)
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



#same as above, but with getExpNtFreq--> getExpNtFreqRNA
getExpDintsRNA=function(seqs,pos1,pos2){
  if (missing(pos1)){
    pos1 = 1
  }
  if (missing(pos2)){
    pos2=length(seqs[[1]])
    
  }
  exp.nt.freq.and.validcount=getExpNtFreqRNA(seqs,pos1,pos2)
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
covar=function(filepath,filename,pos1,pos2){
  seqs= coveIn(filepath,filename)
  if (missing(pos1)){
    pos1=1
  }
  if (missing(pos2)){
    pos2=length(seqs[[1]])
  }
  
  exp.dints=getExpDints(seqs, pos1, pos2)
  bp.count=array(0, c(pos2, pos2, 16))
  obs.dints=multiSeq(seqs, bp.count, pos1, pos2)
  
  bp.chi = array(0, c(pos2, pos2, 16))
  bp.chi = obs.dints - exp.dints
  bp.chi = bp.chi * bp.chi
  bp.chi = bp.chi / exp.dints
  bp.chi[is.na(bp.chi)] = 0
  
  covar.grid = array(0,c(pos2,pos2))
  
  
  for(i in pos1:pos2){
    for(j in i:pos2){
      covar.grid[i,j] = sum(bp.chi[i,j,])
    }
  }
  for (i in pos1:pos2){
    covar.grid[i,i] = 0
  }
  
  covar.grid
}

#same as above but with the method calls to RNA-friendly character mapping ('u') instead of DNA ('t')
#this take either 1. specification of filename by file.choose OR a manual entry passed as args
# this does not name rows and columncs for csv/excel output
covarRNA=function(filepath,filename,pos1,pos2){
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
  
  exp.dints = getExpDintsRNA(seqs,pos1,pos2)
  bp.count = array(0, c(pos2,pos2,16))
  dints.chi = multiSeqRNA(seqs,bp.count,pos1,pos2)
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


#same as above but with the method calls to RNA-friendly character mapping ('u') instead of DNA ('t')
#this takes exclusively ui entered file for analysis
# this does not name rows and columncs for csv/excel output
#meant for use with ui that determines file(s) to be passed as arguments
covarRNAui=function(f,pos1,pos2){
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
  
  exp.dints = getExpDintsRNA(seqs,pos1,pos2)
  bp.count = array(0, c(pos2,pos2,16))
  obs.dints = multiSeqRNA(seqs,bp.count,pos1,pos2)
  
  
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
  
  g = covar.grid[pos1:pos2,pos1:pos2]
  row.names(g) = pos1:pos2
  colnames(g) = pos1:pos2
  g
 }

  ## degrees of freedom are n-1 where n is types of dinucleotide pairs occurring at i,j.




#filepatht="../../../Documents/Syllabi/Grabowskilab/C1_alignments/"
#filenamet= "tRNA.fasta.txt"
