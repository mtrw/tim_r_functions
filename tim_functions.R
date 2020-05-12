library(plyr)
library(dplyr)
library(magrittr)
library(data.table)
library(ggplot2)
library(parallel)
library(Rcpp)
library(colorspace)
library(zoo)



#basically for assigning multiple captures in \\1 \\2 etc in data.table calls like dt[ , c("var1","var2") := .( sub_capturevec("(\\d)_(.*)(\\d+)$",colname,3)]
    #nmatch is needed
    #NEEDS WORK
# sub_capturevec <- function(pattern,string,nmatch){
#   sapply(string,function(st){
#     while( length(grep((sep <- paste0(":",sample(1:1e6,1),":")),st))>0 ){}
#     sub(pattern,paste0("\\",1:nmatch,collapse=sep),st,perl=T) %>% stringi::stri_split_fixed(sep,omit_empty=T) %>% `[[`(1)
#   }) %>% alply(1,function(x) x)
# }
# sub_capturevec(pattern="(.)_(.)_(.)",string=c("a_b_cde","f_g_hij"),nmatch=3)

x <- sample.int(3,20,replace=T)

most_common_thing <- function(x,threshold_prop=0,na.rm=T,draw_out=NA,na_wins_out=NA,threshold_notmet_out=NA){
  if(na.rm){
    x <- x[!is.na(x)]
  }
  if(length(x)==0){
    stop("Length of x is zero")
  }
  x[is.na(x)] <- "NA"
  tbl <- table(x)
  Ma <- names(tbl[order(-tbl)])[1]

  if (length(table)==1){
    as(Ma,class(x))
  }
  else if (tbl[order(-tbl)][1]==tbl[order(-tbl)][2]){
    draw_out
  } else if (Ma=="NA"){
    na_wins_out
  } else if (tbl[order(-tbl)][1]/sum(tbl) < threshold_prop){
    threshold_notmet_out
  } else {
    as(Ma,class(x))
  }
}
#most_common_thing(x=c(1,1,1,2,2,2,3,3,NA,NA,NA,NA,NA),draw_out="DRAW!",na_wins_out="na was the most common",na.rm=T)

#pca is a SNPRelate PCA object. Tables must contain at least a column called "sample.id" (optional col, cex, pch)
plot_pca_snprelate <- function(pca=pca,col_size_table=data.table(sample.id=pca$sample.id,col=1,pch=20,cex=1),pcs=c(1,2),scree=0L){
  if(scree>0L){
    plot(x=1:scree,y=pca$eigenval[1:scree]/sum(pca$eigenval),main="Scree plot",ylab="Proportion of variance explained",xlab="Principal Component",pch=20)
    wait("Please press [ENTER] to continue to the next plot.")
  }
  dt <- data.table(
    sample.id = pca$sample.id,
    x = pca$eigenvect[,pcs[1]] * pca$eigenval[pcs[1]],
    y = pca$eigenvect[,pcs[2]] * pca$eigenval[pcs[2]]
  )
  dt <- col_size_table[dt,on="sample.id"]
  if(!"cex" %in% colnames(dt)) dt[,cex:=1]
  if(!"pch" %in% colnames(dt)) dt[,pch:=20]
  if(!"col" %in% colnames(dt)) dt[,col:=1]
  plot(
    x=dt$x,
    y=dt$y,
    col=dt$col,
    cex=dt$cex,
    pch=dt$pch,
    xlab=paste0("PC",pcs[1]),
    ylab=paste0("PC",pcs[2])
  )
}

#give it "st" as AG or CT or AGT etc, must be sorted alphabetically.
IUPAC <- function(st){
  codes <- c(        "R" ,  "Y" ,  "S" ,  "W" ,  "K" ,  "M" ,  "B" ,   "D" ,   "H" ,   "V" )
  names(codes) <- c( "AG" , "CT" , "CG" , "AT" , "GT" , "AC" , "CGT" , "AGT" , "ACT" , "ACG" )
  codes[st]
}
#IUPAC(c("AG","AGT","CT"))

#minor_allele_frequency(sample(1,20,r=T)) #BY STATE FOR HAPLOIDS, ONLY ... NOT FOR 0,1,2 ENCODING OF DIPLOIDS
minor_allele_frequency <- function(x){
  tbl <- table(x)
  if(length(tbl)==1){ return(as.numeric(0)) }
  min <- tbl[order(tbl)][1]
  if ((min / sum(tbl))==1) { browser() }
  as.numeric(min / sum(tbl))
}


#minor_allele(sample(1:3,20,r=T))
minor_allele <- function(x){
  tbl <- table(x)
  ma <- names(tbl[order(tbl)])[1]
  return(as(ma,class(x)))
}

#great when reading in JPGs, which are often three arrays giving the r g b .. to plot in base
  #rgb is a three item vector
rgb2hex <- function(rgb) sprintf('#%s',paste(as.hexmode(rgb),collapse = ''))

#I hate those empty dots base plot defaults to
par(pch=20)

#just to make the Rcpp compiler use the right C standard.
Sys.setenv("PKG_CXXFLAGS"="-std=gnu++11")

#be tolerant of non-utf-8 characters when grepping
Sys.setlocale('LC_ALL','C')

#create an empty plot with ranges x=c(low,high) and y=ditto
null_plot <- function(x,y,xlab=NA,ylab=NA,...){
  plot(NULL,xlim=range(x,na.rm=T),ylim=range(y,na.rm=T),xlab=xlab,ylab=ylab,...)
}

#turn levels of a factor into colours from a colorspace palette (in the diverge_hcl set)
replace_levels_with_colours <- function(x,palette="Berlin",alpha=1,fun="diverge_hcl",plot=FALSE,newplot=TRUE){
  require(colorspace)
  n <- nu(x)
  cols <- match.fun(fun)(n,palette = palette,alpha = alpha)
  colvec <- swap( x , unique(x) , cols , na.replacement = NA )
  if(plot==FALSE) {
    return(colvec)
  } else {
    # null_plot(y=1:length(cols),x=rep(1,length(cols)),xaxt="n",yaxt="n")
    # text(y=1:length(cols),x=rep(1,length(cols)),labels=unique(x),col=cols)
    if(newplot) {null_plot(x=0,y=0,xaxt="n",yaxt="n",bty="n")}
    legend(x="topleft",legend=unique(x),fill=cols,text.col=cols)
  }
}

#turn levels of a factor into colours from a colorspace palette (in the diverge_hcl set)
replace_scale_with_colours <- function(x,palette="ag_GrnYl",fun="sequential_hcl",alpha=1,plot=FALSE){
  require(colorspace)
  s <- x %>% scale_between(1,100) %>% round
  cols <- match.fun(fun)(100,palette = palette,alpha = alpha)
  if(plot==FALSE) {
    return(cols[s])
  } else {
    plot(y = (1:100)%>%scale_between(min(x,na.rm=T),max(x,na.rm=T)),x=rep(0,100),ylab=NA,xlab=NA,xaxt="n",col=cols)
    return(cols[s])
  }
}

#scan a strong of positions and earmark regions of low density below a threshold (also plots to help you choose)
mark_low_density <- function(vec,bandwidth=NULL,min_run_length=NULL,n_bins=1e6,return_bins_dt=FALSE,plot=TRUE,min_sd=NULL,min_abs=NULL,...){
  if(is.null(bandwidth)){bandwidth <- "nrd0"}
  if((is.null(min_sd) & is.null(min_abs)) | (!is.null(min_sd) & !is.null(min_abs))){ stop("Must give one of min_abs and min_sd") } else if (!is.null(min_abs)){min<-min_abs} else if (!is.null(min_sd)){min<-min_sd}
  binsize <- diff(range(vec))/n_bins
  d <- density(vec,bw=bandwidth,n=n_bins+2,from=min(vec)-binsize,to=max(vec)+binsize,...)
  if(!is.null(min_sd)){min <- mean(d$y)-(sd(d$y)*min_sd)}
  if(plot){
    plot(d$x,d$y,type="l")
    abline(h = min)
  }
  
  filter <- d$y<min #TRUE mean it's too low and should be filtered
  if(!is.null(min_run_length)){
    r <- rle(filter)
    r$values[r$values==TRUE & r$lengths<min_run_length] <- FALSE
    filter <- inverse.rle(r)
  }
  if(return_bins_dt){
    return(data.table(
      x=d$x,
      include=filter
    ))
  }
  dtv <- data.table(vec=vec)
  dtf <- data.table(filter=filter,vec=d$x)
  dtf[dtv,on=.(vec),roll=T]
}
# ex <- data.table(
#   pos={ patches <- runif(50)*1e9; rnorm(1e6,rep(patches,each=1e6/50),rep(runif(50,1e5,1e8),each=1e6/50)) },
#   chr=rep(c(1L,2L),each=0.5e6)
# )
# mark_low_density(vec=ex[chr==1]$pos,bandwidth=NULL,min_sd=0,n_bins=45,min_run_length=12,return_bins_dt=FALSE)
# ex[ , mark_low_density(pos,min_sd=-.2) , by=.(chr)]


#n stats for genomes
Nstat <- function(x,n=50L){
  x <- as.numeric(x[order(x)])
  x[cumsum(x) > (sum(x) * (100-n)/100)][1]
}

#work it out an annotate it
pgrep <- function(search,in_me){
  sapply(search,function(x) in_me[grep(x,in_me)[1]] )
}

#display brewer palletes and remind yourself how to use them
cs <- function(){
  library(colorspace)
  hcl_palettes(plot = TRUE)
  ce("Example: diverge_hcl(30,palette=\"Berlin\")")
  ce("demoplot(x, type = c(\"map\", \"heatmap\", \"scatter\", \"spine\", \"bar\", \"pie\",
\"perspective\", \"mosaic\", \"lines\"), ...)")
}

#random integers
rint <- function(n,from,to,replace=T){
  sample.int(abs(to-from)+1,n,replace=replace) + from - 1L
}

#apply a fun of two vars in every combo, and give the results as a matrix
#grid_ply(1:2,3:1,sum)
grid_ply <- function(rows,cols,FUN) {
  lapply( rows , function(row) lapply(cols , function(col) {
    FUN(row,col)
  }) %>% unlist ) %>% unlist %>%
    matrix(nrow=length(rows),dimnames=list(rows,cols))
}

#apply a fun of two vars in every combo, and give the results as a matrix
#mc_grid_ply(1:2,3:1,sum)
mc_grid_ply <- function(rows,cols,FUN,cores=25,...) {
  mclapply( mc.cores=cores , rows , function(row) lapply(cols , function(col) {
    FUN(row,col)
  }) %>% unlist ) %>% unlist %>%
    matrix(nrow=length(rows),dimnames=list(rows,cols))
}

#scale a list of values to between two points, proportionally spaced as they were originally
#rnorm(100) %>% scale_between(20,29) %>% pd
scale_between <- function(x,lower,upper){
  if(all(x==mean(x,na.rm=T))) return(rep(mean(c(lower,upper),na.rm=T),length(x)))
  ( x - min(x,na.rm=T) ) / (max(x,na.rm=T)-min(x,na.rm=T)) * (upper-lower) + lower
}

#easy way to see the spread of your values. plot the density.
#pd(c(NA,rnorm(500),NA))
pd <- function(x,add=F,...){
  if(!add){
    x %>% density(na.rm=TRUE,...) %>% plot()
  } else {
    x %>% density(na.rm=TRUE,...) %>% lines()
  }
}

#difference between two values but vectorises more intuitively than diff()
#difff(10,-286)
difff <- function(a,b){
  abs(a-b)
}

#z-transform values in a vec
#z_transform(x = ((runif(100))**2)+20  ) %>% pd
z_transform <- function(x){
  if ( sd(x,na.rm=TRUE)==0 ){
    a <- rep(0,times=length(x))
    a[is.na(x)] <- NA
    return(a)
  }
  b <- (mean(x,na.rm=T)-x) / sd(x,na.rm=T)
  b[is.na(x)] <- NA
  b
}

#number of unique entries in vector
#nu(c(1,1,2,3,4,4))
nu <-function(x){
  unique(x) %>% length
}

#fetch the upper triangular matrix in a vector without all that typing
#upper_tri( matrix(1:9,ncol=3) , diag=T )
upper_tri <- function(mat,diag=F){
  mat[upper.tri(mat,diag=diag)]
}

#parallel mean of two numbers. Great for BLAST search (midpoint of hits = pmean(start_vec,end_vec))
#pmean(1:5,2:6,5:1)
pmean <- function(...){
  invecs <- list(...)
  out <- rep(0, times=length(invecs[[1]]) )
  lapply(invecs,function(x) out <<- out+x )
  out/length(invecs)
}
pmean2 <- pmean #legacy reasons


#read a gff3 format file
read_gff <- function(fname){
  fread( cmd=paste("cat",fname,"| grep -v '^#'") ,col.names=c("seqname","source","feature","start","end","score","strand","frame","attribute"))
}

#as above but
#not the parts (eg CDSs ... just the whole genes)
read_gff_genes_pgsb <- function(fname){
  g <- fread( paste("cat",fname,"| grep -v '^#' | awk '$3==", '"gene"' ,"'") ,col.names=c("seqname","source","feature","start","end","score","strand","frame","attribute"))
  g[ , id := sub("^ID=(.*);.*","\\1",attribute) ][]
}
#seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
# source - name of the program that generated this feature, or the data source (database or project name)
# feature - feature type name, e.g. Gene, Variation, Similarity
# start - Start position of the feature, with sequence numbering starting at 1.
# end - End position of the feature, with sequence numbering starting at 1.
# score - A floating point value.
# strand - defined as + (forward) or - (reverse).
# frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
# attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.

#for reading in my fave blast output
#as produiced by ~/bin/blastn_basic.zsh
read_blastn_basic <- function(fname){
  fread(fname,col.names=c("qseqid", "sseqid" , "slength" , "qlength" , "match_len" , "qstart" , "qend" , "sstart" , "send" , "pct_id" , "evalue" , "bitscore" ))
}

#length, because why spend your life typing "ength" all the time?
l <- function(x){
  stop("Naughty length shortcut used! Please replace with call to `length`() ...")
}

#simultaneously change the names of things to a particular thing if they match a particular string.
#name_by_match(vec=iris$Species,matches = c("set","sicol"),names = c("SETOSA","VERSICOLOR"))
name_by_match <- function(vec,matches,names){
  if(length(matches) != length(names)) { stop("Lengths of `matches` and `names` vectors don't match, you massive knob!")  }
  vec %<>% as.character()
  l_ply(1:length(matches) , function(n){
    vec[grep(matches[n],vec)] <<- names[n]
  })
  vec
}

#simultaneously change the names of things to a particular thing if they match (EXACTLY!) a particular string.
#swap(vec=iris$Species,matches = c("virginica"),names = c("XXXXX"))

swap <- function(vec,matches,names,na.replacement=NA){
  orig_vec <- vec
  #if(sum(! matches %in% names ) > 0 ) { stop("Couldn't find all matches in names") }
  if(length(matches) != length(names)) { stop("Lengths of `matches` and `names` vectors don't match, you old bison!") }
  if(is.factor(vec)) { levels(vec) <- c(levels(vec),names,na.replacement) }
  vec[is.na(orig_vec)] <- na.replacement
  l_ply( 1:length(matches) , function(n){
    vec[orig_vec==matches[n]] <<- names[n]
  })
  vec
}

#length of a vector
#(c(3,4))
vec_length <-function(x){
  ( sum( x**2 ) )**(1/2)
}

#euc dists between points described in lists of coords
euc_dist<-function(x1,x2,y1,y2){
  sqrt( ((x1-x2)**2) + ((y1-y2)**2) )
}


#as the title says. relative probs should be positive (duh) and same order as events (duh)
#t <- random_draws_from_any_discreet_distribution(n=50000,events=LETTERS[1:5],relative_probs=1:5) %>% table; t / max(t) * 5 ; rm(t)
random_draws_from_any_discreet_distribution <- function(n=1,events,relative_probs){
  lapply( 1:n,   function(x){
    events[ which( runif(1)*sum(relative_probs) < (cumsum(relative_probs)) )[1] ]    
  }) %>% unlist
}

#n random (non-repeated!) rows of a data frame
#sample_df(iris,20)
sample_df <- function(df,n=10,...){
  df[ sample.int(nrow(df),n,...) ]
}


#random_matrix( 5 , 5 , symmetric = T, all_pos=F )
random_matrix <- function(m,n=m,symmetric=F,all_pos=T){
  if(all_pos){
    M <- matrix(runif(m*n),nrow=m)
  } else {
    M <- matrix(runif(m*n)-0.5,nrow=m)
  }
  if(symmetric){
    M <- M+t(M)
  }
  M
}

#the middle of the range of a dataset
#midpoint(c(0,4,4,4,4,4,4,4,4,4,10,NA,10))
midpoint <- function(x){
  mean( c( min(x, na.rm=T) , max(x, na.rm=T) ) )
}

#the expected value this many SDs from the mean of distribution (assumed normal, taken from a vector)
#SDs_from_mean(rnorm(10000),2) #around 2
SDs_from_mean <- function(x,SDs){
  mean(x) + SDs*sd(x,na.rm=T)
}

#echo to the global environment. good warning messager. still doesn't work in mclappy, hashtag doh
#ce("beans ",list("hello "," goodbye")," whatever")
ce <- function(...){   cat(paste0(...,"\n"), sep='', file=stderr()) %>% eval(envir = globalenv() ) %>% invisible() }

#distance to the nearest element of each element
nearest_in_set <- function(x) { 
  if(length(x)==1) {return(0)} else  { sapply(1:length(x), function(y) { min (  abs( x[y] - x[setdiff(1:length(x),y)] ) ) } ) } 
}
#nearest_in_set(1:5); nearest_in_set((1:5)**2)



#sigmoid curves
sigmoid <- function(x,L=1,k=1,xo=0){
  L/(1+exp(-k*(x-xo)))
}
#sigmoid(-10:10)

interpolate <- function( x , y , points_x=NULL , points_y=NULL  ){
  if((is.null(points_x) & is.null(points_y)) | ((!is.null(points_y) & !is.null(points_x)) )){
    stop("Must enter one of points_x and points_y, not both.")
  }
  if(!is.null(points_y)){
    x -> t
    x <- y
    y <- t
    rm(t)
    points_x <- points_y
    rm(points_y)
  }
  
  len <- length(x)
  if(len!=length(y) | len<2){
    stop("x and y must be equal lengths of at least 2.")
  }
  
  o_x <- order(x)
  o_y <- order(y)
  if((length(o_x != o_y))<0){
    stop("x and y must be monotonic.")
  }
  x <- x[o_x]
  y <- y[o_x]
  
  out <- lapply(points_x,function(pt_x){
    #it hits an exact point
    if(length(xeq <- which(x==pt_x))==1){
      y[xeq]
    } else if(length(xeq)>1) {
      warning("WARNING: Falls within run of identical x-values. Estimating by mean(y).")
      mean(y[xeq])
    } else if(pt_x>max(x)){
        warning("WARNING: Extrapolating at high end.")
      if((x[len]-x[len-1])==0){
        warning("WARNING: Falls beyond run of identical x-values. Estimating by mean(y).")
        mean(y[x==x[len]])
      } else {
        y[len] + (pt_x-x[len])*( (y[len]-y[len-1])/(x[len]-x[len-1]) )
      }
    } else if(pt_x<min(x)){
      warning("WARNING: Extrapolating based on last two points at low end.")
      if((x[2]-x[1])==0){
        warning("WARNING: Falls before run of identical x-values. Estimating by mean(y).")
        mean(y[x==x[1]])
      } else {
          y[1] + (pt_x-x[1])*( (y[2]-y[1])/(x[2]-x[1]) )
        }
    } else {
      nearest_under_idx <- last(which(x<pt_x))
      y[nearest_under_idx] + (pt_x-x[nearest_under_idx])*( (y[nearest_under_idx+1]-y[nearest_under_idx])/(x[nearest_under_idx+1]-x[nearest_under_idx]) )
    }
  }) %>% unlist
  out
}

#just makes things easier
view_matrix <- function(mat,n=10){
  m <- cbind(mat[1:n,1:n],rep("...",n))
  m <- rbind(m,rep("...",n+1))
  m
}




number_runs <- function(vec){
  x <- copy(vec)
  r <- runif(1)
  while( r%in%x==TRUE ) { r <- runif(1) }
  x[ is.na(x) ] <- r
  o <- 1
  n <- 1
  l <- x[1]
  for(i in x[2:length(x)]){
    if( i!=l ){
      n <- n + 1
    }
    o <- append(o,n,length(o))
    l <- i
  }
  o
}
#number_runs( c(5,5,5,NA,1,1,"a","a","b",1,1,1,1,FALSE,TRUE,TRUE) )

#top left
tl <- function(dt,n=10){
  dt[1:n,1:n]
}
iris %>% tl(3)

#top right
tr <- function(dt,n=10){
  dt[1:n,(ncol(dt)-n):ncol(dt)]
}
iris %>% tr(3)

u <- function(...){
  unique(...)
}

`%prop_in%` <- function(a,b){
  ce("Prop(A in B): " , sum( a %in% b ) / length(a))
  ce("Prop(B in A): " , sum( b %in% a ) / length(b))
}
#1:10 %prop_in% 5:15

#a simple ggplot theme
reduced_l <- theme(
  axis.line = element_line(colour = "black"),
  
  panel.grid.major = element_line(size=.1 , colour = "grey"),
  panel.grid.minor = element_blank(),
  panel.background = element_rect( fill = "transparent", colour = "black"),
  strip.background = element_rect(fill = "transparent", colour = "black"),
)

wait <- function(message="Press [enter] to continue"){
  invisible(readline(prompt=message))
}

#violin plots in base. could use tweaking to make various things controllable.
violin_plot <- function(x,y){
  data <- data.table(
    data_x=x,
    data_y=y
  )
  #calculte densities for violin plot
  violin_data <- data[ , .(density=.(density(data_y))) , by=.(data_x) ]
  
  #calculate how wide to plot violins
  violin_width <- data[ order(data_x) , {s <- data_x %>% unique %>% diff %>% summary; s[[4]]*0.4} , ]
  
  #set up plot
  plot(NULL,xlim=range(data$data_x),ylim=range(data$data_y),xlab=NA,ylab=NA)
  
  #plot mean as a line
  data[ , .(mean_y = mean(data_y)) , by=.(data_x) ][, lines(data_x,mean_y) ]
  
  
  #plot points and violins
  lapply( unique(data$data_x) , function(x) {
    #points
    with( data[data_x==x] , points(data_x,data_y,pch=20,cex=.2) )
    #violin
    d <- violin_data[data_x==x]$density[[1]]
    scale_by <- violin_width/max(d$y)
    lines(x+d$y*scale_by,d$x,cex=.2)
    lines(x-d$y*scale_by,d$x,cex=.2)
  })
}

#make some fake data and look at it
# n_y_each_x <- sample(20:50,10,r=T)
# p <- data.table(
#   data_x=(rep(1:10,n_y_each_x) -> t),
#   data_y=rnorm(length(t),t**2,t)
# )
# rm(t)
# 
# violin_plot(p$data_x,p$data_y)


#inverses of the effect of doing letters[indices] and LETTERS[indices]
numbers <- function(x){
  sapply(x,function(l) which(letters==l))
}
NUMBERS <- function(x){
  sapply(x,function(l) which(LETTERS==l))
}

#an n-wheeled enigma machine with arbitrary alphabet
enigma_machine <- function(message,setting=NULL,alphabet="standard"){
  require(plyr)
  require(gtools)
  if(length(alphabet)==1){
    if(alphabet=="standard"){alphabet<-chr(32:126)} else if (alphabet=="ascii"){alphabet<-chr(1:255)} else if (alphabet=="alphanumeric"){alphabet<-c(LETTERS,0:9," ",".") ; message <- toupper(message) }
  }
  nsymbols <- length(alphabet)
  if(is.null(setting)) {setting<-alphabet[1:min(3,nsymbols)]}
  nwheel <- length(setting)
  message_vec <- strsplit(message,"")[[1]]
  if(any(! message_vec %in% alphabet)){ stop("Illegal characters in message. Change alphabet (to \"ascii\"?) or message.") }
  
  #store the wheel settings in enigma. Each wheel is represented twice (once for forwards and once for backwards traversal)
  wheels <- matrix(rep(1:nsymbols,nwheel),ncol=nwheel) #define enigma here, and set "reverse" mappings by rotating them backwards
  cap <- nsymbols:1 #any mapping that always swaps entries! (i.e pos3 == 4 => pos4 == 3). This solution is just an easy option (besides the trivial 1:nsymbols, which makes a symmetrical mapping on either side of the cap and doesn't actually encode your message)
  enigma <- cbind(wheels,cap,wheels)
  
  #instructions are given per wheel, but it changes all columns as it should
  .rotate <- function(col,progress){
    col <- c( col , (2*nwheel+2)-col )
    progress <- progress %% nsymbols
    progress <- c( progress  , nsymbols - progress )
    progress <- progress %% nsymbols
    l_ply(1:length(col),function(i){
      c <- col[i]
      p <- progress[i]
      enigma[,c] <<- enigma[ c((p+1):nsymbols,0:(p)) , c ]
    })
  }
  
  #rotate to initial settings
  .rotate(1:nwheel,sapply(setting,function(s) which(alphabet==s))-1)

  .keystrike <- function(key){
    k <- which(alphabet==key)
    for(i in 1:ncol(enigma)){
      k <- enigma[k,i]
    }
    alphabet[k]
  }
  
  #to turn the wheels as the [en|de]coding takes place
  intervals <- nsymbols**(0:(nwheel-1))
  
  encoded <- character()
  
  #run the machine
  for( i in 1:length(message_vec) ){
    char <- message_vec[i]
    .keystrike(char) -> k
    .rotate(col=which((i %% intervals)==0)->t,progress=rep(1,length(t)))
    encoded[i] <- k
  }
  paste0(encoded,collapse = "")
}


enigma_machine(
                message="hello sailor",
                setting=c("e","s","h"),
                alphabet=c(letters," ")
              ) -> c

#c
enigma_machine(
  message=c,
  setting=c("e","s","h"),
  alphabet=c(letters," ")
)






most_frequent_by_margin <- function(x,f1_vs_f2=2){
  if(nu(x)==1){return(x[1])}
  top2 <- (table(x) %>% sort(decreasing = T))[1:2]
  if(top2[1]/top2[2] >= f1_vs_f2){
    return(as(names(top2)[1],class(x)))
  } else {
    as(NA,class(x))
  }
}
#most_frequent_by_margin(c(0,0,0,1,1,1,1,2,2,3,4,5,6))
#most_frequent_by_margin(c(0,0,1,1,1,1,2,2,3,4,5,6))
