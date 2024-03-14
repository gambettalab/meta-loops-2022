#!/usr/bin/env Rscript


## Copyright (C) 2022 Julien Dorier and UNIL (University of Lausanne).
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or (at
## your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
## General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <https://www.gnu.org/licenses/>.



#############################################
#parse command line options
#############################################
          
args=commandArgs(trailingOnly = TRUE)
library(optparse)
option_list=list( 
    make_option(c("--chunk-size"), type="integer", default=5000,metavar="size",
                help="chunk size (unit: number of bins) [default %default]."),
    make_option(c("--resolution"), type="integer", default=2000,metavar="RES",
                help="use hic data at resolution RES in bp [default %default]."),
    make_option(c("--chrs"), type="character", default="chr2L,chr2R,chr3L,chr3R,chr4,chrX",metavar="CHRS",
                help="comma separated list of chromosomes to consider
\t\t[default %default]."),
    make_option(c("--score-thresh"), type="double", default=30,metavar="score",
                help="select bin pairs with obs/exp>SCORE as candidate
\t\tregion [default %default]."),
    make_option(c("--count-thresh-dist"), type="double", default=100000,metavar="DIST",
                help="filter out bin pairs with normalized count<threshold,
\t\twith thesholds taken as the mean normalized count at
\t\tdistance DIST (in bp, per chromosome) [default %default]."),
    make_option(c("--clustering-distance"), type="double", default=5,metavar="DIST",
                help="group selected bin pairs within DIST (unit: number
\t\tof bins) [default %default]."),
    make_option(c("--min-count"), type="integer", default=10,metavar="N",
                help="filter out clusters with less than N selected bin
\t\tpairs [default %default]."),
    make_option(c("--min-size"), type="integer", default=5,metavar="N",
                help="filter out clusters with width or height below N
\t\t(unit: number of bins) [default %default]."),
    make_option(c("--blur-sigma"), type="double", default=1,metavar="TOL",
                help="standard deviation (unit: number of bins) of the
\t\tgaussian blur filter used for meta-loop detection [default %default]."),
    make_option(c("--watershed-tolerance"), type="double", default=5,metavar="TOL",
                help="watershed tolerance used for meta-loop detection [default %default].
\t\tCorresponds to the minimum height (obs/exp) of the
\t\tmeta-loop between highest point and point where it
\t\tcontacts another object. If the height is smaller than
\t\tthe tolerance, the object will be merged with its
\t\t highest neighbor."),
    make_option(c("--loop-filter-param"), type="double", default=2,metavar="P1",
                help="filter out meta-loops with obs/exp < max(score)/P1 [default %default]."),
    make_option(c("--output"), type="character", default="output.tsv",metavar="FILENAME",
                help="output filename (tsv file format).")
    )


opt=parse_args(OptionParser(option_list=option_list,
                            usage= "usage: %prog [options] MCOOL",
                            description="Positional argument:
\tMCOOL
\t\tinput mcool file generated with cooler."),positional_arguments=1,args=args)



chunk.size=opt$options[["chunk-size"]]
binsize=opt$options[["resolution"]]
threshold_score=opt$options[["score-thresh"]]
threshold_count_dist=opt$options[["count-thresh-dist"]]
threshold_distance_regions=opt$options[["clustering-distance"]]
threshold_min_cluster_count=opt$options[["min-count"]]
threshold_min_cluster_size=opt$options[["min-size"]]
loop.filter.param=opt$options[["loop-filter-param"]]
blur.sigma=opt$options[["blur-sigma"]]
watershed.tolerance=opt$options[["watershed-tolerance"]]
chrs=strsplit(opt$options[["chrs"]],",")[[1]]
outputfile=opt$options[["output"]]
##positional arguments
inputfile=opt$args[1]


library(data.table)
library(EBImage)
library(rhdf5)
library(mlpack)
library(igraph)

options(scipen=10000)
dir.create(dirname(outputfile), showWarnings = FALSE,recursive=TRUE)


#############################################
#simple functions to get bins and hic counts from mcool
#############################################
##return bins for chr
get_hic_bins=function(f,binsize,chr){
    ##note:
    ## schema: see https://cooler.readthedocs.io/en/latest/schema.html
    ## index in h5read is 1-based,  bin_id is 0-based
    
    ##load all bins for chr
    chr.idx=which(h5read(f,paste0("/resolutions/",binsize,"/chroms/name"))==chr)
    chrom_offset=h5read(f,paste0("/resolutions/",binsize,"/indexes/chrom_offset"))[c(chr.idx,chr.idx+1)]
    bins=data.table(bin.index=(chrom_offset[1]):(chrom_offset[2]-1),
                    chr=as.character(h5read(f,paste0("/resolutions/",binsize,"/bins/chrom"),index=list((chrom_offset[1]+1):(chrom_offset[2])))),
                    bin.start=as.integer(h5read(f,paste0("/resolutions/",binsize,"/bins/start"),index=list((chrom_offset[1]+1):(chrom_offset[2])))),
                    bin.end=as.integer(h5read(f,paste0("/resolutions/",binsize,"/bins/end"),index=list((chrom_offset[1]+1):(chrom_offset[2])))),
                    weight=as.numeric(h5read(f,paste0("/resolutions/",binsize,"/bins/weight"),index=list((chrom_offset[1]+1):(chrom_offset[2])))))
    return(bins)
}
##return hic sparse matrix (upper triangle) for chr1 x chr2
get_hic_sparse=function(f,binsize,chr1,chr2){
    ##note:
    ## schema: see https://cooler.readthedocs.io/en/latest/schema.html
    ## index in h5read is 1-based,  bin_id is 0-based

    ##load all bins for chr1
    chr.idx=which(h5read(f,paste0("/resolutions/",binsize,"/chroms/name"))==chr1)
     ## ##chrom_offset: indicates which row in the bin table each chromosome first appears. The last element stores the length of the bin table.
    chrom_offset=h5read(f,paste0("/resolutions/",binsize,"/indexes/chrom_offset"))[c(chr.idx,chr.idx+1)]
    bins1=data.table(chrom=as.character(h5read(f,paste0("/resolutions/",binsize,"/bins/chrom"),index=list((chrom_offset[1]+1):(chrom_offset[2])))),
                     start=as.integer(h5read(f,paste0("/resolutions/",binsize,"/bins/start"),index=list((chrom_offset[1]+1):(chrom_offset[2])))),
                     end=as.integer(h5read(f,paste0("/resolutions/",binsize,"/bins/end"),index=list((chrom_offset[1]+1):(chrom_offset[2])))),
                     weight=as.numeric(h5read(f,paste0("/resolutions/",binsize,"/bins/weight"),index=list((chrom_offset[1]+1):(chrom_offset[2])))),
                     bin_id=(chrom_offset[1]):(chrom_offset[2]-1))

    ##load all bins for chr2
    chr.idx=which(h5read(f,paste0("/resolutions/",binsize,"/chroms/name"))==chr2)
    chrom_offset=h5read(f,paste0("/resolutions/",binsize,"/indexes/chrom_offset"))[c(chr.idx,chr.idx+1)]
    bins2=data.table(chrom=as.character(h5read(f,paste0("/resolutions/",binsize,"/bins/chrom"),index=list((chrom_offset[1]+1):(chrom_offset[2])))),
                     start=as.integer(h5read(f,paste0("/resolutions/",binsize,"/bins/start"),index=list((chrom_offset[1]+1):(chrom_offset[2])))),
                     end=as.integer(h5read(f,paste0("/resolutions/",binsize,"/bins/end"),index=list((chrom_offset[1]+1):(chrom_offset[2])))),
                     weight=as.numeric(h5read(f,paste0("/resolutions/",binsize,"/bins/weight"),index=list((chrom_offset[1]+1):(chrom_offset[2])))),
                     bin_id=(chrom_offset[1]):(chrom_offset[2]-1))

    ##create upper triangular matrix (sparse format)
    bin1_offset=as.integer(h5read(f,paste0("/resolutions/",binsize,"/indexes/bin1_offset")))
    index_first=bin1_offset[bins1[,min(bin_id)]+1]+1
    index_last=bin1_offset[bins1[,max(bin_id)]+2]
    hic.data=data.table(bin1_id=as.integer(h5read(f,paste0("/resolutions/",binsize,"/pixels/bin1_id"),start=index_first,block=index_last-index_first+1,count=1)),
                          bin2_id=as.integer(h5read(f,paste0("/resolutions/",binsize,"/pixels/bin2_id"),start=index_first,block=index_last-index_first+1,count=1)),
                          count=as.numeric(h5read(f,paste0("/resolutions/",binsize,"/pixels/count"),start=index_first,block=index_last-index_first+1,count=1)))
    hic.data=hic.data[bin2_id%in%bins2$bin_id]

    ##add chr, bin.start, bin.end, weight.
    setkey(bins1,bin_id)
    setkey(hic.data,bin1_id)
    hic.data[,weight1:=bins1[hic.data,weight]]
    setkey(bins2,bin_id)
    setkey(hic.data,bin2_id)
    hic.data[,weight2:=bins2[hic.data,weight]]
    ##set invalid bin to NaN
    hic.data[!(is.finite(weight1)&is.finite(weight2)),count:=NaN]
    ##normalize
    hic.data[,count.normalized:=count*weight1*weight2]
    return(hic.data[,.(bin.index1=bin1_id,bin.index2=bin2_id,count.normalized)])
    
}

#############################################
## Check
#############################################

###check resolution
tmp=h5ls(inputfile)
resolutions=sort(as.numeric(tmp[tmp$group=="/resolutions","name"]))
if(!binsize%in%resolutions)
    stop(inputfile,": resolution ",binsize," not found. Valid resolutions: ",paste(resolutions,collapse=", "))


for(res in resolutions)
{
    ##check upper triangular matrix
    metadata=attributes(h5read(inputfile,paste0("/resolutions/",res),read.attributes=TRUE,count=1))
    if(metadata[["storage-mode"]]!="symmetric-upper")
        stop(inputfile,": storage-mode must be symmetric-upper.")

    ##check chrs
    tmp=as.character(h5read(inputfile,paste0("/resolutions/",res,"/chroms/name")))
    if(!all(chrs%in%tmp))
    {
        stop(inputfile,": chromosome ",chrs[!chrs%in%tmp]," not found. Valid chromosome names: ",paste(tmp,collapse=", ")) 
    }
}


#############################################
# Mean count vs distance
#############################################

cat("Evaluating mean count vs distance\n")
count.vs.distance=rbindlist(lapply(chrs,function(ch){
    cat(" ",ch,"\n")
    hic.bins.chr=get_hic_bins(inputfile,binsize,ch)
    setorder(hic.bins.chr,bin.index)
    hic.data.chr=get_hic_sparse(inputfile,binsize,ch,ch)
    setkey(hic.data.chr,bin.index1,bin.index2)
    
    count.vs.distance.tmp=NULL
    
    tmp=rbindlist(lapply(split(hic.bins.chr,hic.bins.chr[,.I%/%chunk.size],sorted=FALSE),function(hic.bins.chunk1){
        rbindlist(lapply(split(hic.bins.chr,hic.bins.chr[,.I%/%chunk.size],sorted=FALSE),function(hic.bins.chunk2){
            ##add missing valid bins (i.e. with finite weight)
            hic.data.chunk=hic.data.chr[CJ(bin.index1=hic.bins.chunk1[is.finite(weight),bin.index],bin.index2=hic.bins.chunk2[is.finite(weight),bin.index])[bin.index1<=bin.index2],on=.(bin.index1,bin.index2)]
            ##replace NA by 0
            hic.data.chunk[is.na(count.normalized),count.normalized:=0]
            
            ##add bin.start
            setkey(hic.bins.chunk1,bin.index)
            setkey(hic.data.chunk,bin.index1)
            hic.data.chunk[,bin.start1:=hic.bins.chunk1[hic.data.chunk,bin.start]]
            setkey(hic.bins.chunk2,bin.index)
            setkey(hic.data.chunk,bin.index2)
            hic.data.chunk[,bin.start2:=hic.bins.chunk2[hic.data.chunk,bin.start]]
            
            hic.data.chunk[,.(count.normalized.sum=sum(count.normalized),n=.N),by=.(distance=bin.start2-bin.start1)]
        }))[,.(count.normalized.sum=sum(count.normalized.sum),n=sum(n)),by=.(distance)]
    }))[,.(count.normalized.sum=sum(count.normalized.sum),n=sum(n)),by=.(distance)]
    count.vs.distance.tmp=tmp[,.(distance,count.normalized.mean=count.normalized.sum/n)]
    ##spline smoothing
    spline.mean=smooth.spline(log(count.vs.distance.tmp[distance>0&count.normalized.mean>0&is.finite(count.normalized.mean),distance]),log(count.vs.distance.tmp[distance>0&count.normalized.mean>0&is.finite(count.normalized.mean),count.normalized.mean]),spar=0.8)
    count.vs.distance.tmp[,count.normalized.mean.smooth:=exp(predict(spline.mean,log(distance))$y)]
    
    ##replace inf and <0 by nan
    count.vs.distance.tmp[!is.finite(count.normalized.mean.smooth),count.normalized.mean.smooth:=NaN]
    count.vs.distance.tmp[count.normalized.mean.smooth<0,count.normalized.mean.smooth:=NaN]
    
    count.vs.distance.tmp[,.(chr=ch,distance,count.normalized.mean.smooth)]
}))
    
#############################################
## coverage
#############################################
cat("Evaluating coverage\n")

binsize.coverage=128000
binsize.coverage=resolutions[which.min(abs(resolutions-binsize.coverage))]
hic.bins=data.table(chr=as.character(h5read(inputfile,paste0("/resolutions/",binsize.coverage,"/bins/chrom"))),
                    start=as.integer(h5read(inputfile,paste0("/resolutions/",binsize.coverage,"/bins/start"))),
                    end=as.integer(h5read(inputfile,paste0("/resolutions/",binsize.coverage,"/bins/end"))))
hic.bins[,bin_id:=(1:.N)-1]

hic.data=data.table(bin1_id=as.integer(h5read(inputfile,paste0("/resolutions/",binsize.coverage,"/pixels/bin1_id"))),
                          bin2_id=as.integer(h5read(inputfile,paste0("/resolutions/",binsize.coverage,"/pixels/bin2_id"))),
                          count=as.numeric(h5read(inputfile,paste0("/resolutions/",binsize.coverage,"/pixels/count"))))

##count per bin (we have upper triangular matrix)
tmp1=hic.data[,.(count=sum(count)),by=.(bin_id=bin1_id)]
tmp2=hic.data[,.(count=sum(count)),by=.(bin_id=bin2_id)]
tmp=rbind(tmp1,tmp2)[,.(count=sum(count)),by=.(bin_id)]
setkey(hic.bins,bin_id)
setkey(tmp,bin_id)
coverage=tmp[hic.bins]
coverage[is.na(count),count:=0]

## flag bins with high coverage
coverage[,keep:=(count>=median(count)-2*mad(count))]

#############################################
# Step 1: find regions with high obs/exp
#############################################

cat("Step 1: find regions with high obs/exp\n")
regions=rbindlist(lapply(chrs,function(ch){        
    cat(" ",ch,"\n",sep="")
    hic.bins.chr=get_hic_bins(inputfile,binsize,ch)
    setorder(hic.bins.chr,bin.index)
    hic.data.chr=get_hic_sparse(inputfile,binsize,ch,ch)
    setkey(hic.data.chr,bin.index1,bin.index2)

    ##exclude low coverage regions
    hic.bins.chr[,excluded:=FALSE]
    for(i in coverage[keep==FALSE&chr==ch,.I])
    {
        hic.bins.chr[bin.start>=coverage[keep==FALSE&chr==ch][i,start]&bin.end<=coverage[keep==FALSE&chr==ch][i,end],excluded:=TRUE]
    }    

    threshold_count=count.vs.distance[chr==ch&distance==round(threshold_count_dist/binsize)*binsize,count.normalized.mean.smooth]
    
    tmp.selected=rbindlist(lapply(split(hic.bins.chr,hic.bins.chr[,.I%/%chunk.size],sorted=FALSE),function(hic.bins.chunk1){
        rbindlist(lapply(split(hic.bins.chr,hic.bins.chr[,.I%/%chunk.size],sorted=FALSE),function(hic.bins.chunk2){
            ##add missing valid bins (i.e. with finite weight)
            hic.data.chunk=hic.data.chr[CJ(bin.index1=hic.bins.chunk1[is.finite(weight),bin.index],bin.index2=hic.bins.chunk2[is.finite(weight),bin.index])[bin.index1<=bin.index2],on=.(bin.index1,bin.index2)]
            ##replace NA by 0
            hic.data.chunk[is.na(count.normalized),count.normalized:=0]
            
            ##add chr, bin.start, bin.end
            setkey(hic.bins.chunk1,bin.index)
            setkey(hic.data.chunk,bin.index1)
            hic.data.chunk[,c("chr1","bin.start1","bin.end1","excluded1"):=hic.bins.chunk1[hic.data.chunk,.(chr,bin.start,bin.end,excluded)]]
            setkey(hic.bins.chunk2,bin.index)
            setkey(hic.data.chunk,bin.index2)
            hic.data.chunk[,c("chr2","bin.start2","bin.end2","excluded2"):=hic.bins.chunk2[hic.data.chunk,.(chr,bin.start,bin.end,excluded)]]

            ##add smoothed mean count at corresponding distance and eval score=obs/exp
            hic.data.chunk[,count.normalized.mean.smooth:=count.vs.distance[chr==ch][hic.data.chunk[,.(distance=bin.start2-bin.start1)],count.normalized.mean.smooth,on=.(distance)]]
            hic.data.chunk[,score:=(count.normalized/count.normalized.mean.smooth)]
            hic.data.chunk[(excluded1==FALSE&excluded2==FALSE)&count.normalized>threshold_count&score>threshold_score,.(bin.index1,chr1,bin.start1,bin.end1,bin.index2,chr2,bin.start2,bin.end2,score)]
        }))
    }))
    ##cluster selected pairs of bins
    if(nrow(tmp.selected)>1)
    {
        ##minimum spanning tree
        mst=emst(as.matrix(tmp.selected[,c("bin.index1","bin.index2")]),verbose=TRUE)$output
        ##keep only edges with lengths < threshold_distance_regions
        mst=mst[mst[,3]<=threshold_distance_regions,,drop=FALSE]
        ##search for connected components
        g=graph_from_data_frame(mst[,1:2,drop=FALSE]+1, directed = FALSE,v=tmp.selected[,.(name=as.numeric(.I))])
        comp=components(g)
        ##assign connected component id to cluster id 
        tmp.selected[,cluster:=comp$membership] #assuming identical(V(g)$name,names(comp$membership))==TRUE and identical(as.numeric(V(g)$name),as.numeric(tmp.selected[,.I]))

        ## ##Note: minimum spanning tree + connected compnent is equivalent (but faster) to hclust followed by cutree:
        ## ct=cutree(hclust.vector(tmp.selected[,c("bin.index1","bin.index2")],method="single",metric="euclidean"),h=threshold_distance_regions)
        ## tmp.selected[,cluster:=ct]
    }
    else
    {
        tmp.selected[,cluster:=1]
    }

    ##clusters bounding boxes
    tmp.regions=tmp.selected[,.(chr1=unique(chr1),start1=min(bin.start1),end1=max(bin.end1),chr2=unique(chr2),start2=min(bin.start2),end2=max(bin.end2),N=.N),by=cluster]    
    cat("  ",tmp.regions[,length(unique(cluster))]," high obs/exp bins clusters\n",sep="")
    ##filter based on cluster size
    tmp.regions=tmp.regions[N>=threshold_min_cluster_count&pmin((end1-start1),(end2-start2))>=threshold_min_cluster_size*binsize]    
    cat("  ",tmp.regions[,length(unique(cluster))]," high obs/exp bins clusters (after filtering on size)\n",sep="")
    tmp.regions    
}))

cat("Found ",nrow(regions)," regions\n",sep="")

#############################################
# Step 2: find meta-loops 
#############################################
cat("Step 2: find meta-loops\n")
metaloops=NULL
##add id to our regions
regions[order(chr1,start1,end1,chr2,start2,end2),id:=1:.N]
setorder(regions,id)
metaloops=NULL
ids.todo=regions$id
expand=20 #pixels (note: 1 pixel=1 bin of size binsize) around region
hic.bins.chr=NULL
hic.data.chr=NULL
while(length(ids.todo)>0){
    id=head(regions[id%in%ids.todo,id],1)
    i=match(id,regions$id)
    
    ch=regions[i,chr1]
    w=regions[i,max(end1-start1,end2-start2)+2*binsize]
    cx=regions[i,(end1+start1)/2]
    cy=regions[i,(end2+start2)/2]
    rx=c(cx-w/2,cx+w/2)
    ry=c(cy-w/2,cy+w/2)

    cat(" ",ch,":",rx[1],"-",rx[2]," x ",ch,":",ry[1],"-",ry[2],"\n",sep="")

    ## find all data fully included in the region
    regions.tmp=regions[chr1==ch&start1>=rx[1]&end1<=rx[2]&chr2==ch&start2>=ry[1]&end2<=ry[2]]
    ids.todo=ids.todo[!ids.todo%in%regions.tmp[,id]]

    ##load matrix
    if(is.null(hic.bins.chr)||!ch%in%hic.bins.chr[,chr])
    {
        hic.bins.chr=get_hic_bins(inputfile,binsize,ch)
        setorder(hic.bins.chr,bin.index)
        hic.data.chr=get_hic_sparse(inputfile,binsize,ch,ch)
        setkey(hic.data.chr,bin.index1,bin.index2)
    }
        
    ##extract matrix (only valid bins, i.e. with finite weight)
    idx=hic.bins.chr[chr==ch&is.finite(weight)][,which(bin.end>=rx[1]&bin.start<=rx[2])]
    if(max(1,min(idx)-expand)<=min(idx)-1)
        idx=c(seq(max(1,min(idx)-expand),min(idx)-1,by=1),idx)
    if(max(idx)+1<=min(max(idx)+expand,nrow(hic.bins.chr[chr==ch&is.finite(weight)])))
        idx=c(idx,seq(max(idx)+1,min(max(idx)+expand,nrow(hic.bins.chr[chr==ch&is.finite(weight)])),by=1))
    hic.bins.x=hic.bins.chr[chr==ch&is.finite(weight)][idx]
    idx=hic.bins.chr[chr==ch&is.finite(weight)][,which(bin.end>=ry[1]&bin.start<=ry[2])]
    if(max(1,min(idx)-expand)<=min(idx)-1)
        idx=c(seq(max(1,min(idx)-expand),min(idx)-1,by=1),idx)
    if(max(idx)+1<=min(max(idx)+expand,nrow(hic.bins.chr[chr==ch&is.finite(weight)])))
        idx=c(idx,seq(max(idx)+1,min(max(idx)+expand,nrow(hic.bins.chr[chr==ch&is.finite(weight)])),by=1))
    hic.bins.y=hic.bins.chr[chr==ch&is.finite(weight)][idx]
    hic.data.tmp=hic.data.chr[CJ(bin.index1=hic.bins.x[,bin.index],bin.index2=hic.bins.y[,bin.index]),on=.(bin.index1,bin.index2)]
    ##set NA to 0
    hic.data.tmp[is.na(count.normalized),count.normalized:=0]

    ##add chr, bin.start, bin.end
    setkey(hic.bins.x,bin.index)
    setkey(hic.data.tmp,bin.index1)
    hic.data.tmp[,c("bin.start1","bin.end1"):=hic.bins.x[hic.data.tmp,.(bin.start,bin.end)]]
    setkey(hic.bins.y,bin.index)
    setkey(hic.data.tmp,bin.index2)
    hic.data.tmp[,c("bin.start2","bin.end2"):=hic.bins.y[hic.data.tmp,.(bin.start,bin.end)]]

    ##add smoothed mean count at corresponding distance and eval score=obs/exp
    hic.data.tmp[,count.normalized.mean.smooth:=count.vs.distance[chr==ch][hic.data.tmp[,.(distance=abs(bin.start2-bin.start1))],count.normalized.mean.smooth,on=.(distance)]]
    hic.data.tmp[,score:=(count.normalized/count.normalized.mean.smooth)]
    ##set diagonal to 0
    hic.data.tmp[bin.index1==bin.index2,score:=0]
    
    ##convert to image
    hic.data.tmp[,bin.index1:=factor(bin.index1,levels=hic.bins.x[,bin.index])]
    hic.data.tmp[,bin.index2:=factor(bin.index2,levels=hic.bins.y[,bin.index])]
    ##image with score (i.e. obs/exp)
    image_score=data.table::dcast(hic.data.tmp,bin.index1~bin.index2,value.var="score",fill=0,drop=FALSE)
    image_score=as.Image(unname(as.matrix(image_score[,-1])))
    ##image with count.normalized 
    image_countnorm=data.table::dcast(hic.data.tmp,bin.index1~bin.index2,value.var="count.normalized",fill=0,drop=FALSE)
    image_countnorm=as.Image(unname(as.matrix(image_countnorm[,-1])))
    
    ##gaussian blur
    r=2*ceiling(3*blur.sigma)+1
    if(r>min(dim(image_score)))
    {
        r=2*ceiling(min(dim(image_score))/2)-1
    }
    image_score_blur=gblur(image_score,blur.sigma,radius=r)
    image_countnorm_blur=gblur(image_countnorm,blur.sigma,radius=r)
    
    mask=watershed(image_score_blur-min(image_score_blur)+1,tolerance=watershed.tolerance,ext=1) ##0 is interpreted as background

    ##find min/max value in each watershed region as well as position of the max
    ids=sort(unique(as.vector(mask)))
    ids=ids[ids>0]
    if(length(ids)==0)return(NULL)
    tmp=data.table(reshape2::melt(imageData(mask),varnames=c("x","y"),value.name="id"))
    tmp[,value:=reshape2::melt(imageData(image_score_blur),varnames=c("x","y"),value.name="value")$value]
    tmp=tmp[,.(x.max=x[which.max(value)],y.max=y[which.max(value)]),by=id]

    ##adjust position of the max to match position of highest non-blurred pixel within radius blur sigma
    m=CJ(dx=seq(-ceiling(blur.sigma),ceiling(blur.sigma),by=1),dy=seq(-ceiling(blur.sigma),ceiling(blur.sigma),by=1))[sqrt(dx^2+dy^2)<=blur.sigma*1.5]
    get.max=function(x,y){
        p=m[,.(x=x+dx,y=y+dy)][x>=1&y>=1&x<=nrow(image_score)&y<=ncol(image_score)];
        ind=which.max(image_score[as.matrix(p)]);
        return(p[ind])}
    tmp[,c("x.max","y.max"):=get.max(x.max,y.max),by=.(x.max,y.max)]
    
    setorder(tmp,id)
    tmp=tmp[id>0]
    tmp[,score.max:=image_score_blur[as.matrix(tmp[,.(x.max,y.max)])]]

    ##count normalized
    tmp[,count.normalized.max:=image_countnorm_blur[as.matrix(tmp[,.(x.max,y.max)])]]

    ##filter 1
    tmp=tmp[score.max>threshold_score] 
    ##filter 2
    tmp=tmp[score.max>max(score.max)/loop.filter.param]
    ##filter 3
    threshold_count=count.vs.distance[chr==ch&distance==round(threshold_count_dist/binsize)*binsize,count.normalized.mean.smooth]
    tmp=tmp[count.normalized.max>threshold_count]
    
    ##add coordinates
    tmp=cbind(tmp,
              hic.bins.x[tmp$x.max,.(chr1=chr,start1=bin.start,end1=bin.end)],
              hic.bins.y[tmp$y.max,.(chr2=chr,start2=bin.start,end2=bin.end)])

    ##filter meta-loops outside of region
    tmp[,keep:=FALSE]
    for(j in 1:nrow(regions.tmp))
    {
        tmp[start1>=regions.tmp[j,start1]-binsize&end1<=regions.tmp[j,end1]+binsize&
            start2>=regions.tmp[j,start2]-binsize&end2<=regions.tmp[j,end2]+binsize,keep:=TRUE]
    }
    tmp=tmp[keep==TRUE]
    
    metaloops=rbind(metaloops,
                    tmp[,.(chr1,start1,end1,chr2,start2,end2)])
}
if(is.null(metaloops))
    metaloops=data.table(chr1=character(0),start1=numeric(0),end1=numeric(0),chr2=character(0),start2=numeric(0),end2=numeric(0))
##remove duplicate
metaloops=unique(metaloops[,.(chr1,start1=as.integer(start1),end1=as.integer(end1),chr2,start2=as.integer(start2),end2=as.integer(end2))])

cat("Found ",nrow(metaloops)," meta-loops\n",sep="")

cat("creating",outputfile,"\n")
fwrite(metaloops[order(chr1,start1,end1,chr2,start2,end2)],outputfile,sep="\t",col.names=TRUE)

cat("Done\n")

