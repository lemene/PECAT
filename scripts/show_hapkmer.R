#!/usr/bin/env Rscript
library(ggplot2)
args <- commandArgs(T)
print(args)

kmertypes <- unlist(strsplit(args[1], ","))
fnames <- unlist(strsplit(args[2], ","))
methods <- unlist(strsplit(args[3], ","))
output <- args[4]

set.seed(1234)

show.kmers <- function (kmertypes, fnames, methods, output) {
  datas = c()
  for (i in 1:length(fnames)) {
    print(fnames[i])
    datas[[i]] <- read.table(fnames[i]); # name length patkmer matkmer okkmer, errkmer type

    datas[[i]]['V1'] =  as.factor(datas[[i]]['V1'])

    datas[[i]][,3] = datas[[i]][,3] / datas[[i]][,2]
    datas[[i]][,4] = datas[[i]][,4] / datas[[i]][,2]

    t <- factor(c(fnames[i]))

    Method = methods[i]
    datas[[i]] <- cbind(datas[[i]], Method)
    
    datas[[i]] <- datas[[i]][sample(nrow(datas[[i]]), 1000),]
  }

  alldata <- do.call("rbind", datas)
  
  p <- ggplot(alldata, aes(x=V3, y=V4,color=Method)) + geom_point(alpha=0.4) + xlab(kmertypes[1]) + ylab(kmertypes[2]) + 
       theme(text = element_text(family = "Times New Roman",size=20))
  ggsave(p, filename=output)
}


show.kmers3 <- function (kmertypes, fnames, methods, output) {
  datas = c()
  for (i in 1:length(fnames)) {
    print(fnames[i])
    datas[[i]] <- read.table(fnames[i]); # name length patkmer matkmer okkmer, errkmer type

    datas[[i]]['V1'] =  as.factor(datas[[i]]['V1'])

    datas[[i]][,3] = datas[[i]][,3] / datas[[i]][,2]
    datas[[i]][,4] = datas[[i]][,4] / datas[[i]][,2]

    t <- factor(c(fnames[i]))

    Method = methods[i]
    datas[[i]] <- cbind(datas[[i]], Method)
    
    datas[[i]] <- datas[[i]][sample(nrow(datas[[i]]), 1000),]
  }

  alldata <- do.call("rbind", datas)
  
  #p <- ggplot(alldata, aes(x=V3, y=V4,color=Method)) + geom_point(alpha=0.4) + xlab(kmertypes[1]) + ylab(kmertypes[2]) + facet_wrap(~ Method, nrow = 1, ncol = 4) +
  #     coord_fixed(ratio=1) + theme(legend.position="None", text = element_text(family = "Times New Roman",size=12))
  p <- ggplot(alldata, aes(x=V3, y=V4,color=Method)) + 
        geom_point(alpha=0.4) + xlab(kmertypes[1]) + ylab(kmertypes[2]) + 
        facet_wrap(~ Method, nrow=1) +
        coord_fixed(ratio=1) + 
        theme_bw()+ theme( legend.position="None", text = element_text(family = "Arial",size=7),  
                        axis.title=element_text(size=8,face="bold"),
                          strip.background = element_rect(fill = "white", size = 0.5, linetype = "solid", colour = "NA"),
                          strip.text = element_text(size=8)
            )
  ggsave(p, filename=output, width=18.3/2.54)
}

show.kmers3(kmertypes, fnames, methods, output)

show.kmers2 <- function(fname0, fname1) {
  
  data0 <- read.table(fname0)
  data1 <- read.table(fname1)
  
  #data0 <- read.table('arab/result_mecat2_1_100M')
  #data1 <- read.table('arab/result_fsa_100M')
  
  #data0 <- read.table('dro/result_fsa_100M')
  #data1 <- read.table('dro/result_hifi_100M')
  
  #data0 <- read.table('yeast/result_fsa_100M')
  #data1 <- read.table('yeast/result_crr2_100M')
  
  data0['V1'] = as.factor(data0['V1'])
  data1['V1'] = as.factor(data1['V1'])
  
  #type <- factor(c(0))
  type <- factor(c(fname0))
  data0 <- cbind(data0, type)
  #data0 <- data0[sample(nrow(data0), 1000),]
  
  #type <- factor(c(1))
  type <- factor(c(fname1))
  data1 <- cbind(data1, type)
  #data1 <- data1[sample(nrow(data1), 1000),]
  
  data0 = data0[which(data0$V4>data0$V5),]
  data1 = data1[which(data1$V4>data1$V5),]
  
  #data0$V2 = data0$V2/(data0$V2+data0$V3+data0$V4+data0$V5)
  #data0$V3 = data0$V3/(data0$V2+data0$V3+data0$V4+data0$V5)
  #data1$V2 = data1$V2/(data1$V2+data1$V3+data1$V4+data1$V5)
  #data1$V3 = data1$V3/(data1$V2+data1$V3+data1$V4+data1$V5)
  
  for (i in 1:length(data0[,1])) { a = sum(data0[i,2:5]); data0[i,2]=data0[i,2]/a; data0[i,3]=data0[i,3]/a; }
  for (i in 1:length(data1[,1])) { a = sum(data1[i,2:5]); data1[i,2]=data1[i,2]/a; data1[i,3]=data1[i,3]/a; }
  
  alldata <- rbind(data0, data1)
  
  ggplot(alldata, aes(x=V2, y=V3, color=type)) + geom_point(alpha=0.4)
}
