#!/usr/bin/env Rscript
library(ggplot2)
args <- commandArgs(T)
print(args)

kmertypes <- unlist(strsplit(args[1], ","))
fnames <- unlist(strsplit(args[2], ","))
methods <- unlist(strsplit(args[3], ","))
output <- args[4]

set.seed(1234)


show.kmers3 <- function (kmertypes, fnames, methods, output) {
  datas = c()
  limxy = 0;
  for (i in 1:length(fnames)) {
    print(fnames[i])
    datas[[i]] <- read.table(fnames[i]); # name length patkmer matkmer okkmer, errkmer type

    datas[[i]]['V1'] =  as.factor(datas[[i]]['V1'])

    datas[[i]][,3] = datas[[i]][,3] / datas[[i]][,2]
    datas[[i]][,4] = datas[[i]][,4] / datas[[i]][,2]

    limxy = max(limxy, datas[[i]][,3],datas[[i]][,4] )
    t <- factor(c(fnames[i]))

    Method = methods[i]
    datas[[i]] <- cbind(datas[[i]], Method)
    
  }
  limxy = limxy*1.1

  alldata <- do.call("rbind", datas)
  
  #p <- ggplot(alldata, aes(x=V3, y=V4,color=Method)) + geom_point(alpha=0.4) + xlab(kmertypes[1]) + ylab(kmertypes[2]) + facet_wrap(~ Method, nrow = 1, ncol = 4) +
  #     coord_fixed(ratio=1) + theme(legend.position="None", text = element_text(family = "Times New Roman",size=12))
  p <- ggplot(alldata, aes(x=V3, y=V4,color=Method)) + 
        geom_point(alpha=0.4) + xlab(kmertypes[1]) + ylab(kmertypes[2]) + 
        facet_wrap(~ Method, nrow=1) +
        coord_fixed(ratio=1) + xlim(0, limxy) +  ylim(0, limxy)+
        theme_bw()+ theme( legend.position="None", text = element_text(family = "Arial",size=7),  
                        axis.title=element_text(size=8),
                          strip.background = element_rect(fill = "white", size = 0.5, linetype = "solid", colour = "NA"),
                          strip.text = element_text(size=8)
            )
  ggsave(p, filename=output, width=18, units="cm", height=6)
}

show.kmers3(kmertypes, fnames, methods, output)
