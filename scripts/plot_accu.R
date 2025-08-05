library(ggplot2)

load_data <- function(fname) {

   data <- read.table(fname)

}
cbpal_4_1 <- c("#377eb8", "#4daf4a", "#ff7f00", "#a65628")

theme_my <- theme(
  legend.position = "right",
  legend.title = element_blank(),
  legend.text = element_text(size=6, family = "Arial"),
  #plot.background = element_blank(),
  text = element_text(size=7, family="Arial"),
  axis.title.y = element_text(size = 6, family = "Arial"),
  axis.title.x = element_text(size = 6, family = "Arial"),
  axis.text=element_text(size=6, family = "Arial"),
  plot.title = element_text(size=8, family = "Arial", hjust=0.5, margin=margin(0,0,5,0)),
  # panel.border = element_blank(),
  strip.background = element_rect(colour="white", fill="white", linewidth=1, linetype="solid"),
  strip.placement="outside",
  strip.text=element_text(size=6, family = "Arial"),
  strip.text.x = element_blank(),
  panel.grid.major = element_line(color = "grey",linewidth = 0.3, linetype = 1),
  panel.grid.minor = element_line(color = "grey",linewidth = 0.3, linetype = 2),
  )

accu0 <- load_data("hg002_ont_60286.txt")
accu1 <- load_data("hg002_ont_106002.txt")
accu <- rbind(accu0, accu1)
accu$V1 <- factor(accu$V1)
print(colnames(accu))
#colnames(accu) = c("Accuracy", "Local minimum accuracy", "Technique")
#colnames(accu) = c("normalized scores")
print(colnames(accu))
png("a.png", width = 5, height = 5, units = "cm", res=500)
#plt <- ggplot(accu,aes(x=V3, fill=V1))+geom_density(alpha=0.2, colour = alpha(0, "black"))+ xlim(0.95,1)
plt <- ggplot(accu, aes(x=V3)) + facet_wrap( ~ V1, ncol=1) + geom_histogram(bins=60) +  theme() + theme_bw() + theme_my +  ylab("Count") + xlab("Normalized score") + xlim(c(-1,1)) + labs(NULL)

plt
dev.off()
