library(data.table)
library(ggplot2)

setwd("~/DocumentsResearch/courses/unix-course/exercises/mouse-go-analysis/")

div = fread("data/01-divergence/divergence.bed")
names(div) = c("gene","num_diff","size","div")

ggplot(div, aes(x=div)) +
  geom_histogram(bins=50)


go = fread("data/02-go/divergence_by_go.txt")
names(go) = c("GOID","DOMAIN","NAME","DIV","NUMGENES")

go[,AVG:=mean(DIV)]
go[,REL_DIV:=DIV-AVG]

p = ggplot(go, aes(x=reorder(NAME,DIV),y=REL_DIV,fill=REL_DIV)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1)) +
  scale_fill_gradient2()

pdf("results/go-enrichment.pdf")
plot(p)
dev.off()
