library(tidyverse)
library(ggplot2)

#setwd("~/DocumentsResearch/courses/unix-course/exercises/mouse-go-analysis/")

read_tsv('data/02-go/divergence_by_go.txt', col_names=F) -> go

names(go) <- c("GOID","DOMAIN","NAME","DIV","NUMGENES")

go %>%
  mutate(AVG = mean(DIV)) %>%
  mutate(REL_DIV = DIV-AVG) -> 
  go2

p = ggplot(go2, aes(x=reorder(NAME,DIV),y=REL_DIV,fill=REL_DIV)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1)) +
  scale_fill_gradient2() +
  ylab("Gene divergence") +
  xlab("GO terms") +
  labs(fill = "Relative\ndivergence") +
  coord_flip()
  

pdf("results/go-enrichment.pdf",w=10,h=8)
plot(p)
dev.off()

jpeg("results/go-enrichment.jpg",width=1000*0.8,height=800*0.8)
plot(p)
dev.off()
