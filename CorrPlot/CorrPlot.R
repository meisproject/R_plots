library(dplyr)
library(ggplot2)

setwd("C:\\Users\\asus\\Desktop\\1\\test2")

geneexp <- read.table("gene_expression.txt", sep = "\t",header = T, fill = TRUE, 
                    quote = NULL, check.names = F, stringsAsFactors = F, comment.char = "")

# get_each_corr()获取两两间的相关系数和p值
get_each_corr <- function(dat_x, dat_y) {
  tmp_corr <- cor.test(dat_x, dat_y)
  cor.dat <- tmp_corr$estimate
  pvalue.dat <- tmp_corr$p.value
  dat.out <- c(cor.dat, pvalue.dat)
  return(dat.out)
}

# 获取Gene_1和其他基因的相关性
# 下面这句apply()里第一个参数选择geneexp中除了基因名和Gene_1外的其他列
# 最后一个参数dat_x选择Gene_1所在列
cor.all <- apply(geneexp[, 3:ncol(geneexp)], 2, get_each_corr, dat_x = geneexp[, 2])
rownames(cor.all) <- c("correlation", "pvalue")

cor.all2 <- data.frame(t(cor.all), Gene = colnames(cor.all))

my_plot <- ggplot(cor.all2) + theme_bw() + geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(aes(x = 0, xend = correlation, y = Gene, yend = Gene)) + 
  geom_point(aes(x = correlation, y = Gene, color = pvalue), size = 6) +
  labs(x = "correlation", y = "") + theme(panel.border = element_blank())

write.table(cor.all2, file = "相关性.txt", sep = "\t", row.names = F, quote = F)
ggsave(my_plot, filename = "相关性.png", width = 5, height = 5)
ggsave(my_plot, filename = "相关性.pdf", width = 5, height = 5)

