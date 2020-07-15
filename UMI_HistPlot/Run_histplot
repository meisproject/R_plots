library(tidyverse)

# 选择工作路径
setwd("C:\\Users\\Desktop\\1")

set.seed(233)

# 构建作图数据
umiinfo <- data.frame(matrix(ncol = 4, nrow = 1450))
colnames(umiinfo) <- c("Cell", "nGene", "nUMI", "mito.percent")
umiinfo[, 1] <- paste0("Cell", seq(1, 1450))
umiinfo[, 2] <- ceiling(runif(1450, min = 1, max = 4523))
umiinfo[, 3] <- ceiling(runif(1450, min = 1, max = 5698))
umiinfo[, 4] <- runif(1450)

# 构建绘图函数
get_HisPlot <- function(need) {
  need.num <- as.numeric(as.matrix(umiinfo[need]))
  dat_median <- median(need.num)
  
  dat <- hist(need.num, breaks=76, plot=FALSE)
  dat <- data.frame(breaks = dat$breaks[-length(dat$breaks)], 
                    counts = dat$counts)
  
  write.table(dat, quote = FALSE, row.names = FALSE, sep = "\t", 
              eol = "\n", file = paste(need, "_区间统计.txt"))
  
  my_plot <- ggplot(dat, aes(x = breaks, y = counts)) + 
    geom_col(fill = "#c6c6c6", orientation = "x") + geom_smooth(colour = "blue", se = F, method = "gam") + 
    geom_vline(xintercept = dat_median, linetype = "dashed", size = 0.8) +
    annotate("text", x = max(dat$breaks), y = max(dat$counts), vjust = 1, hjust = 1,
             label = paste0("Median = ", dat_median), colour = "black", size = 5) + 
    coord_cartesian(expand = F) + xlab(paste0(need, "/Cell")) + theme_classic() +
    theme(axis.title.y = element_blank(), axis.line.y = element_blank(),
          axis.ticks.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.x = element_line(colour = "black", size = 0.8),
          axis.line.x = element_line(size = 0.8),
          axis.ticks.length.x = unit(1.5, "mm"),
          axis.text.x = element_text(colour = "black", size = 13),
          axis.title.x = element_text(size = 20))
  
  ggsave(my_plot, filename = paste0(need, ".png"), width = 5, height = 5)
  ggsave(my_plot, filename = paste0(need, ".pdf"), width = 5, height = 5)
}

type <- colnames(umiinfo)[2:4]
sapply(type, get_HisPlot)
