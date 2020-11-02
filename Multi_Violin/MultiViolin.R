library(dplyr)
library(ggplot2)
library(data.table)
library(gtable)
library(getopt)

command = matrix(c('gene_exp','e',1,"character",   # 表达量表格（行为细胞，列为基因）
                   'summary_cell','s',1,"character",       # 第一列细胞，第二列cluster或细胞类型
                   'cluster_order','l',2,"character", # cluster顺序，只要一列
                   'gene_order','g',2,"character",      # 基因顺序，只要一列
                   'cluster_color','c',2,"character", # cluster的颜色，两列，第一列cluster，第二列颜色
                   'title_size','t',1,"numeric", # 基因名的文字大小
                   'legend_width','w',2,"numeric",  # 图片中legend宽度（cm）
                   'each_plot_width','p',2,"numeric",  # 图片中每张violin的宽度（cm）
                   'plot_height','h',2,"numeric",    # 图片高度（cm）
                   'output', 'o', 1, "character")  #输出路径
                 ,byrow = TRUE, ncol = 4)


args = getopt(command)
options(stringsAsFactors = FALSE)

# --------------
# 读取文件
gene_exp <- fread(args$gene_exp, sep = "\t")
colnames(gene_exp)[1] <- "Cell"
gene_exp <- as.data.frame(gene_exp)

summary_cell <- fread(args$summary_cell, sep = "\t")
colnames(summary_cell) <- c("Cell", "Cluster")

# 如果有自定义的gene_order文件，按照对应的基因进行筛选和排序
if (!is.null(args$gene_order)) {
  gene_order <- read.table(args$gene_order, sep = "\t",header = T, fill = TRUE, quote = NULL,
                           check.names = F, stringsAsFactors = F, comment.char = "")
  gene_exp <- gene_exp[, c("Cell", gene_order[, 1])]
}

# ------------------
# 基于theme_classic进行主题的调整
theme_classic_new <- function(base_size = 11, base_family = "", base_line_size = base_size/22, 
                              base_rect_size = base_size/22)
{
  theme_classic(base_size = base_size, base_family = base_family, 
                base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace%
    theme(axis.ticks = element_line(colour = "black", size = 0.8),
          axis.line = element_line(size = 0.8),
          axis.ticks.length = unit(1.5, "mm"),
          axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(size = args$title_size),
          legend.text = element_text(colour = "black", size = 12),
          legend.title = element_text(colour = "black", size = 12))
}

# -----------------
# 根据Seurat包中的Vlnplot写的画图函数
get_each_plot <- function(feature) {
  # 数据校正
  exp_tmp <- gene_exp[, c("Cell", feature)]
  noise <- rnorm(n = length(x = exp_tmp[, -1]))/1e+05
  if (all(exp_tmp[, feature] == exp_tmp[, feature][1])) {
    warning(paste0("All cells have the same value of ", feature, "."))
  } else {
    exp_tmp[, feature] <- exp_tmp[, feature] + noise
  }
  
  # 根据summary_cell表格和gene_exp筛选交集的细胞
  exp_tmp <- inner_join(exp_tmp, summary_cell, by = "Cell") %>% 
    arrange(desc(Cluster))
  
  if (!is.null(args$cluster_order)) {
    cluster_order <- read.table(args$cluster_order, sep = "\t",header = T, fill = TRUE, quote = NULL,
                                check.names = F, stringsAsFactors = F, comment.char = "")
    colnames(cluster_order)[1] <- "Cell"
    exp_tmp$Cluster <- factor(exp_tmp$Cluster, levels = rev(unique(cluster_order[, 1])))
  } else {
    exp_tmp$Cluster <- factor(exp_tmp$Cluster, levels = unique(as.character(unlist(exp_tmp$Cluster))))
  }
  
  
  fill <- "Cluster"
  x <- "Cluster"
  y <- paste0("`", feature, "`")
  
  plot <- ggplot(data = exp_tmp, mapping = aes_string(x = x, y = y, fill = fill)[c(2, 3, 1)]) + 
    geom_violin(scale = "width", adjust = 1, trim = TRUE) + coord_flip() +
    labs(x = "", y = "", title = feature, fill = NULL) + theme_classic_new() + 
    scale_y_continuous(position = "right") +
    theme(plot.title = element_text(hjust = 0.5),  axis.text.y = element_blank(),
          plot.margin = margin(0,0,4,0,"pt"), legend.position = "none")
  
  if (!is.null(args$cluster_color)) {
    cluster_color <- read.table(args$cluster_color, sep = "\t",header = T, fill = TRUE, quote = NULL,
                                check.names = F, stringsAsFactors = F, comment.char = "")
    colnames(cluster_color)[1:2] <- c("Cluster", "Color")
    
    if (!is.null(args$cluster_order)) {
      cluster_final_color <- data.frame(Cluster = rev(cluster_order[, 1])) %>% 
        left_join(cluster_color, by = "Cluster")
    }else {
      cluster_final_color <- data.frame(Cluster = unique(as.integer(unlist(exp_tmp$Cluster)))) %>% 
        left_join(cluster_color, by = "Cluster")
    }
    
    plot <- plot + scale_fill_manual(values = cluster_final_color$Color)
  }
  
  plot2 <- ggplotGrob(plot)
  plot2$widths[3] <- unit(4, "pt")
  
  return(plot2)
}

# ---------------
# plots里存储所有基因的violin图
plots <- lapply(colnames(gene_exp)[-1], get_each_plot)

# --------------
# 利用geom_point()伪造一个legend
dat_for_legend <- inner_join(summary_cell, gene_exp, by = "Cell")
dat_for_legend <- data.frame(Cluster = unique(dat_for_legend$Cluster)) %>% 
  arrange(desc(Cluster))
if (!is.null(args$cluster_order)) {
  cluster_order <- read.table(args$cluster_order, sep = "\t",header = T, fill = TRUE, quote = NULL,
                              check.names = F, stringsAsFactors = F, comment.char = "")
  colnames(cluster_order)[1] <- "Cell"
  dat_for_legend$Cluster <- factor(dat_for_legend$Cluster, levels = rev(unique(cluster_order[, 1])))
} else {
  dat_for_legend$Cluster <- factor(dat_for_legend$Cluster, levels = unique(as.character(unlist(dat_for_legend$Cluster))))
}

plot_legend <- ggplot(dat_for_legend, aes(x = Cluster, y = 1, color = Cluster)) + 
  geom_point(size = 5) + coord_flip() + labs(x = "", y = "", title = "  ") +
  scale_y_continuous(position = "right") + theme_classic_new() + 
  theme(legend.position = "none", axis.line.x = element_line(color = "white"),
        axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = margin(0,0,4,0,"pt"))

if (!is.null(args$cluster_color)) {
  cluster_color <- read.table(args$cluster_color, sep = "\t",header = T, fill = TRUE, quote = NULL,
                              check.names = F, stringsAsFactors = F, comment.char = "")
  colnames(cluster_color)[1:2] <- c("Cluster", "Color")
  
  if (!is.null(args$cluster_order)) {
    cluster_final_color <- data.frame(Cluster = rev(cluster_order[, 1])) %>% 
      left_join(cluster_color, by = "Cluster")
  }else {
    cluster_final_color <- data.frame(Cluster = unique(as.integer(unlist(dat_for_legend$Cluster)))) %>% 
      left_join(cluster_color, by = "Cluster")
  }
  
  plot_legend <- plot_legend + scale_color_manual(values = cluster_final_color$Color)
}

plot_legend <- ggplotGrob(plot_legend)

# ---------------
# 计算legend的宽度
if (!is.null(args$legend_width)) {
  legend_width <- args$legend_width
}else {
  legend_width <- max(nchar(unique(summary_cell$Cluster)))*0.23 + 1.5
}
# 计算每张violin的宽度
if (!is.null(args$each_plot_width)) {
  each_plot_width <- args$each_plot_width
} else {
  each_plot_width <- 2.5
}
# 计算整张图的高度
if (!is.null(args$plot_height)) {
  plot_height <- args$plot_height
} else {
  plot_height <- 0.7*(length(unique(summary_cell$Cluster))) + 1.5
}
# 计算整张图的宽度
plot_width <- max(nchar(unique(summary_cell$Cluster)))*0.25 + 1 + each_plot_width*length(plots)

# ---------------
# 构建一个空的gtable，将所有的图合并到一张图上，并输出
fig_combined <- gtable(widths = unit(c(legend_width, rep(each_plot_width, length(plots))), "cm"),
                       heights = unit(plot_height, "cm"))

fig_combined <- gtable_add_grob(fig_combined, grobs = plot_legend, t = 1, l = 1)

for (i in 1:length(plots)) {
  fig_combined <- gtable_add_grob(fig_combined, grobs = plots[[i]], t = 1, l = i + 1)
}

setwd(args$output)
ggsave(fig_combined, filename = "ViolinPlot.png", width = plot_width + 1, height = plot_height + 1, units = "cm")
ggsave(fig_combined, filename = "ViolinPlot.pdf", width = plot_width + 1, height = plot_height + 1, units = "cm")
