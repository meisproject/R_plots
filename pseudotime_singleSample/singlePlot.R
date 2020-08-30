# 因为monocle的拟时序图会把所有样本的结果以分面的形式展示
# 希望能在一张图上单独展示单个样本的结果
# 大致思路是先用monocle的函数画一个图，再用ggbuild格式提取对应的数据，再用ggplot2重新画

library(grid)
library(gtable)
# library(gridExtra)
library(ggplot2)
library(monocle)
library(dplyr)
library(scales)

setwd("C:\\Users\\asus\\Desktop\\1\\test")

seuCDS <- readRDS("Pseudotime.monocle.rds")
gg =  plot_cell_trajectory(seuCDS, color_by = "Cluster") + 
  facet_grid(.~Sample) + 
  theme(axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
        axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))+
  scale_color_discrete(na.value="#FFFFFF00")

# 转化为ggplot_built
ggbuild <- ggplot_build(gg)

# 提取拟时序线的数据
line.data <- ggbuild$data[[1]] %>% 
  select(-PANEL) %>% 
  unique()

base.plot <- ggplot() + 
  geom_segment(data = line.data, aes(x = x, y = y, xend = xend, yend = yend), size = 0.75, linetype = "solid") +
  labs(x = "Component 1", y = "Component 2") + theme_classic() + 
  theme(axis.title = element_text(size = 25), 
        axis.text = element_text(size = 15))

# 提取拟时序节点数据
dot.data <- ggbuild$data[[4]] %>% 
  select(-PANEL) %>% 
  unique()

base.plot <- base.plot 

# 提取拟时序上细胞的数据
point.data <- ggbuild$plot$data

# 定义cluster对应的颜色，保证每个样本图里cluster的颜色一致
hex_codes <- data.frame(Cluster = unique(point.data$Cluster), Color = hue_pal()(length(unique(point.data$Cluster))))
point.data <- left_join(point.data, hex_codes, by = "Cluster") %>% 
  arrange(Cluster)
point.data$Cluster <- factor(point.data$Cluster)
xmax <- ceiling(max(point.data$data_dim_1))
xmin <- floor(min(point.data$data_dim_1))
ymax <- ceiling(max(point.data$data_dim_2))
ymin <- floor(min(point.data$data_dim_2))

for (sample in unique(point.data$Sample)){
  point.data.tmp <- filter(point.data, Sample == sample)
  mono.plot <- base.plot + 
    geom_point(data = point.data.tmp, aes(x = data_dim_1, y = data_dim_2, color = Cluster), size = 1.5) +
    scale_color_manual(values = unique(point.data$Color), labels = unique(point.data$Cluster), drop = F) +
    labs(title = sample) + theme(legend.position = "top") + 
    geom_point(data = dot.data, aes(x = x, y = y), size = 5, colour = "black") + 
    geom_text(data = dot.data, aes(x = x, y = y, label = label), size = 4, colour = "white" ) +
    xlim(xmin, xmax) + ylim(ymin, ymax)
  ggsave(mono.plot, filename = paste0(sample, ".png"), width = 6, height = 6)
  ggsave(mono.plot, filename = paste0(sample, ".pdf"), width = 6, height = 6)
}

