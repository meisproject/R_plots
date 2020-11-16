library(ggplot2)
library(dplyr)
library(scales)
library(ggnewscale)

# 自己新建一个ggplot2主题
theme_classic_new <- function(base_size = 11, base_family = "", base_line_size = base_size/22, 
                              base_rect_size = base_size/22)
{
  theme_classic(base_size = base_size, base_family = base_family, 
                base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace%
    theme(axis.ticks = element_line(colour = "black", size = 0.8),
          axis.line = element_line(size = 0.8),
          axis.ticks.length = unit(1.5, "mm"),
          axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(size = 20),
          legend.text = element_text(colour = "black", size = 12),
          legend.title = element_text(colour = "black", size = 12))
}

# 基于scales包中的alpha函数新建一个明度函数
luminance <- function (colour, luminance = NA) 
{
  if (length(colour) != length(luminance)) {
    if (length(colour) > 1 && length(luminance) > 1) {
      stop("Only one of colour and luminance can be vectorised")
    }
    if (length(colour) > 1) {
      luminance <- rep(luminance, length.out = length(colour))
    }
    else {
      colour <- rep(colour, length.out = length(luminance))
    }
  }
  hcl <- farver::decode_colour(colour, to = "hcl")
  hcl[!is.na(luminance), 2] <- luminance[!is.na(luminance)]
  farver::encode_colour(hcl, from = "hcl")
}


# 自己的工作路径
setwd("C:\\Users\\asus\\Desktop\\1\\test")

# 创建一个示例数据
dat <- data.frame(groupA = ceiling(runif(50, min = 0, max = 30)), groupB = ceiling(runif(50, min = 0, max = 50))) %>% 
  mutate(plus = groupB + groupA) %>% 
  mutate(division = groupB/groupA)

# 分别计算x轴和y轴的比例
x.scale <- -0.05*max(dat$groupA)
y.scale <- -0.05*max(dat$groupB)

# 按照数值，分别赋予不同明度的颜色
x.max <- max(dat$groupA)
y.max <- max(dat$groupB)

dat.l <- filter(dat, division >= 1) %>% 
  mutate(lumi = floor(groupB*100/y.max)) %>% 
  mutate(color = luminance("pink", lumi))
dat.r <- filter(dat, division < 1) %>% 
  mutate(lumi = floor(groupA*100/x.max)) %>% 
  mutate(color = luminance("turquoise", lumi))

dat.m <- rbind(dat.l, dat.r)

# theme_classic_new()是导入的new_theme.R中新建的主题，基于theme_classic()修改
# 去除原有的坐标轴，以annotate的方式，加上新的坐标轴，使得从0开始，不做延伸
my_plot <- ggplot(dat.m, aes(x = groupA, y = groupB)) + theme_classic_new() +
  scale_x_continuous(limits = c(x.scale, max(dat$groupA)*1.1), n.breaks = 8, expand = c(0,0)) + 
  scale_y_continuous(limits = c(y.scale, max(dat$groupB)*1.1), n.breaks = 8, expand = c(0,0)) +
  annotate("segment", x = x.scale, xend = x.scale, y = 0, yend = max(dat$groupB), size = 1.5) +
  annotate("segment", x = 0, xend = max(dat$groupA), y = y.scale, yend = y.scale, size = 1.5) +
  theme(axis.line = element_blank(), legend.position = "none")

my_plot2 <- my_plot +  
  geom_point(aes(color = color), size = 4) + 
  scale_color_manual(values = dat.m$color, limits = dat.m$color) +
  geom_point(data = filter(dat.m, groupA > (x.max/4) & groupB > (y.max/4)), aes(alpha = plus), color = "blue", size = 4)

ggsave(my_plot2, filename = "plot.png", width = 7, height = 7)
       
