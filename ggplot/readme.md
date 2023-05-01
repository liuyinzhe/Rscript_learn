
####R语言ggplot2绘图显示不全

- 提高分辨率：ggsave函数，指定分辨率来保存图像

 ```
 ggsave("top10_feature_importance.png",plot=p, height=8, width=12, units="in",dpi=300)
 ggsave("top10_feature_importance.pdf",plot=p, height=8, width=12, units="in",dpi=300)
 ```
 
- 调整图像大小：coord_cartesian函数，调整图像的大小
- 
```
coord_cartesian(
  xlim = NULL,
  ylim = NULL,
  expand = TRUE,
  default = FALSE,
  clip = "on"
)
```

- 调整坐标轴限制：scale_x_continuous和scale_y_continuous


