override.color<- c('darkgray','darkgray','black','black',
                   'black','darkgray','darkgray','black')
override.shape <- c(4,3,15,16,4,15,16,NA)
override.lty <- c(1,1,1,1,1,1,1,2)
combined_melted_plus_horizontal$is_horizontal = (combined_melted_plus_horizontal$Method== 'FDR = 0.1')

p_1_Combined_h_line_gs = ggplot(combined_melted_plus_horizontal,
                             aes(x=Effect,
                                 y = Value,
                                 color = Method,
                                 shape = Method,
                                 lty = is_horizontal))+
  geom_line(lwd=0.4)+geom_point(size = 1.7)+
  facet_grid(Type~m1,scales = "free_y") +
  theme_bw(base_size = BASE_FONT_SIZE_VERTICAL) +
  xlab(TeX("$\\\\Delta_{ml}$"))+
  xlim(c(0,3))+
  scale_color_manual(values=override.color)+
  scale_shape_manual(values = override.shape)+
  scale_alpha_manual(values = c(rep(0.3,7),1))+
  ylab("")+ 
  theme(panel.spacing = unit(0, "lines")) +
  scale_linetype(guide='none') +
  #scale_shape(guide = 'none') +
  #scale_color(guide = 'none') +
  guides(color = guide_legend(override.aes = list(shape = override.shape,
                                                    linetype = override.lty,
                                                    color = override.color)))
p_1_Combined_h_line_gs
ggsave(p_1_Combined_h_line_gs,filename = paste0('../../Results/AOAS_combined_greyscale.pdf'),width = 7*1.5,height = 4.5)
