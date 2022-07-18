## do a facet wrap style plot with facet grid style facet labels
## plot data with p1 = facet wrap and p2 = facet grid
## for p1 set   theme(strip.background = element_blank(),
## strip.text = element_blank())
## for p2 set     theme(strip.background = element_blank()
plot_facets <- function(p1,p2){
  library(grid)
  library(gtable) 
  gt1 = ggplot_gtable(ggplot_build(p1))
  gt2 = ggplot_gtable(ggplot_build(p2))
  gt1$grobs[grep('strip-t.+1$', gt1$layout$name)] = gt2$grobs[grep('strip-t', gt2$layout$name)]
  grid.draw(gt1)
  
  gt.side1 = gtable_filter(gt2, 'strip-r-1')
  gt.side2 = gtable_filter(gt2, 'strip-r-2')
  #gt.side3 = gtable_filter(gt2, 'strip-r-3')
  
  gt1 = gtable_add_cols(gt1, widths=gt.side1$widths[1], pos = -1)
  gt1 = gtable_add_grob(gt1, zeroGrob(), t = 1, l = ncol(gt1), b=nrow(gt1))
  
  panel_id <- gt1$layout[grep('panel-.+1$', gt1$layout$name),]
  gt1 = gtable_add_grob(gt1, gt.side1, t = panel_id$t[1], l = ncol(gt1))
  gt1 = gtable_add_grob(gt1, gt.side2, t = panel_id$t[2], l = ncol(gt1))

  return(gt1)
}