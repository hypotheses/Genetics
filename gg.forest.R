credplot.gg <- function(d){
  # http://www.r-bloggers.com/forest-plots-using-r-and-ggplot2/
  # d is a data frame with 4 columns
  # d$x gives variable names
  # d$y gives center point
  # d$ylo gives lower limits
  # d$yhi gives upper limits
  require(ggplot2)
  p <- ggplot(d, aes(x=x, y=y, ymin=ylo, ymax=yhi))+geom_pointrange()+
    coord_flip() + geom_hline(aes(x=1), lty=2)+ ylab('Odd Ratios')+xlab('Variable')+ scale_y_log()
  return(p)
}
credplot.gg(d)

