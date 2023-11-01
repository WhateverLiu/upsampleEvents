

plotCoor <- function(x, y)
{
  rst = par('usr')
  x = (rst[2] - rst[1]) * x + rst[1]
  y = (rst[4] - rst[3]) * y + rst[3]
  c(x, y)
}








