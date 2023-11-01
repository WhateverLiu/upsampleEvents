

num2str <- function(x, Nround, regularSpace = F, fillLeft = T, Nchar = 5)
{
  y = paste0(round(x, Nround), "")
  if(!grepl("[.]", y) & Nround > 0)
    rst = paste0(y, ".", paste0(rep("0", Nround), collapse = ""))
  else
  {
    Nzero = max(0L, Nround - nchar(strsplit(y, "[.]")[[1]][2]))
    rst = paste0(y, paste0(rep("0", Nzero), collapse = ""))
  }
  if(regularSpace)
  {
    nfill = max(0L, Nchar - nchar(rst))
    theFill = paste0(rep(" ", nfill), collapse = "")
    if (fillLeft) rst = paste0(theFill, rst)
    else rst = paste0(rst, theFill)
  }; rst
}
