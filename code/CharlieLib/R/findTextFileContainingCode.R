

findTextFileContainingCodeLine <- function(
  direct, codeLine, filetype = c("R", "r"),
  findAll = F, recursive = T, readCharNchar = 1e7L)
{
  flist = list.files(direct, full.names = T, recursive = recursive)
  ext = tools::file_ext(flist)
  # print(ext)
  flist = flist[ext %in% filetype]
  if(findAll)
  {
    rst = unlist(lapply(flist, function(x)
    {
      tmp = readChar(x, nchars = readCharNchar)
      grepl(codeLine, tmp)
    })); gc()
    flist[rst]
  }
  else
  {
    found = F
    for(x in flist)
    {
      tmp = readChar(x, nchars = readCharNchar)
      if(grepl(codeLine, tmp))
      {
        found = T
        break
      }
    }
    if(found) x else ""
  }
}
















