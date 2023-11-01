
load("data/eventLosses.RData")
load("data/yearEvent.RData")
load("data/countyInfo.RData")
load("data/source/stateBorderLine.RData")
load("data/regionIDs.RData")
load("data/eventTlosses.RData")
load("data/trueEPs.RData")
source('code/R/rfuns.R')
invisible(sapply(list.files('code/CharlieLib/R', full.names = T), 
                 function(x) source(x)))
eventMapping = data.table::fread('data/source/M21_upsampled_100K.m__') # Previous work.


# Plot loss footprint.
if (T)
{
  
  
  colnames(eventMapping) = c('original', 'mapped')
  xlim = range(countyInfo$lon)
  ylim = range(countyInfo$lat)
  clrs = unique(colorRampPalette(c(
    'white', 'yellow', 'orange', 'red', 'darkred', 'black'))(10000))
  euc = new.env()
  l1 = new.env()
  lhalf = new.env()
  load("data/Euclidean.RData", envir = euc)
  load("data/L1.RData", envir = l1)
  load("data/L-half.RData", envir = lhalf)
  cutOff = which(eventTlosses$year > 10000L)[1]
  
  
  if (T)
  {
    
    set.seed(42)
    e1iV = sample(cutOff:length(eventLosses), 100)
    u = 1L
    while (u <= length(e1iV)) 
    {
      e1i = e1iV[u]
      eventID = eventTlosses$event[e1i]
      e2i = which(eventMapping$mapped[e1i] == eventMapping$original)
      previousWorkReplacedID = eventTlosses$event[e2i]
      e3i = which(euc$emap$mapped[e1i] == eventMapping$original)
      eucID = eventTlosses$event[e3i]
      e4i = which(l1$emap$mapped[e1i] == eventMapping$original)
      l1ID = eventTlosses$event[e4i]
      e5i = which(lhalf$emap$mapped[e1i] == eventMapping$original)
      lhalfID = eventTlosses$event[e5i]
      
      
      eis = c(e1i, e2i, e3i, e4i, e5i)
      eiIDs = c(eventID, previousWorkReplacedID, eucID ,l1ID, lhalfID)
      tils = c('Truth', 'Previous work', 'New, Euclidean', 'New, L1', 'New, L-0.5') 
      X = lapply(eis, function(x)
      {
        getEventLossFootprint(
          eventLosses[[x]], countyInfo = countyInfo, regionIDs = regionIDs)
      })
      
      
      polygonColors = getColors( lapply(X, function(x) x$loss), clrs, method = 'gau' )
      
      
      dname = 'figure/lossFootprint5eventsOnePic'
      dir.create(dname, showWarnings = F)
      fname = paste0('e-', paste0(eiIDs, collapse = '-'), '-.png')
      fname = paste0(dname, '/', fname)
      
      
      png(filename = fname, width = 16, height = 16 * 0.45, unit = "in", res = 120)
      layout(rbind(c(1, 2, 2, 3, 3, 4),
                   c(5, 5, 6, 6, 7, 7)))
      par(mar = c(1, 1, 1, 1), family = "serif")
      
      
      plot(x = numeric(length(clrs)) + 0.25, y = seq(0,1,len = length(clrs)), col = clrs, 
           pch = 15, bty = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", 
           cex = 2, xlim = c(0, 1), ylim = c(0, 1))
      text(x = 0.5, y = 0.05, labels = "Low loss", cex = 2)
      text(x = 0.5, y = 0.95, labels = "High loss", cex = 2)
      
      
      for (k in 1:length(eis))
      {
        # if (k == 1) plot.new()
        e = X[[k]]
        plot(0, xlim = xlim, ylim = ylim, bty = "n", xaxt = 'n', yaxt = 'n', 
             xaxs = "i", yaxs = "i")
        invisible(lapply(bdLine, function(x) lines(x, col = scales::alpha('blue', 0.5))))
        i = 1L
        invisible(lapply(e$polygon, function(x)
        {
          polygon(x[,1], x[,2], border = NA, col = polygonColors[[1]][i])
          i <<- i + 1L
        }))
        
        
        pc = plotCoor(x = 0.7, y = 0.25)
        text(x = pc[1], y = pc[2], labels = tils[k], cex = 2, pos = 4)
        pc = plotCoor(x = 0.25, y = 0.05)
        text(x = pc[1], y = pc[2], labels = paste0(c(paste0('Event ', eiIDs[k]),
          paste0('Total loss = ', simpDollar(eventTlosses$loss[eis[k]]))), 
          collapse = '\n'), cex = 2, pos = 4)
        
        
        if (k == 2) plot.new()
      }  
      dev.off()
      
      
      if (file.info(fname)$size > 100) u = u + 1L
      else unlink(fname, force = T)
    }
    
    
    tmp = list.files('figure/lossFootprint5eventsOnePic',full.names = T)
    tmp = tmp[tools::file_ext(tmp) == 'png']
    putFirst = as.integer(c(723953, 797622, 1064481, 417818, 1041464, 
                            1327685, 1343451, 751587, 1089792, 1737467, 2327042))
    tmp2 = as.integer(unlist(lapply(strsplit(basename(tmp), '-'), function(x) x[2])))
    tmp2 = match(putFirst, tmp2)
    tmp = tmp[c(tmp2, setdiff(1:length(tmp), tmp2))]
    magick::image_write(
      magick::image_read(tmp), 
      "figure/lossFootprint5eventsOnePic.pdf", format = "pdf") 
  }
  
  
}




# Plot loss footprint, Euclidean, L1, Euclidean on expanded dim, i.e. 
#   with extra dimensions that are nearest 8, 16, 32, 64 counties for each
#   county.
if (T)
{
  colnames(eventMapping) = c('original', 'mapped')
  xlim = range(countyInfo$lon)
  ylim = range(countyInfo$lat)
  clrs = unique(colorRampPalette(c(
    'white', 'yellow', 'orange', 'red', 'darkred', 'black'))(10000))
  euc = new.env()
  l1 = new.env()
  eucDimExt = new.env()
  L1DimExt = new.env()
  load("data/Euclidean.RData", envir = euc)
  load("data/L1.RData", envir = l1)
  # load("data/L-half.RData", envir = lhalf)
  load("data/EuclideanExpanded-8-16-32-64.RData", envir = eucDimExt)
  load("data/L1Expanded-8-16-32-64.RData", envir = L1DimExt)
  cutOff = which(eventTlosses$year > 10000L)[1]
  
  
  if (T)
  {
    
    #   unlink('figure/lossFootprint6eventsOnePic-dimExtended', recursive = T); dir.create('figure/lossFootprint6eventsOnePic-dimExtended', showWarning = F)
    
    
    set.seed(42)
    e1iV = sample(cutOff:length(eventLosses), 100)
    u = 1L
    
    
    while (u <= length(e1iV))
    {
      e1i = e1iV[u]
      eventID = eventTlosses$event[e1i]
      
      
      e2i = which(eventMapping$mapped[e1i] == eventMapping$original)
      previousWorkReplacedID = eventTlosses$event[e2i]
      
      
      e3i = which(euc$emap$mapped[e1i] == eventMapping$original)
      eucID = eventTlosses$event[e3i]
      
      
      e4i = which(l1$emap$mapped[e1i] == eventMapping$original)
      l1ID = eventTlosses$event[e4i]
      
      
      e5i = which(eucDimExt$emap$mapped[e1i] == eventMapping$original)
      eucExtID = eventTlosses$event[e5i]
      
      
      e6i = which(L1DimExt$emap$mapped[e1i] == eventMapping$original)
      L1ExtID = eventTlosses$event[e6i]
      
      
      eis = c(e1i, e2i, e3i, e5i, e4i, e6i)
      eiIDs = c(eventID, previousWorkReplacedID, eucID, eucExtID ,l1ID, L1ExtID)
      tils = c('Truth', 'Previous work', 'New, Euclidean', 
               'New, Euclidean\ndim extended', 'New, L1', 
               'New, L1\ndim extended')
      X = lapply(eis, function(x)
      {
        getEventLossFootprint(
          eventLosses[[x]], countyInfo = countyInfo, regionIDs = regionIDs)
      })
      
      
      polygonColors = getColors( lapply(X, function(x) x$loss), clrs, method = 'gau' )
      
      
      dname = 'figure/lossFootprint6eventsOnePic-dimExtended'
      dir.create(dname, showWarnings = F)
      fname = paste0('e-', paste0(eiIDs, collapse = '-'), '-.png')
      fname = paste0(dname, '/', fname)
      
      
      png(filename = fname, width = 20, height = 20 * 0.33, unit = "in", res = 120)
      # layout(rbind(c(1, 2, 2, 3, 3, 4, 4, 5),
      #              c(6, 6, 7, 7, 8, 8, 9, 9)))
      layout(rbind(c(1, 2, 3, 4),
                   c(5, 6, 7, 8)))
      par(mar = c(1, 1, 1, 1), family = "serif")
      
      
      plot(x = numeric(length(clrs)) + 0.25, y = seq(0,1,len = length(clrs)),
           col = clrs, pch = 15, bty = "n", xlab = "", ylab = "",
           xaxt = "n", yaxt = "n", cex = 2, xlim = c(0, 1), ylim = c(0, 1))
      text(x = 0.4, y = 0.05, labels = "Low loss", cex = 2)
      text(x = 0.4, y = 0.95, labels = "High loss", cex = 2)
      
      # u = 1L
      # for (k in 1:length(eis))
      k = 1L
      for (uu in 2:8)
      {
        # print(k)
        if (uu %in% c(4) ) { plot.new(); next }
        # cat(uu, ', ')
        e = X[[k]]
        plot(0, xlim = xlim, ylim = ylim, bty = "n", xaxt = 'n', yaxt = 'n', 
             xaxs = "i", yaxs = "i")
        invisible(lapply(bdLine, function(x) lines(x, col = scales::alpha('blue', 0.5))))
        i = 1L
        invisible(lapply(e$polygon, function(x)
        {
          polygon(x[,1], x[,2], border = NA, col = polygonColors[[1]][i])
          i <<- i + 1L
        }))
        
        
        pc = plotCoor(x = 0.65, y = 0.20)
        text(x = pc[1], y = pc[2], labels = tils[k], cex = 2, pos = 4)
        
        
        pc = plotCoor(x = 0.24, y = 0.05)
        text(x = pc[1], y = pc[2], labels = paste0(
          c(paste0('Event ', eiIDs[k]), 
            paste0('Total loss = ', simpDollar(eventTlosses$loss[eis[k]]))), 
          collapse = '\n'), cex = 2, pos = 4)
        
        
        k = k + 1L
      }
      
      
      dev.off()
      
      
      if (file.info(fname)$size > 100) u = u + 1L
      else unlink(fname, force = T)
    }
    
    
    
    
    tmp = list.files('figure/lossFootprint6eventsOnePic-dimExtended', full.names = T)
    tmp = tmp[tools::file_ext(tmp) == 'png']
    putFirst = as.integer(c(723953, 797622, 1064481, 417818, 1041464, 
                            1327685, 1343451, 751587, 1089792, 1737467, 2327042))
    tmp2 = as.integer(unlist(lapply(strsplit(basename(tmp), '-'), function(x) x[2])))
    tmp2 = match(putFirst, tmp2)
    tmp = tmp[c(tmp2, setdiff(1:length(tmp), tmp2))]
    magick::image_write(magick::image_read(tmp), 
      "figure/lossFootprint6eventsOnePic-dimExtended.pdf", format = "pdf")
    
    
  }
  
  
}




# Plot loss footprint, Euclidean, L1, Euclidean on expanded dim, i.e. 
#   with extra dimensions that are nearest 8, 16, 32, 64 counties for each
#   county. But, no country and state dimensions.
if (T)
{
  
  
  colnames(eventMapping) = c('original', 'mapped')
  xlim = range(countyInfo$lon)
  ylim = range(countyInfo$lat)
  clrs = unique(colorRampPalette(c(
    'white', 'yellow', 'orange', 'red', 'darkred', 'black'))(10000))
  euc = new.env()
  l1 = new.env()
  eucDimExt = new.env()
  L1DimExt = new.env()
  load("data/Euclidean.RData", envir = euc)
  load("data/L1.RData", envir = l1)
  # load("data/L-half.RData", envir = lhalf)
  load("data/EuclideanExpanded-8-16-32-64-noCountryStates.RData", envir = eucDimExt)
  load("data/L1Expanded-8-16-32-64-noCountryStates.RData", envir = L1DimExt)
  cutOff = which(eventTlosses$year > 10000L)[1]
  
  
  if (T)
  {
    
    #   unlink('figure/lossFootprint6eventsOnePic-dimExtended-noCountryStates', recursive = T); dir.create('figure/lossFootprint6eventsOnePic-dimExtended-noCountryStates', showWarning = F)
    
    
    set.seed(42)
    e1iV = sample(cutOff:length(eventLosses), 100)
    u = 1L
    
    
    while (u <= length(e1iV))
    {
      e1i = e1iV[u]
      eventID = eventTlosses$event[e1i]
      
      
      e2i = which(eventMapping$mapped[e1i] == eventMapping$original)
      previousWorkReplacedID = eventTlosses$event[e2i]
      
      
      e3i = which(euc$emap$mapped[e1i] == eventMapping$original)
      eucID = eventTlosses$event[e3i]
      
      
      e4i = which(l1$emap$mapped[e1i] == eventMapping$original)
      l1ID = eventTlosses$event[e4i]
      
      
      e5i = which(eucDimExt$emap$mapped[e1i] == eventMapping$original)
      eucExtID = eventTlosses$event[e5i]
      
      
      e6i = which(L1DimExt$emap$mapped[e1i] == eventMapping$original)
      L1ExtID = eventTlosses$event[e6i]
      
      
      eis = c(e1i, e2i, e3i, e5i, e4i, e6i)
      eiIDs = c(eventID, previousWorkReplacedID, eucID, eucExtID ,l1ID, L1ExtID)
      tils = c('Truth', 'Previous work', 'New, Euclidean', 
               'New, Euclidean\ndim extended\nCountry,States\nexcluded', 'New, L1', 
               'New, L1\ndim extended\nCountry,States\nexcluded')
      X = lapply(eis, function(x)
      {
        getEventLossFootprint(
          eventLosses[[x]], countyInfo = countyInfo, regionIDs = regionIDs)
      })
      
      
      polygonColors = getColors( lapply(X, function(x) x$loss), clrs, method = 'gau' )
      
      
      dname = 'figure/lossFootprint6eventsOnePic-dimExtended-noCountryStates'
      dir.create(dname, showWarnings = F)
      fname = paste0('e-', paste0(eiIDs, collapse = '-'), '-.png')
      fname = paste0(dname, '/', fname)
      
      
      png(filename = fname, width = 20, height = 20 * 0.33, unit = "in", res = 120)
      # layout(rbind(c(1, 2, 2, 3, 3, 4, 4, 5),
      #              c(6, 6, 7, 7, 8, 8, 9, 9)))
      layout(rbind(c(1, 2, 3, 4),
                   c(5, 6, 7, 8)))
      par(mar = c(1, 1, 1, 1), family = "serif")
      
      
      plot(x = numeric(length(clrs)) + 0.25, y = seq(0,1,len = length(clrs)),
           col = clrs, pch = 15, bty = "n", xlab = "", ylab = "",
           xaxt = "n", yaxt = "n", cex = 2, xlim = c(0, 1), ylim = c(0, 1))
      text(x = 0.4, y = 0.05, labels = "Low loss", cex = 2)
      text(x = 0.4, y = 0.95, labels = "High loss", cex = 2)
      
      # u = 1L
      # for (k in 1:length(eis))
      k = 1L
      for (uu in 2:8)
      {
        # print(k)
        if (uu %in% c(4) ) { plot.new(); next }
        # cat(uu, ', ')
        e = X[[k]]
        plot(0, xlim = xlim, ylim = ylim, bty = "n", xaxt = 'n', yaxt = 'n', 
             xaxs = "i", yaxs = "i")
        invisible(lapply(bdLine, function(x) lines(x, col = scales::alpha('blue', 0.5))))
        i = 1L
        invisible(lapply(e$polygon, function(x)
        {
          polygon(x[,1], x[,2], border = NA, col = polygonColors[[1]][i])
          i <<- i + 1L
        }))
        
        
        pc = plotCoor(x = 0.68, y = 0.13)
        text(x = pc[1], y = pc[2], labels = tils[k], cex = 2, pos = 4)
        
        
        pc = plotCoor(x = 0.24, y = 0.05)
        text(x = pc[1], y = pc[2], labels = paste0(
          c(paste0('Event ', eiIDs[k]), 
            paste0('Total loss = ', simpDollar(eventTlosses$loss[eis[k]]))), 
          collapse = '\n'), cex = 2, pos = 4)
        
        
        k = k + 1L
      }
      
      
      dev.off()
      
      
      if (file.info(fname)$size > 100) u = u + 1L
      else unlink(fname, force = T)
    }
    
    
    tmp = list.files('figure/lossFootprint6eventsOnePic-dimExtended-noCountryStates', full.names = T)
    tmp = tmp[tools::file_ext(tmp) == 'png']
    putFirst = as.integer(c(723953, 797622, 1064481, 417818, 1041464, 
                            1327685, 1343451, 751587, 1089792, 1737467, 2327042))
    tmp2 = as.integer(unlist(lapply(strsplit(basename(tmp), '-'), function(x) x[2])))
    tmp2 = match(putFirst, tmp2)
    tmp = tmp[c(tmp2, setdiff(1:length(tmp), tmp2))]
    magick::image_write(
      magick::image_read(tmp), 
      "figure/lossFootprint6eventsOnePic-dimExtended-noCountryStates.pdf", format = "pdf")
    
    
  }
  
  
}




# Plot loss footprint, Euclidean, Jaccard, Euclidean on expanded dim, i.e. 
#   with extra dimensions that are nearest 8, 16, 32, 64 counties for each
#   county. But, no country and state dimensions, and is binary
if (T)
{
  
  
  colnames(eventMapping) = c('original', 'mapped')
  xlim = range(countyInfo$lon)
  ylim = range(countyInfo$lat)
  clrs = unique(colorRampPalette(c(
    'white', 'yellow', 'orange', 'red', 'darkred', 'black'))(10000))
  euc = new.env()
  eucDimExt = new.env()
  jaccardDimExt = new.env()
  load("data/Euclidean.RData", envir = euc)
  load("data/EuclideanExpanded-8-16-32-64-noCountryStates.RData", 
       envir = eucDimExt)
  load("data/EuclideanExpanded-8-16-32-64-noCountryStates-binary.RData", 
       envir = jaccardDimExt)
  # load("data/L1.RData", envir = l1)
  # load("data/L-half.RData", envir = lhalf)
  cutOff = which(eventTlosses$year > 10000L)[1]
  
  
  if (T)
  {
    
    #   unlink('figure/lossFootprint5eventsOnePic-dimExt-Jaccard', recursive = T); dir.create('lossFootprint5eventsOnePic-dimExt-Jaccard', showWarning = F)
    
    
    set.seed(42)
    e1iV = sample(cutOff:length(eventLosses), 100)
    u = 1L
    while (u <= length(e1iV)) 
    {
      e1i = e1iV[u]
      eventID = eventTlosses$event[e1i]
      e2i = which(eventMapping$mapped[e1i] == eventMapping$original)
      previousWorkReplacedID = eventTlosses$event[e2i]
      
      
      e3i = which(euc$emap$mapped[e1i] == eventMapping$original)
      eucID = eventTlosses$event[e3i]
      
      
      e4i = which(eucDimExt$emap$mapped[e1i] == eventMapping$original)
      eucDimExtID = eventTlosses$event[e4i]
      
      
      e5i = which(jaccardDimExt$emap$mapped[e1i] == eventMapping$original)
      jaccardDimExtID = eventTlosses$event[e5i]
      
      
      eis = c(e1i, e2i, e3i, e4i, e5i)
      eiIDs = c(eventID, previousWorkReplacedID, eucID ,eucDimExtID, jaccardDimExtID)
      tils = c('Truth', 'Previous work', 'New, Euclidean', 
               'New, Euclidean\ndim extended\nCountry,States\nExcluded', 
               'New, Jaccard\ndim extended\nCountry,States\nExcluded') 
      X = lapply(eis, function(x)
      {
        getEventLossFootprint(
          eventLosses[[x]], countyInfo = countyInfo, regionIDs = regionIDs)
      })
      
      
      polygonColors = getColors( lapply(X, function(x) x$loss), clrs, method = 'gau' )
      
      
      dname = 'figure/lossFootprint5eventsOnePic-dimExt-Jaccard'
      dir.create(dname, showWarnings = F)
      fname = paste0('e-', paste0(eiIDs, collapse = '-'), '-.png')
      fname = paste0(dname, '/', fname)
      
      
      png(filename = fname, width = 16, height = 16 * 0.45, unit = "in", res = 120)
      layout(rbind(c(1, 2, 2, 3, 3, 4),
                   c(5, 5, 6, 6, 7, 7)))
      par(mar = c(1, 1, 1, 1), family = "serif")
      
      
      plot(x = numeric(length(clrs)) + 0.25, y = seq(0,1,len = length(clrs)), col = clrs, 
           pch = 15, bty = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", 
           cex = 2, xlim = c(0, 1), ylim = c(0, 1))
      text(x = 0.5, y = 0.05, labels = "Low loss", cex = 2)
      text(x = 0.5, y = 0.95, labels = "High loss", cex = 2)
      
      
      for (k in 1:length(eis))
      {
        # if (k == 1) plot.new()
        e = X[[k]]
        plot(0, xlim = xlim, ylim = ylim, bty = "n", xaxt = 'n', yaxt = 'n', 
             xaxs = "i", yaxs = "i")
        invisible(lapply(bdLine, function(x) lines(x, col = scales::alpha('blue', 0.5))))
        i = 1L
        invisible(lapply(e$polygon, function(x)
        {
          polygon(x[,1], x[,2], border = NA, col = polygonColors[[1]][i])
          i <<- i + 1L
        }))
        
        
        pc = plotCoor(x = 0.7, y = 0.15)
        text(x = pc[1], y = pc[2], labels = tils[k], cex = 2, pos = 4)
        pc = plotCoor(x = 0.25, y = 0.05)
        text(x = pc[1], y = pc[2], labels = paste0(c(
          paste0('Event ', eiIDs[k]), paste0(
            'Total loss = ', simpDollar(eventTlosses$loss[eis[k]]))), 
          collapse = '\n'), cex = 2, pos = 4)
        
        
        if (k == 2) plot.new()
      }  
      dev.off()
      
      
      if (file.info(fname)$size > 100) u = u + 1L
      else unlink(fname, force = T)
    }
    
    
    tmp = list.files('figure/lossFootprint5eventsOnePic-dimExt-Jaccard',full.names = T)
    tmp = tmp[tools::file_ext(tmp) == 'png']
    putFirst = as.integer(c(723953, 797622, 1064481, 417818, 1041464, 
                            1327685, 1343451, 751587, 1089792, 1737467, 2327042))
    tmp2 = as.integer(unlist(lapply(strsplit(basename(tmp), '-'), function(x) x[2])))
    tmp2 = match(putFirst, tmp2)
    tmp = tmp[c(tmp2, setdiff(1:length(tmp), tmp2))]
    magick::image_write(
      magick::image_read(tmp), 
      "figure/lossFootprint5eventsOnePic-dimExt-Jaccard.pdf", format = "pdf") 
    
    
  }
  
  
}















