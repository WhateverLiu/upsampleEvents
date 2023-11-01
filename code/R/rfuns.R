

# ==============================================================================
# eventMapping is a 2-col data.frame.
# ==============================================================================
generateEPs = function(eventMapping, regionIDs, eventLosses, yearEvent, 
                       useCppAccelerator = TRUE)
{
  origEvents = eventMapping[[1]]
  if (length(origEvents) != length(yearEvent$event))
    stop("Event mapping file wrong.")
  if (!all(eventMapping[[1]] == yearEvent$event))
    stop("Event mapping file wrong.")
  mappedEvents = eventMapping[[2]]
  newLosses = eventLosses[match(mappedEvents, yearEvent$event)]
  
  
  if (useCppAccelerator ) 
    makeEPs(newLosses, yearEvent$year)
  else
  {
    NregionsHitByEvent = unlist(lapply(newLosses, function(x) length(x[[1]])))
    yearIDs = rep(yearEvent$year, NregionsHitByEvent)
    
    
    dat = data.table::data.table(
      year = yearIDs,
      region = unlist(lapply(newLosses, function(x) x[[1]])),
      loss = unlist(lapply(newLosses, function(x) x[[2]])))
    
    
    aggLoss = dat[, .(loss = sum(loss)), by = .(region, year)][
      order(region, -loss), c('region', 'loss')][
        , .(loss = list(loss)), by = region]$loss
    
    
    occLoss = dat[, .(loss = max(loss)), by = .(region, year)][
      order(region, -loss), c('region', 'loss')][
        , .(loss = list(loss)), by = region]$loss
    
    
    list(aggLoss = aggLoss, occLoss = occLoss)
  }
  
  
}


# ==============================================================================
# eventLoss is an element in eventLosses.
# ==============================================================================
getEventLossFootprint = function(eventLoss, countyInfo, regionIDs)
{
  rid = regionIDs[eventLoss[[1]] + 1L]
  validInd = rid %% 1000L != 0L
  lss = eventLoss[[2]][validInd]
  data.table::data.table(
    countyInfo[eventLoss[[1]][validInd] + 1L, 
               c('FIPS', 'lon', 'lat', 'polygon'), ], loss = lss)
}


# ==============================================================================
# Get colors
# ==============================================================================
getColors = function(vals, colrs, method = "gau")
{
  if (!is.list(vals))
  {
    if (method == 'gau')
    {
      vals = rank(vals)
      r = vals / (length(vals) + 1)
      qr = qnorm(r)
    }
    else if (method == "rank")
    {
      vals = rank(vals)
      qr = (vals - 1) / (length(vals) - 1)
    }
    else if (method == "original")
    {
      qr = (vals - min(vals)) / diff(range(vals))
    }
    else stop("Coloring method not implemented.")
    ind = as.integer(round((
      qr - min(qr)) / diff(range(qr)) * (length(colrs) - 1L) + 1L))
    colrs[ind]
  }
  else
  {
    ns = cumsum(c(0L, unlist(lapply(vals, function(x) length(x)))))
    clrs = getColors(unlist(vals), colrs, method)
    rst = list()
    for (i in 2:length(ns))
      rst[[length(rst) + 1L]] = clrs[(ns[i-1L] + 1L):ns[i]]
    rst
  }
  
}


epToPlot = function(x, catSize = 100000L)
{
  data.table::data.table(ep = c(1:length(x), catSize) / catSize, loss = x)
}


simpDollar = function(x)
{
  if (x < 1000) paste0('$', round(x))
  else if (x < 1e6) paste0('$', round(x / 1e3), 'K')
  else if (x < 1e9) paste0('$', round(x / 1e6), 'M')
  else if (x < 1e12) paste0('$', round(x / 1e9), 'B')
  else paste0('$', round(x / 1e9), 'T')
}


# Return FIPS of nearest neighbors for each county.
# Return is an integer matrix.
getCountyNN = function(countyFIPS, countyLon, countyLat, K = 64)
{
  D = as.matrix(dist(data.table::data.table(countyLon, countyLon)))
  apply(D, 2, function(x)
  {
    countyFIPS[order(x)[2:(K+1L)]]
  })
}


expandDim = function(eventLosses, countyInfo, newDimStart,
                     NNs = as.integer(c(8, 16, 32, 64)), 
                     takeAvg = TRUE, useCppAccelerator = TRUE)
{
  
  # newDimStart = max(unlist(lapply(eventLosses, function(x) max(x[[1]])))) + 1L
  D = as.matrix(dist(countyInfo[, c('lon', 'lat')]))
  Nmat = apply(D, 2, function(x) 
  {
    countyInfo$countyZbRegionInd[order(x)[1:max(NNs)]]
  })
  
  
  if (useCppAccelerator)
  {
    makeExpandedEvents(
      eventLosses = eventLosses, nmat = Nmat, dimStart = newDimStart, 
      NNs = NNs, takeAvg = takeAvg, maxCore = 10000) 
  }
  else
  {
    f = function(eloss, nmat, dimInd)
    {
      ind = eloss[[1]]
      loss = eloss[[2]]
      if (takeAvg) l = apply(nmat, 2, function(x) { sum(loss[ind %in% x]) / nrow(nmat) })
      else l = apply(nmat, 2, function(x) { sum(loss[ind %in% x]) })
      nonzero = l > 1e-10
      list(dimInd[nonzero], l[nonzero])
    }
    
    
    # for (k in NNs)
    expansion = list()
    for (i in 1:length(NNs))
    {
      k = NNs[i]
      nmat = Nmat[1:k, , drop = F]
      dimInd = newDimStart:(newDimStart + nrow(countyInfo) - 1L) 
      # print(str(dimInd))
      expansion[[i]] = lapply(eventLosses, function(e)
      {
        f(e, nmat, dimInd)
      })
      newDimStart = dimInd[length(dimInd)] + 1L
    }
    
    
    X = c(list(eventLosses), expansion)
    
    
    # list(nmat = Nmat, losses = 
    rst = lapply(1:length(eventLosses), function(i)
    {
      Y = lapply(X, function(x) x[[i]])
      list(unlist(lapply(Y, function(y) y[[1]])), 
           unlist(lapply(Y, function(y) y[[2]])))
    })
    for (i in 1:length(rst)) names(rst[[i]][[2]]) = NULL
    rst
    
    # )    
  }
}








# r = 0.999 means using 0.999 * min distance between counties as the length
#   of the edge of a pixel.
getPixelLen = function(subareaCoords, r = 0.999)
{
  min(dist(subareaCoords)) * r
}


getImgForColoring = function(X, maxPr = 0.995)
{
  ux = sort(unique(X))
  X = matrix(match(X, ux), nrow = nrow(X))
  X = (X - min(X)) / diff(range(X))
  X = X * maxPr + (1 - maxPr) / 2
  # return(X)
  qnorm(X)
}










# # ==============================================================================
# # Find candidate events.
# # poolEventLosses: the total event losses of the first N events. The first
# #   N events serve as replacements for the remaining events. One does not need
# #   to use total loss as the criterion.
# # lossesOfEventsToBeReplaced: the total losses of the remaining M events.
# # Return: a list of size M. rst[[i]] is an integer vector, the 0-based index
# #   vector pointing to the first N events.
# # ==============================================================================
# findCandidates001 = function(poolEventLosses, lossesOfEventsToBeReplaced,
#                          K = 1000L, relativeDiff = TRUE)
# {
#   if (K > length(poolEventLosses)) stop("K is too high.")
#   lapply(lossesOfEventsToBeReplaced, function(x)
#   {
#     err = abs((x - poolEventLosses) / (x + poolEventLosses))
#     order(err)[1:K] - 1L
#   })
# }


codtest = function(X, Y, catSize = 100000L)
{
  yunlist = unlist(Y)
  EY = sum(yunlist) / (catSize * length(Y))
  EY2 = sum(yunlist * yunlist) / (catSize * length(Y))
  rst = sum(mapply(function(x, y)
  {
    x = c(x, numeric(catSize - length(x)))
    y = c(y, numeric(catSize - length(y)))
    sum((x - y) ^ 2)
  }, X, Y)) / (catSize * length(Y))
  1 - rst / (EY2 - EY * EY)
}












