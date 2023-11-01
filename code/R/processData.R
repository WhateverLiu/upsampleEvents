

# ==============================================================================
# Process raw input and generate lossInfo.csv.
# ==============================================================================
if (F)
{
  dat = data.table::fread(
    'data/source/TouchstoneRe_Stochastic100kWSST_Catalog_M021_20230503.l__')
  colnames(dat) = c('modelNumber', 'year', 'event', 'day', 'country-FIPS', 
                    'unknown',
                    'lineOfBusiness1loss', 'lineOfBusiness2loss', 
                    'lineOfBusiness3loss', 'lineOfBusiness4loss')
  lossInfo = dat[, c('year', 'event'), ]
  lossInfo$FIPS = as.integer(substring(dat$`country-FIPS`, first = 4))
  lossInfo$loss = dat$lineOfBusiness1loss + dat$lineOfBusiness2loss + 
    dat$lineOfBusiness3loss + dat$lineOfBusiness4loss
  data.table::fwrite(lossInfo, file = "data/lossInfo.csv")
  eventTlosses = lossInfo[FIPS %% 1000L != 0L, list(loss = sum(loss)), 
                          by = list(year, event)]
  save(eventTlosses, file = 'data/eventTlosses.RData')
}


# ==============================================================================
# Read data and generate eventLosses.RData, yearEvent.RData, regionIDs.RData.
# ==============================================================================
if (T)
{
  lossInfo = data.table::fread("data/lossInfo.csv")
  regionIDs = sort(unique(lossInfo$FIPS))
  lossInfo$zbRegionInd = match(lossInfo$FIPS, regionIDs) - 1L # Zero-based index for loss regions.
  lossInfo = lossInfo[order(lossInfo$event, lossInfo$zbRegionInd), , ]
  yearEvent = unique(lossInfo$event)
  yearEvent = data.table::data.table(
    year = lossInfo$year[match(yearEvent, lossInfo$event)], 
    event = yearEvent)
  tmp = data.table::data.table(
    ind = 1:length(lossInfo$event), event = lossInfo$event)[
      , .(ind = list(ind)), by = event]
  eventLosses = lapply(tmp$ind, function(x)
  {
    list(zbRegionInd = lossInfo$zbRegionInd[x], loss = lossInfo$loss[x])
  })
  rm(tmp, lossInfo); gc()
  save(eventLosses, file = "data/eventLosses.RData")
  save(yearEvent, file = "data/yearEvent.RData")
  save(regionIDs, file = "data/regionIDs.RData")
}


# Generate countyInfo.
if (T)
{
  countyInfo = rvest::html_table(rvest::read_html(
    'https://en.wikipedia.org/wiki/User:Michael_J/County_table'))[[1]]
  countyInfo = as.data.frame(countyInfo)
  data.table::setDT(countyInfo)
  countyInfo = countyInfo[, c('FIPS', 'State', 'County [2]', 'Longitude', 'Latitude'), ]
  colnames(countyInfo) = c('FIPS', 'state', 'county', 'lon', 'lat')
  countyInfo$lon = gsub('°', '', countyInfo$lon)
  countyInfo$lon = gsub('–', '-', countyInfo$lon)
  countyInfo$lat = gsub('°', '', countyInfo$lat) 
  countyInfo$lon = as.numeric(countyInfo$lon)
  countyInfo$lat = as.numeric(countyInfo$lat)
  shdat = raster::shapefile('data/cb_2018_us_county_20m/cb_2018_us_county_20m.shp')
  countyPolygons = lapply(shdat@polygons, function(x)
  {
    rst = x@Polygons[[1]]@coords
    dimnames(rst) = NULL; rst
  })
  tmp = sort(intersect(countyInfo$FIPS, as.integer(shdat$GEOID)))
  countyInfo = countyInfo[match(tmp, countyInfo$FIPS), ]
  countyInfo$polygon = countyPolygons[match(tmp, as.integer(shdat$GEOID))]
  tmp = regionIDs[regionIDs %% 1000L != 0L]
  if (length(setdiff(tmp, countyInfo$FIPS)) != 0)
      stop("Some county indices are not present in countyInfo")
  countyInfo = countyInfo[match(tmp, countyInfo$FIPS), ]
  countyInfo$countyZbRegionInd = match(countyInfo$FIPS, sort(unique(lossInfo$FIPS))) - 1L
  save(countyInfo, file = "data/countyInfo.RData")  
  rm(tmp); gc()
}


# Generate state line
if (F)
{
  
  
  shdat = raster::shapefile('data/cb_2018_us_state_20m/cb_2018_us_state_20m.shp')
  statePolygons = lapply(shdat@polygons, function(x)
  {
    rst = x@Polygons[[1]]@coords
    dimnames(rst) = NULL; rst
  })
  
  
  shdat = raster::shapefile('data/cb_2022_us_nation_20m/cb_2022_us_nation_20m.shp')
  tmp = lapply(shdat@polygons[[1]]@Polygons, function(x) {u = x@coords; dimnames(u) = NULL; u})
  statePolygons = c(statePolygons, tmp)
  plot(0, col = "white", xlim = xlim, ylim = ylim)
  invisible( lapply(statePolygons, function(x) lines(x)))
  
  
  
  plot(0, xlim = xlim, ylim = ylim, col = "white")
  invisible(lapply(statePolygons, function(x) lines(
    x, xlim = xlim, ylim = ylim)))
  
  
  
  
  
  
  
  tmp = sort(intersect(countyInfo$FIPS, as.integer(shdat$GEOID)))
  countyInfo = countyInfo[match(tmp, countyInfo$FIPS), , ]
  countyInfo$polygon = countyPolygons[match(tmp, as.integer(shdat$GEOID))]
  tmp = regionIDs[regionIDs %% 1000L != 0L]
  if (length(setdiff(tmp, countyInfo$FIPS)) != 0)
    stop("Some county indices are not present in countyInfo")
  countyInfo = countyInfo[match(tmp, countyInfo$FIPS), , ]
  save(countyInfo, file = "data/countyInfo.RData")  
  rm(tmp); gc()
  
  
}



# Generate true EPs for all loss regions.
if (T)
{
  source('code/R/rfuns.R')
  trueEPs = generateEPs(eventMapping = list(
    yearEvent$event, yearEvent$event), regionIDs = regionIDs, 
    eventLosses = eventLosses, yearEvent = yearEvent)
  save(trueEPs, file = "data/trueEPs.RData") 
}


# Read all saved data. Run some tests.
if (F)
{
  
  
  # LD_PRELOAD="/lib64/libasan.so.5 /lib64/libubsan.so.1"  R
  
  
  load("data/eventLosses.RData")
  load("data/yearEvent.RData")
  load("data/countyInfo.RData")
  load("data/source/stateBorderLine.RData")
  load("data/regionIDs.RData")
  load("data/eventTlosses.RData")
  load("data/trueEPs.RData")
  
  
  eventMapping = data.table::fread('data/source/M21_upsampled_100K.m__')
  
  
  source('code/R/rfuns.R')
  
  
  # trueEPs = generateEPs(eventMapping = data.table::data.table(
  #   yearEvent$event, yearEvent$event), regionIDs = regionIDs, 
  #   eventLosses = eventLosses, yearEvent = yearEvent)
  # save(trueEPs, file = "data/trueEPs.RData")
  
  source("code/CharlieLib/R/CharlieSourceCpp.R")
  CharlieSourceCpp('code/cpp/cppFuns.cpp', 
                   cacheDir = 'tempFiles/CharlieRcpp')
  
  
  
  tmp = generateEPs(eventMapping = eventMapping, regionIDs = regionIDs, 
                    eventLosses = eventLosses, yearEvent = yearEvent)
  countryInd = 1L
  stateInd = which(regionIDs %% 1000L == 0L)[-1]
  countyInd = which(regionIDs %% 1000L != 0L) 
  cod(newEPs = tmp$aggLoss, trueEPs = trueEPs$aggLoss)
  cod(newEPs = tmp$aggLoss[countryInd], trueEPs = trueEPs$aggLoss[countryInd])
  cod(newEPs = tmp$aggLoss[stateInd], trueEPs = trueEPs$aggLoss[stateInd])
  cod(newEPs = tmp$aggLoss[countyInd], trueEPs = trueEPs$aggLoss[countyInd])
  
  
  tmp = sample(length(eventLosses), 1)
  tmp = getEventLossFootprint(
    eventLosses[[tmp]], countyInfo = countyInfo, regionIDs = regionIDs)
  xlim = range(countyInfo$lon)
  ylim = range(countyInfo$lat)
  lss = tmp$loss
  clrs = unique(colorRampPalette(c('white', 'yellow', 'orange',
                                   'red', 'darkred', 'black'))(10000))
  polygonColors = getColors( tmp$loss, clrs, method = 'gau' )
  
  
  png(filename = 'figure/tmp.png', width = 1000, height = 618, unit = "px")
  par(mar = c(0, 0, 0, 0))
  plot(0, xlim = xlim, ylim = ylim, bty = "n", xaxt = 'n', yaxt = 'n', 
       xaxs = "i", yaxs = "i")
  invisible(lapply(bdLine, function(x) lines(x, col = scales::alpha('blue', 0.5))))
  i = 1L
  invisible(lapply(tmp$polygon, function(x)
  {
    polygon(x[,1], x[,2], border = NA, col = polygonColors[i])
    i <<- i + 1L
  }))
  grid(nx = round(diff(xlim) / 4), ny = round(diff(ylim) / 4), col = "black", lty = 2 )
  dev.off()
  
  
  
  
  
  colnames(eventMapping) = c('original', 'mapped')
  xlim = range(countyInfo$lon)
  ylim = range(countyInfo$lat)
  clrs = unique(colorRampPalette(c('white', 'yellow', 'orange',
                                   'red', 'darkred', 'black'))(10000))
  
  # A not so good one: e1i = 143581
  e1i = sample(length(eventLosses), 1)
  e2i = which(eventMapping$mapped[e1i] == eventMapping$original)
  e1 = getEventLossFootprint(
    eventLosses[[e1i]], countyInfo = countyInfo, regionIDs = regionIDs)
  e2 = getEventLossFootprint(
    eventLosses[[e2i]], countyInfo = countyInfo, regionIDs = regionIDs)
  polygonColors = getColors( list(e1$loss, e2$loss), clrs, method = 'gau' )
  
  
  # png(filename = 'figure/tmp2.png', width = 1000, height = 618 * 2, unit = "px")
  par(mar = c(0, 0, 0, 0), mfrow = c(2, 1))
  plot(0, xlim = xlim, ylim = ylim, bty = "n", xaxt = 'n', yaxt = 'n', 
       xaxs = "i", yaxs = "i")
  invisible(lapply(bdLine, function(x) lines(x, col = scales::alpha('blue', 0.5))))
  i = 1L
  invisible(lapply(e1$polygon, function(x)
  {
    polygon(x[,1], x[,2], border = NA, col = polygonColors[[1]][i])
    i <<- i + 1L
  }))
  
  
  plot(0, xlim = xlim, ylim = ylim, bty = "n", xaxt = 'n', yaxt = 'n', 
       xaxs = "i", yaxs = "i")
  invisible(lapply(bdLine, function(x) lines(x, col = scales::alpha('blue', 0.5))))
  i = 1L
  invisible(lapply(e2$polygon, function(x)
  {
    polygon(x[,1], x[,2], border = NA, col = polygonColors[[2]][i])
    i <<- i + 1L
  }))
  # dev.off()
  
  
  
  
  
  
  LD_PRELOAD="/lib64/libasan.so.5 /lib64/libubsan.so.1"  R
  
  
  load("data/eventLosses.RData")
  load("data/yearEvent.RData")
  load("data/countyInfo.RData")
  load("data/source/stateBorderLine.RData")
  load("data/regionIDs.RData")
  load("data/eventTlosses.RData")
  load("data/trueEPs.RData")
  
  
  source('code/R/rfuns.R')
  source("code/CharlieLib/R/CharlieSourceCpp.R")
  
  
  CharlieSourceCpp('code/cpp/cppFuns.cpp', cacheDir = 'tempFiles/CharlieRcpp', 
                   optFlag = '-Ofast', sanitize = F)
  tmp = makeEPs(eventLosses, yearEvent$year, maxCore = 10000)
  
  
  range(unlist(lapply(tmp$aggLoss, function(x) length(x))) - unlist(
    lapply(trueEPs$aggLoss, function(x) length(x))) )
  
  
  data.table::setDT(tmp)
  colnames(tmp) = c('region', 'year', 'loss')
  tmp = tmp[, .(loss = sum(loss)), by = .(region, year)][
    order(-loss), list(loss = list(loss)), by = region]$loss
  
  
  
  
  x = runif(300)
  y = runif(300)
  K = 300L
  tmp = findCandidates(x, y, K = K, relativeDiff = TRUE, maxCore = 1000)
  tmp2 = matrix(unlist(findCandidates001(x, y, K = K, relativeDiff = TRUE)), nrow = K)
  all(tmp == tmp2)
  
  
  cands = findCandidates(eventTlosses$loss[eventTlosses$year <= 10000L],
                         eventTlosses$loss[eventTlosses$year > 10000L], 
                         K = 1000L, relativeDiff = TRUE, maxCore = 1000)
  
  
  i = which(eventTlosses$year > 10000L)[1] - 1L
  tmp = apply(cands, 2, function(x)
  {
    i <<- i + 1L
    eventTlosses$loss[x + 1L] / eventTlosses$loss[i] - 1
    # 2 * (eventTlosses$loss[x + 1L] - eventTlosses$loss[i]) / 
    #       (eventTlosses$loss[x + 1L] + eventTlosses$loss[i])
  })
  range(tmp)
  # This shows cands are more than enough to be the pool for event selection.
  
  
  # tmp = matrix(unlist(findCandidates001(
  #   eventTlosses$loss[eventTlosses$year <= 10000L], 
  #   eventTlosses$loss[eventTlosses$year > 10000L], 
  #   K = 1000L, relativeDiff = TRUE)), nrow = 1000L)
  CharlieSourceCpp('code/cpp/knns.cpp', cacheDir = 'tempFiles/CharlieRcpp', 
                   optFlag = '-Ofast', sanitize = F)
  
  
  lossInfo = data.table::fread('data/lossInfo.csv')
  regionTloss = lossInfo[, list(loss = sum(loss)), by = FIPS]
  regionTloss = regionTloss[order(FIPS), , ]
  regionW = 1.0 / sqrt(regionTloss$loss)
  
  
  regionW = rep(1.0, length(regionIDs))
  regionW[regionIDs %% 1000L != 0L] = 10
  
  
  neis = knns(eventLosses = eventLosses, candidateMat = cands, keep = 100, 
             distance = "L-half", LPp = 0.5,
             maxCore = 1000, verbose = T, 
             regionW = numeric(0)
             # , regionW = regionW
             )
  emap = data.table::data.table(
    original = eventTlosses$event, mapped =  c(
      eventTlosses$event[1:(nrow(eventTlosses) - ncol(neis))],
      eventTlosses$event[neis[1,] + 1L]))
  save(emap, file = "data/L-half.RData")
  newEPs = generateEPs(emap, regionIDs, eventLosses, yearEvent)
  countryInd = 1L
  stateInd = which(regionIDs %% 1000L == 0L)[-1]
  countyInd = which(regionIDs %% 1000L != 0L) 
  cod(newEPs = newEPs$aggLoss, trueEPs = trueEPs$aggLoss)
  cod(newEPs = newEPs$aggLoss[countryInd], trueEPs = trueEPs$aggLoss[countryInd])
  cod(newEPs = newEPs$aggLoss[stateInd], trueEPs = trueEPs$aggLoss[stateInd])
  cod(newEPs = newEPs$aggLoss[countyInd], trueEPs = trueEPs$aggLoss[countyInd])
  
  
# Previous work: *  0.9884478 0.9907376 0.9857348 0.9740161
# Euc:  * 0.9911265 0.991896 0.9902841 0.9796691
# L-0.4: * 0.9781116 0.9725955 0.9788104 0.9709848
# L-0.5:*  0.984943 0.981857 0.9857877 0.9757144
# L-0.6:*  0.9872506 0.9854407 0.9876303 0.9772735
# L-0.7: * 0.9885603 0.9874199 0.9887387 0.9781605
# L-0.75: * 0.989108 0.9882642 0.9892612 0.9782453
# L-0.8:*  0.9895497 0.9889812 0.9895702 0.9785726
# L-0.9:*  0.9901632 0.9899922 0.9899837 0.9790214
# Euc with weights of 1 / sqrt(region total Loss):* 0.9886138 0.987849 0.9880154 0.9795434
# Euc with weights that make county dim 2x: * 0.9904274 0.9908131 0.989506 0.97999
# Euc with weights that make county dim 5x: * 0.9891761 0.9889699 0.9881758 0.9799347
# Euc with weights that make county dim 10x:*  0.9882902 0.9876384 0.9872829 0.9798294
# L1: *  0.9904578 0.9904791 0.990129 0.9794268
# dotProduct:*  0.3099962 0.3407905 0.1888816 -0.2544753
# Cubic root: * 0.9721636 0.9649208 0.972709 0.9653204
# quarter root:* 0.9572002 0.9466973 0.9567128 0.9488642
# sym cross entropy:  *0.8366758 0.8562635 0.7958815 0.6955693
  tmp = readLines('tempFiles/new 1.csv')
  tmp = strsplit(tmp, split = '[*]')
  for (i in 1:length(tmp))
    tmp[[i]][1] = gsub('[#] ', '', tmp[[i]][1])
  for (i in 1:length(tmp))
  {
    tmp2 = as.numeric(strsplit(tmp[[i]][2], split = ' ')[[1]])
    tmp[[i]][2] = paste0(tmp2[!is.na(tmp2)], collapse = ',')
  }
  tmp = unlist(lapply(tmp, function(x) paste0(x, collapse = ',')))
  writeLines(tmp, 'tempFiles/tmp.csv')
  
  
  tmp = t(matrix(paste0(unlist(tmp), collapse = ','), nrow = 2))
  
  
    tmp[[i]][2] = paste0(as.numeric(strsplit(tmp[[i]][2], split = ' |  |   |    |     ')[[1]]), collapse = ',')
  
  
  
  tmp = sapply(tmp, function(x) gsub('[#] ', '', x))
  
  
  
  eucRst[1,]
  
  
  
  
  knnTest = function(X, candidateMat, keep = 5)
  {
    mxi = max(unlist(lapply(X, function(x) max(x[[1]])))) + 1L
    X = lapply(X, function(x)
    {
      y = numeric(mxi)
      y[x[[1]] + 1L] = x[[2]]; y
    })
    i = length(eventLosses) - ncol(candidateMat) + 1L
    apply(candidateMat, 2, function(x)
    {
      e = X[[i]]
      ds = sapply(x, function(u)
      {
        f = X[[u + 1L]]
        sum((e - f) ^ 2)
      })
      i <<- i + 1L
      x[order(ds)[1:keep]]
    })
  }
  
  
  tmp = knnTest(eventLosses, cands[, 200:210], keep = 100)
  tmp2 = knns(eventLosses = eventLosses, candidateMat = cands[, 200:210], 
              keep = 100, distance = "Euclidean", maxCore = 1000, 
              verbose = T) 
  tmp3 = knns(eventLosses = eventLosses, candidateMat = cands[, 200:210], 
              keep = 100, distance = "dotProduct", maxCore = 1000, 
              verbose = T) 
  
  
  
  
  
  source("code/R/rfuns.R")
  tmp = findCandidates(eventTlosses$loss[eventTlosses$year <= 10000L],
                       eventTlosses$loss[eventTlosses$year > 10000L], 
                       K = 1000L, relativeDiff = TRUE)
  save(tmp, file = "tempFiles/tmp.Rdata")
  
  
  source("code/CharlieLib/R/CharlieSourceCpp.R")
  CharlieSourceCpp('code/cpp/cppFuns.cpp', cacheDir = 'tempFiles/CharlieRcpp')
  tmp2 = findCandidates(eventTlosses$loss[eventTlosses$year <= 10000L],
                        eventTlosses$loss[eventTlosses$year > 10000L], 
                        K = 1000L, relativeDiff = TRUE)
  
  
  tmp3 = findCandidates001(eventTlosses$loss[eventTlosses$year <= 10000L],
                        eventTlosses$loss[eventTlosses$year > 10000L], 
                        K = 1000L, relativeDiff = TRUE)
  
  
  
  
  
  for (i in 2000:length(tmp2))
  {
    if (!all(tmp2[[i]] == tmp3[[i]])) break
  }
  str(tmp2[[i]])
  str(tmp3[[i]])
  sum(tmp2[[i]] != tmp3[[i]])
  
  
  
  
  load("data/eventLosses.RData")
  load("data/yearEvent.RData")
  load("data/countyInfo.RData")
  load("data/source/stateBorderLine.RData")
  load("data/regionIDs.RData")
  load("data/eventTlosses.RData")
  load("data/trueEPs.RData")
  lossInfo = data.table::fread('data/lossInfo.csv')
  lossInfoCountyOnly = lossInfo[FIPS %% 1000L != 0L, ]
  lossInfoCountyOnly$countyInd = match(
    lossInfoCountyOnly$FIPS, sort(unique(lossInfoCountyOnly$FIPS))) - 1L
  # range(diff(lossInfoCountyOnly$event)) # Is sorted
  eventLossesCountyOnly = lossInfoCountyOnly[, c('event', 'countyInd', 'loss')][
    order(event, countyInd), ][
      , list(list(list(zbCountyInd = countyInd, loss = loss))), by = event]
  eventLossesCountyOnly = eventLossesCountyOnly$V1
  save(eventLossesCountyOnly, file = 'data/eventLossesCountyOnly.RData')
  
  
  countyFIPS = countyInfo$FIPS[countyInfo$FIPS %% 1000L != 0L]
  eventLossesCountyOnly = lapply(eventLosses, function(x)
  {
    u = match( countyInfo$FIPS[x[[1]] + 1L], countyFIPS )
    ind = !is.na(u)
    list(countyIndZb = u[ind] - 1L, loss = x[[2]][ind])
  })
  
  
  source('code/R/rfuns.R')
  source("code/CharlieLib/R/CharlieSourceCpp.R")
  clrs = unique(colorRampPalette(c(
    'white', 'yellow', 'orange', 'red', 'darkred', 'black'))(10000))
  
  
  # find min distance between two counties' centers.
  source('code/R/rfuns.R')
  CharlieSourceCpp('code/cpp/cppFuns.cpp', optFlag = '-Ofast', sanitize = F)
  pxLen = getPixelLen(countyInfo[, c('lon', 'lat')], r = 0.999)
  subareaPxIndicesAndImageDim = getSubareaPxIndicesAndImageDim(
    countyInfo$lon, countyInfo$lat, pxLen)
  tmpImg = makeImage(eventLossesCountyOnly[[1]], subareaPxIndicesAndImageDim)
  kern = makeGauKernel(1001)
  tmpImgConv11 = imgConv(tmpImg, kern)
  tmpImgConv11viaSparse = testSparseConv(tmpImg, kern)
  range(tmpImgConv11viaSparse - tmpImgConv11)
  
  
  tmpClr = getColors(tmpImgConv11, clrs, method = "gau")
  tmpImgConv11ForPlot = getImgForColoring(tmpImgConv11)
  par(mar = c(0, 0, 0, 0))
  image( t(tmpImgConv11ForPlot) [, nrow(tmpImgConv11ForPlot):1], 
         col = tmpClr, xaxt = "n", yaxt = "n", bty  = "n")
  
  
  # makeImage(List eventLoss, List subareaPxIndicesAndImageDim)
  
  
  tmp = expandDim(eventLosses = eventLosses, countyInfo = countyInfo,
                  NNs = as.integer(c(8, 16, 32, 64)), takeAvg = TRUE)
  
  
  tmp = getCountyNN(countyInfo$FIPS, countyInfo$lon, countyInfo$lat, K = 64)
  
  
  
  
  
  
  # For each county, find its 8, 16, 32, 64 nearest neighbors.
  Cdist = dist(countyInfo[, c('lon', 'lat')])
  Cdist = as.matrix(Cdist)
  countyNNindZb = apply(Cdist, 2, function(x)
  {
    order(x)[2:257] - 1L
  })
  
  
  
  
  LD_PRELOAD="/lib64/libasan.so.5 /lib64/libubsan.so.1"  R
  
  
  load("data/yearEvent.RData")
  load("data/source/stateBorderLine.RData")
  load("data/regionIDs.RData")
  load("data/eventTlosses.RData")
  load("data/trueEPs.RData")
  
  
  load("data/countyInfo.RData")
  load("data/eventLosses.RData")
  source('code/R/rfuns.R')
  source("code/CharlieLib/R/CharlieSourceCpp.R")
  CharlieSourceCpp('code/cpp/expandDim.cpp', sanitize = F,  optFlag = '-Ofast')
  newDimStart = max(unlist(lapply(eventLosses, function(x) max(x[[1]])))) + 1L
  
  
  tmp = expandDim(
    eventLosses = eventLosses[1:10], countyInfo = countyInfo, 
    newDimStart = newDimStart, NNs = as.integer(c(8, 16, 32, 64)), 
    takeAvg = TRUE, useCppAccelerator = F)
  
  
  NNexpandedEventLosses = expandDim(
    eventLosses = eventLosses, countyInfo = countyInfo, 
    newDimStart = newDimStart, NNs = as.integer(c(8, 16, 32, 64)), 
    takeAvg = TRUE, useCppAccelerator = TRUE)
  save(NNexpandedEventLosses, file = "data/NNexpandedEventLosses-8-16-32-64.RData")
  
  
  CharlieSourceCpp('code/cpp/knns.cpp', sanitize = F, optFlag = '-Ofast')
  neis = knns(eventLosses = NNexpandedEventLosses, candidateMat = cands, keep = 100, 
              distance = "Euclidean", LPp = 0.5,
              maxCore = 1000, verbose = T, 
              regionW = numeric(0)
              # , regionW = regionW
  )
  emap = data.table::data.table(
    original = eventTlosses$event, mapped =  c(
      eventTlosses$event[1:(nrow(eventTlosses) - ncol(neis))],
      eventTlosses$event[neis[1,] + 1L]))
  save(emap, file = "data/EuclideanExpanded-8-16-32-64.RData")
  newEPs = generateEPs(emap, regionIDs, eventLosses, yearEvent)
  countryInd = 1L
  stateInd = which(regionIDs %% 1000L == 0L)[-1]
  countyInd = which(regionIDs %% 1000L != 0L) 
  cod(newEPs = newEPs$aggLoss, trueEPs = trueEPs$aggLoss)
  cod(newEPs = newEPs$aggLoss[countryInd], trueEPs = trueEPs$aggLoss[countryInd])
  cod(newEPs = newEPs$aggLoss[stateInd], trueEPs = trueEPs$aggLoss[stateInd])
  cod(newEPs = newEPs$aggLoss[countyInd], trueEPs = trueEPs$aggLoss[countyInd])
  
  
  
  
  
  
  
  
  
  
  
  # Test things without country and state dimensions.
  LD_PRELOAD="/lib64/libasan.so.5 /lib64/libubsan.so.1"  R
  
  
  load("data/yearEvent.RData")
  load("data/source/stateBorderLine.RData")
  load("data/regionIDs.RData")
  load("data/eventTlosses.RData")
  load("data/trueEPs.RData")
  
  
  load("data/countyInfo.RData")
  load("data/eventLosses.RData")
  source('code/R/rfuns.R')
  source("code/CharlieLib/R/CharlieSourceCpp.R")
  CharlieSourceCpp('code/cpp/expandDim.cpp', sanitize = F,  optFlag = '-Ofast')
  newDimStart = max(unlist(lapply(eventLosses, function(x) max(x[[1]])))) + 1L
  
  
  # tmp = expandDim(
  #   eventLosses = eventLosses[1:10], countyInfo = countyInfo, 
  #   newDimStart = newDimStart, NNs = as.integer(c(8, 16, 32, 64)), 
  #   takeAvg = TRUE, useCppAccelerator = F)
  
  load("data/NNexpandedEventLosses-8-16-32-64.RData")
  dimStart = max(unlist(lapply(eventLosses, function(x) max(x[[1]])))) + 1L
  NNexpandedEventLossesNoStatesCountry = lapply(
    NNexpandedEventLosses, function(x)
  {
    b = x[[1]] >= dimStart | x[[1]] %% 1000L != 0L
    list(x[[1]][b], x[[2]][b])
  })
  
  
  # NNexpandedEventLosses = expandDim(
  #   eventLosses = eventLosses, countyInfo = countyInfo, 
  #   newDimStart = newDimStart, NNs = as.integer(c(8, 16, 32, 64)), 
  #   takeAvg = TRUE, useCppAccelerator = TRUE)
  # save(NNexpandedEventLosses, file = "data/NNexpandedEventLosses-8-16-32-64.RData")
  # load("data/NNexpandedEventLosses-8-16-32-64.RData")
  
  
  CharlieSourceCpp('code/cpp/knns.cpp', sanitize = F, optFlag = '-Ofast')
  neis = knns(eventLosses = NNexpandedEventLossesNoStatesCountry, 
              candidateMat = cands, keep = 100, 
              distance = "L1", LPp = 0.5,
              maxCore = 1000, verbose = T, 
              regionW = numeric(0)
              # , regionW = regionW
  )
  emap = data.table::data.table(
    original = eventTlosses$event, mapped =  c(
      eventTlosses$event[1:(nrow(eventTlosses) - ncol(neis))],
      eventTlosses$event[neis[1,] + 1L]))
  save(emap, file = "data/L1Expanded-8-16-32-64-noCountryStates.RData")
  newEPs = generateEPs(emap, regionIDs, eventLosses, yearEvent)
  countryInd = 1L
  stateInd = which(regionIDs %% 1000L == 0L)[-1]
  countyInd = which(regionIDs %% 1000L != 0L) 
  cod(newEPs = newEPs$aggLoss, trueEPs = trueEPs$aggLoss)
  cod(newEPs = newEPs$aggLoss[countryInd], trueEPs = trueEPs$aggLoss[countryInd])
  cod(newEPs = newEPs$aggLoss[stateInd], trueEPs = trueEPs$aggLoss[stateInd])
  cod(newEPs = newEPs$aggLoss[countyInd], trueEPs = trueEPs$aggLoss[countyInd])
  
  
  NNexpandedEventLossesNoStatesCountryBinary = lapply(
    NNexpandedEventLossesNoStatesCountry, function(x)
    {
      b = x[[2]] > 1e-10
      list(x[[1]], as.numeric(b))
    })
  CharlieSourceCpp('code/cpp/knns.cpp', sanitize = F, optFlag = '-Ofast')
  # For L1 and Euclidean under binary setting, the results will be the same.
  neis = knns(eventLosses = NNexpandedEventLossesNoStatesCountry, 
              candidateMat = cands, keep = 100, 
              distance = "binary", LPp = 0.5,
              maxCore = 1000, verbose = T, 
              regionW = numeric(0)
              # , regionW = regionW
  )
  emap = data.table::data.table(
    original = eventTlosses$event, mapped =  c(
      eventTlosses$event[1:(nrow(eventTlosses) - ncol(neis))],
      eventTlosses$event[neis[1,] + 1L]))
  save(emap, file = "data/L1Expanded-8-16-32-64-noCountryStates-binary.RData")
  CharlieSourceCpp('code/cpp/cppFuns.cpp')
  newEPs = generateEPs(emap, regionIDs, eventLosses, yearEvent)
  countryInd = 1L
  stateInd = which(regionIDs %% 1000L == 0L)[-1]
  countyInd = which(regionIDs %% 1000L != 0L) 
  cod(newEPs = newEPs$aggLoss, trueEPs = trueEPs$aggLoss)
  cod(newEPs = newEPs$aggLoss[countryInd], trueEPs = trueEPs$aggLoss[countryInd])
  cod(newEPs = newEPs$aggLoss[stateInd], trueEPs = trueEPs$aggLoss[stateInd])
  cod(newEPs = newEPs$aggLoss[countyInd], trueEPs = trueEPs$aggLoss[countyInd])
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  err = unlist(NNexpandedEventLosses) / unlist(tmp) - 1
  range(err[!is.nan(err)])
  
  
  
  
  
  k = sample(length(eventLosses), 10)
  tmp = expandDim(eventLosses = eventLosses[k],
                  countyInfo = countyInfo, newDimStart = newDimStart, 
                  NNs = c(8, 16, 32, 64))
  dimnames(tmp[[1]]) = NULL
  for (i in 1:length(tmp[[2]])) names(tmp[[2]][[i]][[2]]) = NULL
  nmat = tmp[[1]]; tmp = tmp[[2]]
  
  
  NNexpandedEventLosses = makeExpandedEvents(
    eventLosses = eventLosses, nmat = nmat, dimStart = newDimStart, 
    NNs = as.integer(c(8, 16, 32, 64)), takeAvg = T, maxCore = 10000) 
  
  
  err = unlist(tmp) / unlist(NNexpandedEventLosses) - 1
  err = (err[!is.nan(err)])
  range(err)
  
  
  
  
  
  
  for (i in 1:length(tmp))
  {
    l = -(1:length(eventLosses[[i]][[1]]))
    tmp[[i]][[1]] = tmp[[i]][[1]][l]
    tmp[[i]][[2]] = tmp[[i]][[2]][l]
    tmp2[[i]][[1]] = tmp2[[i]][[1]][l]
    tmp2[[i]][[2]] = tmp2[[i]][[2]][l]
  }
  tmp = tmp[[1]]; tmp2 = tmp2[[1]]
  str(tmp); str(tmp2)
  
  
  # Manual checking
  if (F)
  {
    dimStart = max(unlist(lapply(eventLosses, function(x) max(x[[1]])))) + 1L
    e = eventLosses[[1]]
    tmpLossesForAllCounties = apply(nmat, 2, function(x)
    {
      sum(e$loss[e$zbRegionInd %in% x]) / nrow(nmat)
    })
    tmpind = which(tmpLossesForAllCounties > 1e-10)
    tmpind = tmpind + dimStart - 1L
    
    
  }
  
  
  
  
  tmp = dist(countyInfo[, c('lon', 'lat'), ])
  tmpMin = min(tmp)
  pixelEdgeLen = tmpMin / 2 * 0.999
  xlim = range(countyInfo$lon)
  ylim = range(countyInfo$lat)
  npixelX = as.integer(round(diff(xlim) / pixelEdgeLen) + 1)
  npixelY = as.integer(round(diff(ylim) / pixelEdgeLen) + 1)
  
  
  
  
  
  
  
  
  
  
}























