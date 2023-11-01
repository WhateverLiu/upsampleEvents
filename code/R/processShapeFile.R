





# Download shape file from this site:
# https://www.census.gov/geographies/mapping-files/time-series/geo/carto-boundary-file.html
# file.copy(from = '/finance_develop/Charlie/backup/Desktop20200130/countyCor/usStateBorderCoor.Rdata', to = 'data')
# emstreeR::ComputeMST(usStateBorderCoor, verbose = F)
# Plotting shape files are too slow. Seek simplified method.
load('data/stateBdLine.RData')
# Plot the states borders and check
if (T)
{
  xlim = range(lapply(bdLine, function(x) range(x[,1])))
  ylim = range(lapply(bdLine, function(x) range(x[,2])))
  plot(0, col = "white", xlim = xlim, ylim = ylim, bty = "n")
  tmp = lapply(bdLine, function(x) lines(x))
}


countyInfo = rvest::html_table(rvest::read_html(
  'https://en.wikipedia.org/wiki/User:Michael_J/County_table'))[[1]]
countyInfo = as.data.frame(countyInfo)
countyInfo = countyInfo[c('FIPS', 'State', 'County [2]', 'Longitude', 'Latitude')]
colnames(countyInfo) = c('FIPS', 'state', 'county', 'lon', 'lat')
countyInfo$lon = gsub('°', '', countyInfo$lon)
countyInfo$lon = gsub('–', '-', countyInfo$lon)
countyInfo$lat = gsub('°', '', countyInfo$lat) 
countyInfo$lon = as.numeric(countyInfo$lon)
countyInfo$lat = as.numeric(countyInfo$lat)




save(countyInfo, file = 'data/countyInfo.RData')


dat = data.table::setDF(data.table::fread(
  'data/TouchstoneRe_Stochastic100kWSST_Catalog_M021_20230503.l__'))
colnames(dat) = c('modelNumber', 'year', 'event', 'day', 'country-FIPS', 
                  'unknown',
                  'lineOfBusiness1loss', 'lineOfBusiness2loss', 
                  'lineOfBusiness3loss', 'lineOfBusiness4loss')
lossInfo = dat[c('year', 'event')]
lossInfo$FIPS = as.integer(substring(dat$`country-FIPS`, first = 4))
lossInfo$loss = dat$lineOfBusiness1loss + dat$lineOfBusiness2loss + 
  dat$lineOfBusiness3loss + dat$lineOfBusiness4loss
lossInfo$dimIDzb = match(lossInfo$FIPS, sort(unique(lossInfo$FIPS))) - 1L
lossInfo = lossInfo[order(lossInfo$event, lossInfo$dimIDzb), , drop = F]
yearEvent = unique(lossInfo$event)
yearEvent = data.frame(year = lossInfo$year[match(yearEvent, lossInfo$event)], 
                       event = yearEvent)
tmp = aggregate(list(ind = 1:length(lossInfo$event)), 
                list(event = lossInfo$event), function(x) x)
eventLosses = lapply(tmp$ind, function(x)
{
  list(ind = lossInfo$dimIDzb[x], loss = lossInfo$loss[x])
})
tmp = aggregate(list(event = lossInfo$event), list(dimIDzb = lossInfo$dimIDzb), 
                function(x) x)
regionEvents = tmp; rm(tmp)
saveRDS(regionEvents, file = "data/regionEvents.rds")
saveRDS(lossInfo, file = "data/lossInfo.rds")
saveRDS(eventLosses, file = "data/eventLosses.rds")
saveRDS(yearEvent, file = "data/yearEvent.rds")





# lossInfo$dimIDzb = sort(unique(lossInfo$FIPS))
  

# raster::plot(UScounties, xlim = c(-125, 49), ylim = c(25.07, 49.35))
raster::plot(UScounties, xlim = c(-125, -67), ylim = c(25.07, 49.35), bty = "L")
raster::plot(USstates, xlim = c(-125, -67), ylim = c(25.07, 49.35), bty = "L")





aggInd = aggregate(list(ind = 1:length(UScounties$STATEFP)),
                   by = list(state = UScounties$STATEFP), function(x) x, 
                   simplify = F)
aggInd = as.list(aggInd)
aggInd$countyName = lapply(aggInd$ind, function(x) UScounties$NAME[x])
aggInd$geoID = lapply(aggInd$ind, function(x) UScounties$GEOID[x])
aggInd$polygon = lapply(aggInd$ind, function(x) 
{ 
  lapply(UScounties@polygons[x], function(u) 
  {
    rst = u@Polygons[[1]]@coords; dimnames(rst) = NULL; rst
  })
})
agg = aggInd; rm(aggInd)
agg$ind = NULL
states = new.env()






tmp = tigris::counties('Texas', cb = T)
plot(tmp$geometry[[1]][[1]][[1]], type = 'l')


library(tigris)
library(ggplot2)

me <- counties("Maine", cb = TRUE)
plot(me$geometry[[1]][[1]][[1]])

gg <- ggplot()
gg <- gg + geom_sf(data = me, color="black",
                   fill="white", size=0.25)
gg

## End(Not run)







