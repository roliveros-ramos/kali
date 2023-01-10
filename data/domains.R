domains = NULL
# domains$name = list(x=c(lonWest, lonEast), y=c(latSouth, latNorth))
domains$world  = list(x=c(-180,+180),  y=c(-90, +90))
domains$world2  = list(x=c(0, 360),  y=c(-90, +90))
# America
domains$peru      = list(x=c(-90,-70),   y=c(-20,0)) 
domains$peru2     = list(x=c(-88,-70),   y=c(-20,-2))
domains$peru3     = list(x=c(-93,-70),   y=c(-20,6))
domains$peruS     = list(x=c(-75,-70),   y=c(-20,-15))
domains$peruN     = list(x=c(-86,-78),   y=c(-10, -2))
domains$peruC     = list(x=c(-80,-74),   y=c(-16,-10))
domains$peruNC    = list(x=c(-87,-73),   y=c(-16,-2))
domains$peps      = list(x=c(-100,-70),  y=c(-40,10))
domains$ESPacific = list(x=c(-100,-65),  y=c(-50,10))
domains$pacific   = list(x=c(+110, -70), y=c(-55,70))
domains$EPO       = list(x=c(-150,-70),  y=c(-50,50))
domains$CentralA  = list(x=c(-95,-75),  y=c(6, 16))
domains$central_america  = list(x=c(-95,-75),  y=c(6, 16))
# Africa
domains$southern_africa  = list(x=c(4, 40), y=c(-45, -17))
domains$southern_benguela  = list(x=c(14, 28), y=c(-38, -28))
# Europe
domains$mediterranean_sea  = list(x=c(-7, 37), y=c(30, 46))
domains$medsea = domains$mediterranean_sea

# OSMOSE domains 
domains$ben = domains$southern_benguela
domains$gol = list(x=c(2.5,8.5), ylim=c(40.5,44.5))
