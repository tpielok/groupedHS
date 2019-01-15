# Stackoverflow
map2color <- function(x, pal, limits = range(x)){
  pal[findInterval(x, seq(limits[1], limits[2], length.out = length(pal) + 1), 
                   all.inside=TRUE)]
}