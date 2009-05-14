require( "RMTstat" )
require( "ggplot2" )
require( "plyr" )

tw <- {
    n    <- 200
    xmin <- -5
    xmax <-  5
    beta <- 1
    x    <- seq( xmin, xmax, len=n )
    d    <- dtw( x, beta=beta )
    data.frame( beta=beta, quantile=x, density=d )
}

p <- ( ggplot( tw, aes( quantile, density ) )
       + geom_line( aes( group=factor( beta ) ) )
       + geom_line( colour="red" )
       + labs( x="Quantile", y="Density" )
       + theme_bw()
     )
p

pdf( file="../plots/tw-law.pdf", width=5, height=4 )
print( p )
dev.off()
