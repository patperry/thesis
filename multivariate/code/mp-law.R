require( "RMTstat" )
require( "plyr" )
require( "ggplot2" )

gamma <- c( 0.25, 1.00, 4.00 )

mp <- ldply( gamma, 
          function( g, ... ) { 
              a <- ( 1 - sqrt( 1/g ) )^2; 
              b <- ( 1 + sqrt( 1/g ) )^2; 
              x <- seq(a, b, len=100 ); 
              y <- dmp( x, svr=g ); 
              data.frame( concentration=g, quantile=x, density=y ) } )
         
p <- ( ggplot( mp, aes( quantile, density, color=concentration ) )
       + labs( x="Quantile", y="Density", colour="Concentration" ) 
       + geom_line( aes( group=concentration ) )
       + scale_colour_gradient( breaks=sort( gamma, decreasing=TRUE ) )
       + theme_bw()
     )
p

pdf( file="../plots/mp-law.pdf", width=6, height=4 )
print( p )
dev.off()
