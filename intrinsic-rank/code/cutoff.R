
require( "ggplot2" )
require( "plyr" )


loss.F.factor <- function( spike, concentration, var=1 ) { 
	mu     <- spike
	gamma  <- concentration
	sigma2 <- var
	
	loss <- ifelse( mu > 1/sqrt( gamma ),
				( sigma2/( gamma * mu^2 ) )*( 3*sigma2 + ( gamma+1 )*mu ),
				1 + ( sigma2/mu )*( 1 + 1/sqrt( gamma ) )^2 )
	loss
}

loss.F.cutoff <- function( concentration, var=1 ) {
	gamma  <- concentration
	sigma2 <- var
	
	a      <- 0.5*( 1 + 1/gamma )
	cutoff <- sigma2*( a + sqrt( a^2 + 3/gamma ) )
	cutoff
}

gamma <- c( 1/16, 1/4, 1.0, 4, 16 )
loss  <- adply( gamma, 1, 
	 		function( gamma ) {
	 			mu    <- seq( from=0, to=25, length=200)
	 			alpha <- loss.F.factor( mu, gamma )
	  			data.frame( gamma=gamma, mu=mu, alpha=alpha )
			}
		)[,-1]

p1 <- ( ggplot( loss, aes( mu, alpha, colour=gamma ) )
        + labs( x="Signal Strength", y="Frobenius Loss Penalty", colour="Aspect Ratio" ) 
        + theme_bw()
        + geom_line( aes( group=gamma ) )
        + scale_colour_gradient( breaks=sort( gamma, decreasing=TRUE ) )
        + scale_y_log2()
      )
p1
pdf( file="../plots/frobenius-loss-penalty.pdf", width=6, height=4 )
print( p1 )
dev.off()



cutoff <- {
    n     <- 100
    gamma <- exp( seq( from=log( 0.04 ), to=log( 25 ), length=n ) )
    mu    <- c( loss.F.cutoff( gamma ), 1/sqrt( gamma ) )
    type  <- rep( c( "include", "detect" ), each=n )
    data.frame( gamma=rep( gamma, 2 ), mu=mu, type=type )
}
p2 <- ( ggplot( cutoff, aes( gamma, mu, colour=type ) )
        + labs( x="Aspect Ratio", y="Signal Strength", colour="Threshold" )
        + scale_colour_hue( breaks=c("include", "detect"), labels=c("Inclusion", "Detection" ) )
        + geom_line()
        + scale_x_log2()
        + theme_bw() )
p2
pdf( file="../plots/frobenius-loss-cutoff.pdf", width=6, height=4 )
print( p2 )
dev.off()

