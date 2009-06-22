
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



gamma  <- exp( seq( from=log( 0.04 ), to=log( 25 ), length=100 ) )
cutoff <- data.frame( gamma=gamma, cutoff=loss.F.cutoff( gamma ) )
p2 <- ( ggplot( cutoff, aes( gamma, cutoff ) )
        + labs( x="Aspect Ratio", y="Critical Signal Strength" )
        + geom_line( colour="red" )
        + scale_x_log2()
        + theme_bw() )
pdf( file="../plots/frobenius-loss-cutoff.pdf", width=4, height=4 )
print( p2 )
dev.off()

