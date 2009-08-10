
require( "ggplot2" )
require( "plyr" )
source( "spikesim.R" )

make_plot <- function( n=100, p=50, maxrank=20, 
                       spike=5:10, factors="gaussian", 
                       main="Weak Gaussian Factors", 
                       factor.seed=1, noise.seed=10001 ) {
    sim.white <- spikesim( spike, n, p, maxrank, factors, factors, "white",
                           factor.seed=factor.seed, noise.seed=noise.seed )
    sim.colored <- spikesim( spike, n, p, maxrank, factors, factors, "colored",
                             factor.seed=factor.seed, noise.seed=noise.seed )
    sim.heavy <- spikesim( spike, n, p, maxrank, factors, factors, "heavy",
                           factor.seed=factor.seed, noise.seed=noise.seed )

    cv <- rbind( data.frame( summary( sim.white )$cv, noise="white" ),
                 data.frame( summary( sim.colored )$cv, noise="colored" ),
                 data.frame( summary( sim.heavy )$cv, noise="heavy" ) )

    truth <- rbind( data.frame( summary( sim.white )$truth, noise="white" ),
                    data.frame( summary( sim.colored )$truth, noise="colored" ),
                    data.frame( summary( sim.heavy )$truth, noise="heavy" ) )

    plt <- ( ggplot( cv, aes( rank, msep, colour=rotation ) ) 
           + theme_bw() 
           + xlab( "Rank" )
           + ylab( "Prediction Error" )
           + scale_colour_hue( "Rotation", formatter=firstUpper )
           + opts( legend.position="none" )
           + opts( title=main )
           + facet_grid( noise ~ type,
                         labeller=function(variable, value) {
                                      if( variable == "noise" )
                                          paste( firstUpper( value ), "Noise" )
                                      else
                                          paste( firstUpper( value ), "CV" ) } )
           + geom_line( aes( rank, msem ), data=truth, 
                        colour=I("black"), linetype=I("dashed") )
           + geom_point( aes( rank, msem ), data=truth, 
                         colour=I("black"), size=I(1.5) )
           + geom_line() 
           + geom_errorbar( aes( x=rank, ymin=msep-se, ymax=msep+se ) ) 
           + geom_point( size=I(1.5) ) 
         )
    plt
}

n            <- 100
p            <- 50
maxrank      <- 20
spike.weak   <- 5:10
spike.strong <- 10*spike.weak

plt.wg <- make_plot( n, p, maxrank, spike.weak, "gaussian", 
                     main="Weak Gaussian Factors" )

plt.sg <- make_plot( n, p, maxrank, spike.strong, "gaussian", 
                     main="Strong Gaussian Factors" )

plt.ws <- make_plot( n, p, maxrank, spike.weak, "sparse", 
                     main="Weak Sparse Factors" )

plt.ss <- make_plot( n, p, maxrank, spike.strong, "sparse", 
                     main="Strong Sparse Factors" )
 
width  <- 4
height <- 6

pdf( "../plots/cvsvd-weak-gauss.pdf", width=width, height=height )
print( plt.wg )
dev.off()

pdf( "../plots/cvsvd-weak-sparse.pdf", width=width, height=height )
print( plt.ws )
dev.off()

pdf( "../plots/cvsvd-strong-gauss.pdf", width=width, height=height )
print( plt.sg )
dev.off()

pdf( "../plots/cvsvd-strong-sparse.pdf", width=width, height=height )
print( plt.ss )
dev.off()
