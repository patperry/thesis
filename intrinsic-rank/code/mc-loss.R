
require( "ggplot2" )
require( "plyr" )
source( "loss-sim.R" )
source( "loss-theory.R" )

set.seed( 0, "Mersenne-Twister" )

conc  <- data.frame( conc=c( 4, 1, 1/4 ) )
size  <- data.frame( size=c( 144, 400, 1600, 4900 ) )
scale <- c(4, 2, 1, 1/2, 1/4)
kmax  <- 6

n.rep <- 500
seed  <- data.frame( seed=seq_len( n.rep ) )


threshold <- function( concentration=n/p, n=NA, p=NA, var=1 ) {
    gamma   <- concentration
    detect  <- 1/sqrt( gamma )
    include <- ((1 + 1/gamma)/2 
                + sqrt( ( (1 + 1/gamma)/2 )^2 + 3/gamma ))
    list( detect=var*detect, include=var*include )
}

n.conc <- nrow( conc )
n.size <- nrow( size )

theory <- ddply( conc, .(conc), splat( function( conc, ... ) {
              spike <- threshold( conc=conc )$include * scale
              loss  <- loss.theory( kmax, spike, conc=conc )
              ddply( loss, .(k), splat( function( frob2, spec2, ... ) 
                  data.frame( type=factor( c("frob2", "spec2") ),
                              loss=c(frob2, spec2) ) 
              ) )
          } ) )


params <- ddply( conc, .(conc), splat( function( conc, ... )
              ddply( size, .(size), splat( function( size, ... ) 
                  data.frame( n=round( sqrt( size*conc ) ),
                              p=round( sqrt( size/conc ) ) )
              ) )
          ) )

rawdata <- adply( params, 1, splat( function( n, p, ... ) {
               spike <- threshold( n=n, p=p )$include * scale
               
               ddply( seed, .(seed),  splat( function( seed, ... ) 
                   loss.sim( kmax, spike, n, p, seed=seed, 
                             left="uniform",
                            right="uniform",
                            noise="white" ) ) 
               )
           } ) )

sim <- ddply( rawdata, .(conc,size,n,p,k), splat( 
            function( n, p, k, frob2, spec2, ... ) {
                spike <- threshold( n=n[1], p=p[1] )$include * scale
                loss  <- loss.theory( k[1], spike, n[1], p[1] )[k[1]+1,]

                data.frame( type=factor( c("frob2", "spec2") ),
                          theory=c( loss$frob2, loss$spec2 ),
                            loss=c( mean( frob2 ), mean( spec2 ) ),
                         sd.loss=c( sd( frob2 ),     sd( spec2 ) ) )
            } ) )

combined <-
    rbind( cbind(theory, desc="theory", sd.loss=NA),
           ddply( sim, .(conc,k,type), splat( function( size, loss, sd.loss, ... )
                data.frame( desc=paste( "size", size, sep="" ),
                            loss=loss,
                         sd.loss=loss )
           ) ) )

#sim


plot_sim <- function( sim, type=c("frob2", "spec2") ) {
    t <- match.arg( type )

    sim <- subset( sim, type == t )
    # ggplot2 doesn't like us to have a column called "size"
    sim$np   <- sim$size
    sim$size <- NULL

    ylab <- switch(t, frob2="Squared Frobenius Loss",
                      spec2="Squared Spectral Loss" )

    p <- ( ggplot( sim, aes( k, loss, colour=I("red") ) ) 
           + theme_bw()
           + xlab( "Rank" )
           + ylab( ylab )
           + facet_grid( conc ~ np, scales="free_y" )
           + geom_line( aes( k, theory ), colour=I("black"), lty=2 )
           + geom_point()
           + geom_line()
           + geom_errorbar( aes( ymin=loss-sd.loss, 
                                 ymax=loss+sd.loss ), 
                            width=0.25 ) 
           + opts( legend.position="none" )
         )
    p
}

p1 <- plot_sim( sim, "frob2" )
p2 <- plot_sim( sim, "spec2" )

pdf( file="../plots/frob2-loss-sim.pdf", width=6, height=4 )
print( p1 )
dev.off()

pdf( file="../plots/spec2-loss-sim.pdf", width=6, height=4 )
print( p2 )
dev.off()

#p2 <- ( ggplot( subset( combined, type == "frob2" ), 
#                aes( k, loss, colour=desc ) ) 
#        + theme_bw()
#        + facet_grid( conc ~ ., scales="free_y" )
#        + geom_line()
#      )
#p2


    