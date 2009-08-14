
require( "bcv", "~/lib/R" )
require( "ggplot2" )
require( "plyr" )

motor <- read.csv( "../data/motor.csv" )

num_neuron    <- 49
num_condition <- 27
num_time      <- 2105/5 + 1

A <- matrix( NA, num_condition*num_time, num_neuron )

for( r in 1:nrow(motor) ) {
    neuron    <- motor$neuron[r]
    condition <- motor$condition[r]
    time      <- motor$time[r]
    response  <- motor$response[r]
    
    i      <- ( condition-1 )*num_time + time/5 + 1
    j      <- neuron
    A[i,j] <- response
}


T <- matrix( NA, num_condition, num_neuron*num_time )

for( r in 1:nrow(motor) ) {
    neuron    <- motor$neuron[r]
    condition <- motor$condition[r]
    time      <- motor$time[r]
    response  <- motor$response[r]
    
    i      <- condition
    j      <- ( neuron-1 )*num_time + time/5 + 1 
    T[i,j] <- response
}

set.seed( 0 )
cvw_5 <- cv.svd.wold( scale( A, center=TRUE, scale=FALSE ), 5, maxrank=48 )

set.seed( 0 )
cvg_2_2 <- cv.svd.gabriel( scale( A, center=TRUE, scale=FALSE ), 2, 2 )

set.seed( 0 )
cvg_3_3 <- cv.svd.gabriel( scale( A, center=TRUE, scale=FALSE ), 3, 3 )

set.seed( 0 )
cvg_4_4 <- cv.svd.gabriel( scale( A, center=TRUE, scale=FALSE ), 4, 4 )

set.seed( 0 )
cvg_5_5 <- cv.svd.gabriel( scale( A, center=TRUE, scale=FALSE ), 5, 5 )

set.seed( 0 )
cvg_2_49 <- cv.svd.gabriel( scale( A, center=TRUE, scale=FALSE ), 2, 49 )

cvw_5_sum    <- summary( cvw_5 )
cvg_2_2_sum  <- summary( cvg_2_2 )
cvg_2_49_sum <- summary( cvg_2_49 )

type <- c( rep( "wold-5",       cvw_5_sum$maxrank + 1 ),
           rep( "gabriel-2-2",  cvg_2_2_sum$maxrank + 1 ),
           rep( "gabriel-2-49", cvg_2_49_sum$maxrank + 1 ) )

rank <- c( 0:cvw_5_sum$maxrank,
           0:cvg_2_2_sum$maxrank,
           0:cvg_2_49_sum$maxrank )

pe <- c( cvw_5_sum$msep.mean,
         cvg_2_2_sum$msep.mean,
         cvg_2_49_sum$msep.mean )

pe.se <- c( cvw_5_sum$msep.se,
            cvg_2_2_sum$msep.se,
            cvg_2_49_sum$msep.se )

pe.se <- pe.se / pe[ 1 ]
pe <- pe / pe[ 1 ]


motor.cv <- data.frame( type=as.factor( type ), 
                        rank=rank, 
                          pe=pe, 
                       pe.se=pe.se )

plt <- ( ggplot( motor.cv, aes( rank, pe, colour=type ) )
         + theme_bw()
         + xlab( "Rank" )
         + ylab( "Estimated Prediction Error" )
         + scale_colour_hue( name="CV Method", formatter=function( x ) {
                g.2.2  <- match( "gabriel-2-2", x )
                g.2.49 <- match( "gabriel-2-49", x )
                w.5    <- match( "wold-5", x )
                y <- rep( NA, 3 )
                y[ g.2.2 ]  <- "Gabriel (2,2)-Fold"
                y[ g.2.49 ] <- "Gabriel (2,49)-Fold"
                y[ w.5 ]    <- "Wold 5-Fold"
                y
            })
         + geom_line()
         + geom_point()
         + geom_errorbar( aes( x=rank, ymin=pe-pe.se, ymax=pe+pe.se ) )
       )
plt

pdf( "../plots/motor-cv-est.pdf", 6, 4)
print( plt )
dev.off()
