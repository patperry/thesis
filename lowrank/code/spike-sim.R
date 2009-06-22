
source( "../../common/code/spiked-data.R" )
require( "RMTstat" )

# This demonstrates that below the threshold, the sample principal components
# are not very correlated with the population principal components.

make.noise <- function( n, p, reps, ... ) {
    es <- array( NA, c( n, p, reps ) )
    for( r in 1:reps ) {
        es[ ,,r ] <- noise( n, p, ... )
    }
    es
}

spike.sim <- function( spike, es ) {
    n    <- dim( es )[ 1 ]
    p    <- dim( es )[ 2 ]
    reps <- dim( es )[ 3 ]
    
    u <- factors( n, 1, "basis" )
    v <- factors( p, 1, "basis" )
    
    spike.est <- array( NA, c( reps ) )
    udot2.est <- array( NA, c( n, reps ) )
    vdot2.est <- array( NA, c( p, reps ) )
    
    for( r in 1:reps ) {
        e <- es[ ,,r ]
        x <- u %*% ( sqrt( spike*n ) * t(v) ) + e
        x.svd <- svd( x )
        
        u.est <- x.svd$u
        v.est <- x.svd$v
        udot2.est[ ,r ] <- ( t(u) %*% u.est )^2
        vdot2.est[ ,r ] <- ( t(v) %*% v.est )^2
        spike.est[ r  ] <- x.svd$d[ 1 ]^2 / n
    }
    
    list( spike=spike, spike.est=spike.est,
          udot2.est=udot2.est, vdot2.est=vdot2.est )
}

n <- 100
p <- n/2
spike <- .5*sqrt( p/n )
reps <- 500
set.seed( 0, "Mersenne-Twister" )
es  <- make.noise( n, p, reps, "white" )
sim <- spike.sim( spike, es )

# centering <- ( spike + 1 )*( 1 + p/(n*spike) )
# scaling   <- sqrt( 2*( 2*spike + 1 + p/n )*( 1 - p/( n*spike^2 ) )/n )
par( mfrow=c(3,1) )
for( i in 1:3 ) hist( sim$vdot2.est[i,], xlim=c(0, 1) )
