
require( "plyr" )
source( "../../common/code/spiked-data.R" )


loss.sim <- function( kmax, spike, n, p, var=1, seed=NULL, ... ) {
    if( !is.null( seed ) )
        set.seed( seed, "Mersenne-Twister" )
    
    k     <- seq_len( kmax+1 ) - 1
    spike <- sort( spike, decreasing=TRUE )
    sim   <- spiked.data( spike, n, p, var=var, ... )
                            
    resid <- sim$signal
    spec2 <- rep( NA, kmax+1 )
    frob2 <- rep( NA, kmax+1 )
    
    if( kmax > 0 ) {
        spec2[1] <- svd( resid, nu=0, nv=0 )$d[1]^2
        frob2[1] <- mean( resid^2 )
    }
    
    for( i in seq_len( min( kmax, n, p ) ) ) {
        u <- sim$u.est[,i,drop=FALSE]
        v <- sim$v.est[,i,drop=FALSE]
        d <- sim$d.est[i]
        
        resid <- resid - (d * u) %*% t(v)
        
        spec2[i+1] <- svd( resid, nu=0, nv=0 )$d[1]^2
        frob2[i+1] <- mean( resid^2 )
    }
    
    spec2 <- spec2 / n
    frob2 <- frob2 * p
    
    data.frame( k=k, spec2=spec2, frob2=frob2 )
}