
require( "RMTstat" )

source( "rank.est.R" )
source( "noise.est.R" )


noise.est.eigs.kn <- function( ell, k, n, p, maxiter=30, tol=1e-5 ) {
    if( length( ell ) < p )
        ell <- c( ell, rep( 0, p-length( ell ) ) )
    
    sigma2.0 <- noise.est.eigs.fbk( ell, k, n, p )
    
    for( iter in seq_len( maxiter ) ) {
        
        # solve quadratic equation for rho, given sigma and eigenvalues
        tmp <- ell[ 1:k ] + sigma2.0 * ( ( n-p+k )/n )
        if( min( tmp^2 - 4 * ell[ 1:k ] * sigma2.0 ) < 0  ) {
            break # otherwise get complex valued solutions
        }
        rho.est <- ( tmp + sqrt( tmp^2 - 4 * ell[ 1:k ] * sigma2.0 ) ) / 2

        if( min( ell[ 1:k ] - rho.est ) < 0 )
            warning( "Error getting consistent noise estimante." )
     
        delta.l.rho <- pmax( ell[ 1:k ] - rho.est, 0 )
        sigma2.new  <- ( sum( ell[ (k+1):p ] ) + sum( delta.l.rho ) )/( p-k )

        if( abs( sigma2.new - sigma2.0 )/sigma2.0 < tol ) {
            break
        } else {
            sigma2.0 <- sigma2.new
        }

        if( iter == maxiter )
            warning( "Did not converge after ", maxiter, " iterations." )
    }
    
        
    sigma2 <- sigma2.0
    sigma2
}
noise.est.kn <- noise.est( noise.est.eigs.kn )

rank.est.eigs.kn <- function( ell, n, p, maxrank, beta=1, alpha=0.001, 
                              maxiter=30, tol=1e-5 ) {
    if( length( ell ) < p )
        ell <- c( ell, rep( 0, p-length( ell ) ) )
    
    s.Wishart <- qtw( alpha, beta, lower.tail=FALSE )     

    sigma2.est.arr <- rep( NA, maxrank )

    rank <- 0
    for( k in seq_len( maxrank ) ) {
        par      <- WishartMaxPar( n, p-k, beta=beta )
        mu.np    <- par$centering
        sigma.np <- par$scaling
        
        sigma2.est.k <- suppressWarnings(
                            noise.est.eigs.kn( ell, k, n, p, maxiter, tol ) )
        sigma2.est.arr[ k ] <- sigma2.est.k
        at.least.k.signals  <-  ( ell[ k ] 
                                  > 
                                  sigma2.est.k
                                  * ( mu.np + s.Wishart * sigma.np ) )
                                 
        if( at.least.k.signals ) {
            rank <- k
        } else {
            break
        }
    }

    rank
}
rank.est.kn <- rank.est( rank.est.eigs.kn )
