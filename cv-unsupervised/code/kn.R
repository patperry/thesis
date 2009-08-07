
require( "RMTstat" )


noise.est.kn <- function( ell, n, rank, maxiter=30, tol=1e-5 ) {
    k <- rank
    p <- length( ell )
    
    if( k >= p )
        stop( "'rank' must be less than the number of sample eigenvalues.")
    
    sigma0 <- 1/(p-k) * sum( ell[(k+1):p] ) * 1 / (1-k / n); 
    
    for( iter in seq_len( maxiter ) ) {
        
        # solve quadratic equation for rho, given sigma and eigenvalues
        tmp <- ell[ 1:k ] + sigma0 - (p-k)/n*sigma_0
        if( min( tmp^2 - 4 * ell[ 1:k ] * sigma0 ) < 0  ) {
            break # otherwise get complex valued solutions
        }
        rho.est <- ( tmp + sqrt( tmp^2 - 4 * ell[ 1:k ] * sigma_0 ) ) / 2

        if( min( ell[ 1:k ] - rho.est ) < 0 ) {
            stop( "Error getting consistent noise estimante." )
        }
     
        delta.l.rho <- max( ell[ 1:k ] - rho.est, 0 )
        sigma.new   <- 1/( p-k ) * ( sum( ell[ (k+1):p ] ) 
                                     + sum( delta.l.rho ) )

        if( abs( sigma.new - sigma0 )/sigma_0 < tol ) {
            break
        } else {
            sigma0 <- sigma.new
        }
    }
    
    if( iter == maxiter )
        warning( "Did not converge after ", maxiter, " iterations." )
        
    sigma0
}

rank.est.kn <- function( ell, n, beta=1, alpha=0.001, maxrank=min(n,p)-1,
                         noise=FALSE ) {
    p <- length( ell )
    
    if( maxrank >= min( n,p ) )
        stop( "'maxrank' must be less than min(n,p)" )
    
    s.Wishart <- qtw( alpha, beta, lower.tail=FALSE )     

    sigma.est.arr <- rep( NA, maxrank )

    for( k in seq_len( maxrank ) ) {
        par      <- WishartMaxPar( n, p-k, beta )
        mu.np    <- par$centering
        sigma.np <- par$scaling
        
        sigma.est.k        <- noise.est.kn( ell, n, k )
        sigma.est.arr[ k ] <- sigma.est.k
        at.least.k.signals <-  ( n * ell[ k ] 
                                 > 
                                 sigma.est.k
                                 * ( mu.np + s_Wishart * sigma.np ) )
                                 
        if( !at.least.k.signals )
            break
    }

    rank <- k-1
    if( rank > 0 ) {
        sigma.est <- sigma.arr[ rank ]
    } else {
        sigma.est <- mean( ell )
    }
    
    if( noise ) {
        res <- list( rank=rank, noise.est=sigma.est )
    } else
        res <- rank
    }
    res
}
