
source( "kritchman-nadler.R" )


minmax.cutoff <- function( n, p, sigma2, lambda.min=sqrt(p/n)+n^(-1/3),
                           cost.include=1, cost.exclude=1, beta=1 ) {
    tw    <- WishartMaxPar( n, p, sigma2, beta )
    spike <- WishartSpikePar( lambda.min, n, p, sigma2, beta )
    
    T <- uniroot( interval=c( -10, 6 ),
             function( t ) { 
                 ( ( cost.include
                     * (1 - ptw( ( t - tw$centering )/tw$scaling, beta ) ) )
                   -
                   ( cost.exclude
                     * pnorm( ( t - spike$centering )/spike$scaling ) ) )
             })$root
    T
}

        T[ i ] <- minmax.cutoff( n, p, sigma2, lambda.min, 
                                 cost.include, cost.exclude*(n-i+1) )
    }
    
    T
}


select.minmax <- function( n, N, r0, sigma2, l.signal, l.noise ) {
    cutoffs <- minmax.cutoffs( n, N, sigma2 )

    if( r0 > 0 && l.signal < cutoffs[ r0 ] ) {
        r <- r0 - 1
    } else if( r0 < n-1 && l.noise > cutoffs[ r0+1 ] ) {
        r <- r0 + 1
    } else {
        r <- r0
    }
    
    r    
}

rank.est.eigs.pw <- function( ell, n, p, maxrank, beta=1, ... ) {
    if( length( ell ) < p )
        ell <- c( ell, rep( 0, p-length( ell ) ) )

    rank <- 0
    for( k in seq_len( maxrank ) ) {
        sigma2.est <- suppressWarnings( 
                          noise.est.eigs.kn( ell, k, n, p, ... ) )
        cost.inc   <- 1
        cost.exc   <- p - k + 1
        cutoff     <- minmax.cutoff( n, p, sigma2.est, 
                                     cost.include=cost.inc,
                                     cost.exclude=cost.exc,
                                             beta=beta )
        
        if( ell[ k ] > cutoff ) {
            rank <- k
        } else {
            break
        }
    }

    rank    
}
rank.est.pw <- rank.est( rank.est.eigs.pw )
