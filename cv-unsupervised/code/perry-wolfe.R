
source( "kritchman-nadler.R" )


minmax.cutoff <- function( n, p, sigma2, 
                           lambda.min=sigma2*( sqrt(p/n)+n^(-1/3) ),
                           cost.include=1, cost.exclude=1, beta=1 ) {
    tw    <- WishartMaxPar( n, p, sigma2, beta )
    spike <- WishartSpikePar( lambda.min, n, p, sigma2, beta )
    
    lower <- -10 * tw$scaling + tw$centering
    upper <-   6 * tw$scaling + tw$centering
    
    T <- uniroot( interval=c(lower, upper),
             function( t ) { 
                 ( ( cost.include
                     * (1 - ptw( ( t - tw$centering )/tw$scaling, beta ) ) )
                   -
                   ( cost.exclude
                     * pnorm( ( t - spike$centering )/spike$scaling ) ) )
             })$root
    T
}

rank.est.eigs.pw <- function( noise.est.eigs ) {
    function( ell, n, p, maxrank, beta=1, ... ) {
        if( length( ell ) < p )
            ell <- c( ell, rep( 0, p-length( ell ) ) )

        rank <- 0
        for( k in seq_len( maxrank ) ) {
            sigma2.est <- suppressWarnings( 
                              noise.est.eigs( ell, k, n, p, ... ) )
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
}

rank.est.pw.naive <- rank.est( rank.est.eigs.pw( noise.est.eigs.naive ) )
rank.est.pw.fbk   <- rank.est( rank.est.eigs.pw( noise.est.eigs.fbk ) )
rank.est.pw.kn    <- rank.est( rank.est.eigs.pw( noise.est.eigs.kn ) )
