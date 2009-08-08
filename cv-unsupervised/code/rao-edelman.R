source( "rank.est.R" )

rank.est.eigs.aic <- function( ell, n, p, maxrank, beta=1 ) {
    if( length( ell ) < p )
        ell <- c( ell, rep( 0, p-length( ell ) ) )
    
    t   <- rep( NA, maxrank+1 )
    aic <- rep( NA, maxrank+1 )
    
    for( k in seq( from=0, to=maxrank, by=1 ) ) {
        t[ k+1 ] <- ( p * ( ( ( p-k ) * sum( ell[ (k+1):p ]^2 )
                              / ( ( sum( ell[ (k+1):p ] ) )^2 ) )
                            - ( 1 + p/n ) )
                      -
                      ( 2/beta - 1 )*( p/n ) )
        aic[ k+1 ] <- ( beta/4 )*( n/p )^2 * t[ k+1 ]^2 + 2*( k+1 )
    }
    
    rank <- which.min( aic ) - 1
    rank
}    
rank.est.aic <- rank.est( rank.est.eigs.aic )
