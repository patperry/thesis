
rank.est.wk <- function( ell, n ) {
    p   <- length( ell )
    np  <- min( n, p )
    
    a   <- rep( NA, np )
    mdl <- rep( NA, np )

    for( k in seq( from=0, to=( np-1 ), by=1 ) ) {
        a[ k+1 ]   <- mean( ell[ (k+1):p ] )
        mdl[ k+1 ] <- -( p-k ) * n * ( mean( log( ell[ (k+1):p ] ) ) 
                                       - log( a[ k+1 ] ) )
        mdl[ k+1 ] <- mdl[ k+1 ] + 1/2 * k * ( 2*p - k ) * log( n )                             
    }
    
    rank <- which.min( mdl ) - 1
    rank
}    
