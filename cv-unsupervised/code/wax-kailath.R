
source( "rank.est.R" )

# MDL estimator for number of signals proposed by Wax & Kailath in 
# "Detection of signals by information theoretic criteria", IEEE
# Transactions on Acoustics, Speech and Signal Processing 33 (2) (1985)
# 387--392.
rank.est.eigs.mdl <- function( ell, n, p, maxrank ) {
    if( length( ell ) < p )
        ell <- c( ell, rep( 0, p-length( ell ) ) )
    
    mdl <- rep( NA, maxrank+1 )

    for( k in seq( from=0, to=maxrank, by=1 ) ) {
        mdl[ k+1 ] <- ( -( p-k ) * n * ( mean( log( ell[ (k+1):p ] ) ) 
                                         - log( mean( ell[ (k+1):p ] ) ) )
                        + 1/2 * k * ( 2*p - k ) * log( n ) )                           
    }
    
    rank <- which.min( mdl ) - 1
    rank
}    
rank.est.mdl <- rank.est( rank.est.eigs.mdl )
