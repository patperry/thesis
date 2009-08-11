
source( "rank.est.R" )

rank.est.eigs.fk <- function( ell, n, p, maxrank, alpha=0.05, beta=1 ) {
    
    ell     <- ell
    q       <- min( n, p )
    maxrank <- min( maxrank, q-1 )
    
    for( k in seq( from=maxrank, to=0, by=-1 ) ) {
        if( k == 0 ) {
            break
        } else {
            nu1  <- n * WishartMaxPar( n-k, p-k, beta=beta )$centering
            nu2  <- ( n-k+1 )*( p-k+1 ) - nu1
            F    <- ( ell[ k ]/sum( ell[ -(1:k) ] ) )*( nu2/nu1 )

            if( nu2 >= 0 
                && pf( F, nu1, nu2, lower.tail=FALSE ) <= alpha ) {
                    
                break
            }
        }
    }

    rank <- k
    rank
}    
rank.est.fk <- rank.est( rank.est.eigs.fk )
