
source( "rank.est.R" )

rank.est.eigs.mF <- function( ell, n, p, maxrank, alpha=0.05 ) {

    q <- min( n, p )
    F <- rep( NA, maxrank )
    
    for( k in seq( from=1, to=maxrank, by=1 ) ) {
        j    <- seq( from=k+1, to=p, by=1 )
        jq   <- seq( from=k+1, to=q )
        F[k] <- ( ( ell[k] * sum( ( n-j+1 )*( p-j+1 ) ) )
                  / ( sum(  ell[ jq ] )*( n-k+1 )*( p-k+1 ) ) )
    }

    for( k in seq( from=maxrank, to=0, by=-1 ) ) {
        if( k == 0 )
            break
        
        if( pf( F[k], 1, q-k, lower.tail=FALSE ) <= alpha )
            break
    }
    
    rank <- k
    rank
}    
rank.est.mF <- rank.est( rank.est.eigs.mF )
