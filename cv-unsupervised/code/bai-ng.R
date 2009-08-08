
source( "rank.est.R" )

info.crit.bn <- function( ell, n, p, maxrank, type=c("1", "2", "3") ) {
    type <- match.arg( type )

    if( length( ell ) < maxrank )
        ell <- c( ell, rep( 0, maxrank-length( ell ) ) )
        
    v    <- ( sum( ell ) - c( 0, cumsum( ell ) ) )/( n*p )
    c.np <- min( sqrt( n ), sqrt( p ) )
    k    <- 0:maxrank
    
    
    if( type == "1" ) {
        ic <- log( v ) + k * ( (n+p)/(n*p) ) * log( ( n*p )/( n+p ) )
    } else if( type == "2" ) {
        ic <- log( v ) + k * ( (n+p)/(n*p) ) * log( c.np^2 )
    } else {
        ic <- log( v ) + k * ( log( c.np^2 )/( c.np^2 ) )
    }
    
    ic
}

rank.est.eigs.bic <- function( type ) {
    function( ell, n, p, maxrank ) {
        ic   <- info.crit.bn( ell, n, p, maxrank, type )
        rank <- which.min( ic ) - 1
        rank
    }
}

rank.est.bic1 <- rank.est( rank.est.eigs.bic( "1" ) )
rank.est.bic2 <- rank.est( rank.est.eigs.bic( "2" ) )
rank.est.bic3 <- rank.est( rank.est.eigs.bic( "3" ) )
