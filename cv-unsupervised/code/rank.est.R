
rank.est <- function( rank.est.eigs ) {
    function( x, n=NA, p=NA, center=TRUE, scale=FALSE, maxrank=min(n,p)-1, ... ) {
        if( is.matrix( x ) || is.data.frame( x ) ) {
            if( is.character( x ) || is.factor( x ) )
                stop( "'x' cannot be a 'charactor' or 'factor' type." )
            if( !is.na( n ) )
                warning( "Parameter 'n' is ignored." )
            if( maxrank < 0 )
                stop( "'maxrank' must be non-negative.")
            if( !is.na( p ) )
                warning( "Parameter 'p' is ignored." )
            
            x   <- scale( as.matrix( x ), center, scale )
            n   <- nrow( x )
            p   <- ncol( x )

            d   <- svd( x )$d
            ell <- d^2 / n
            
            rank.est.eigs( ell, n, p, maxrank, ... )
        } else {
            if( is.character( x ) || is.factor( x ) )
                stop( "'x' cannot be a 'charactor' or 'factor' type." )
            if( is.na( n ) )
                stop( "Must supply a value for 'n'." )
            if( is.na( p ) )
                stop( "Must supply a value for 'p'." )
            if( n <= length( x ) )
                stop( "'n' must be greater than or equal to the length of 'x'.")
            if( p <= length( x ) )
                stop( "'p' must be greater than or equal to the length of 'x'.")
            if( maxrank < 0 )
                stop( "'maxrank' must be non-negative.")
            if( !is.missing( center ) )
                warning( "Parameter 'center' is ignored." )
            if( !is.missing( scale ) )
                warning( "Parameter 'scale' is ignored." )
                
            ell <- as.numeric( x )
            
            rank.est.eigs( ell, n, p, maxrank, ... )
        }
    }
}
