
noise.est <- function( noise.est.eigs ) {
    function( x, k, n=NA, p=NA, center=TRUE, scale=FALSE, ... ) {
        if( is.matrix( x ) || is.data.frame( x ) ) {
            if( is.character( x ) || is.factor( x ) )
                stop( "'x' cannot be a 'charactor' or 'factor' type." )
            if( k < 0 )
                stop( "'k' must be non-negative.")
            if( !is.na( n ) )
                warning( "Parameter 'n' is ignored." )
            if( !is.na( p ) )
                warning( "Parameter 'p' is ignored." )

            x   <- scale( as.matrix( x ), center, scale )
            n   <- nrow( x )
            p   <- ncol( x )

            if( ( center && k >= min( n - 1, p  ) )
                || ( !center && k >= min( n, p ) ) )
                stop( "'k' must be less than the number of nonzero",
                      " sample covariance eigenvalues." )

            d   <- svd( x )$d
            ell <- d^2 / n
            
            noise.est.eigs( ell, k, n, p, ... )

        } else {
            if( is.character( x ) || is.factor( x ) )
                stop( "'x' cannot be a 'charactor' or 'factor' type." )
            if( k < 0 )
                stop( "'k' must be non-negative.")
            if( k >= min( n, p ) )
                stop( "'k' must be less than the number of nonzero",
                      " sample covariance eigenvalues." )
            if( is.na( n ) )
                stop( "Must supply a value for 'n'." )
            if( is.na( p ) )
                stop( "Must supply a value for 'p'." )
            if( n <= length( x ) )
                stop( "'n' must be greater than or equal to the length of 'x'.")
            if( p <= length( x ) )
                stop( "'p' must be greater than or equal to the length of 'x'.")
            if( !is.missing( center ) )
                warning( "Parameter 'center' is ignored." )
            if( !is.missing( scale ) )
                warning( "Parameter 'scale' is ignored." )
                
            ell <- sort( as.numeric( x ), decreasing=TRUE )
            
            noise.est.eigs( ell, k, n, p, ... )
        }
    }
}

noise.est.eigs.naive <- function( ell, k, n, p ) {
    if( k == 0 ) {
        sigma2 <- mean( ell )
    } else {
        sigma2 <- mean( ell[ -seq_len( k ) ] )
    }
    sigma2
}
noise.est.naive <- noise.est( noise.est.eigs.naive )

# N.M. Faber, L.M.C. Buydens, and G. Kateman. "Aspects of pseudorank 
#    estimation methods based on the eigenvalues of principal component 
#    analysis of random matrices." Chemometrics and Intelligent Laboratory 
#    Systems, 25:203–226, 1994. 
noise.est.eigs.fbk <- function( ell, k, n, p ) {
    sigma2.0 <- noise.est.eigs.naive( ell, k, n, p )
    sigma2   <- sigma2.0/( 1 - k/n )
    sigma2
}
noise.est.fbk <- noise.est( noise.est.eigs.fbk )
