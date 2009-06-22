
factors <- function( p, k, type=c("basis", "uniform", "gaussian") ) {
    factors.basis <- function( p, k ) {
        u <- diag( 1, p, k )
        u
    }
    
    factors.uniform <- function( p, k ) {
        z <- matrix( rnorm( p*k ), p, k )
        q <- qr.Q( qr( z ) )
        s <- sample( c(-1,1), k, replace=TRUE )
        u <- q %*% diag( s, k, k )
        u
    }

	factors.gaussian <- function( p, k ) {
		u <- matrix( rnorm( p*k, sd=1/p ), p, k )
		u
	}
    
    type <- match.arg( type )
    u <- switch( type,
             basis=factors.basis( p, k ),
           uniform=factors.uniform( p, k ),
	      gaussian=factors.gaussian( p, k )
         )
    u
}

noise <- function( n, p, type=c("white"), var=1 ) {
    noise.white <- function( n, p, var ) {
        e <- matrix( rnorm( n*p, sd=sqrt( var ) ), n, p )
        e
    }
    
    type <- match.arg( type )
    e <- switch( type,
             white=noise.white( n, p, var )
         )
    e
}

spiked.data <- function( spikes, n, p, 
                         left=c("basis", "uniform", "gaussian"), 
                        right=c("basis", "uniform", "gaussian"),
                        noise=c("white"),
                          var=1,
                  compute.svd=TRUE ) {
    k <- length( spikes )
    u <- factors( n, k, left )
    v <- factors( p, k, right )
    e <- noise( n, p, noise, var )

    if( k > 0 ) {
        x <- u %*% diag( sqrt( n*spikes ), k, k ) %*% t( v ) + e
    } else {
        x <- e
    }
    
    res <- list( n=n, p=p, k=k, u=u, spikes=spikes, v=v, e=e, x=x )
    
    if( compute.svd ) {
        x.svd <- svd( x )
        res <- c( res, list( u.est=x.svd$u, spikes.est=(x.svd$d)^2/n, v.est=t( x.svd$v ) ) )
    }
    
    res
}
