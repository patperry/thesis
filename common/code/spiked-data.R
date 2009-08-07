
factors <- function( p, k, type=c("basis", "uniform", "gaussian", "sparse"), 
                     ... ) {
    factors.basis <- function( p, k, ... ) {
        u <- diag( 1, p, k )
        u
    }
    
    factors.uniform <- function( p, k, ... ) {
        z <- matrix( rnorm( p*k ), p, k )
        q <- qr.Q( qr( z ) )
        s <- sample( c(-1,1), k, replace=TRUE )
        u <- q %*% diag( s, k, k )
        u
    }

	factors.gaussian <- function( p, k, ... ) {
		u <- matrix( rnorm( p*k, sd=1/sqrt( p ) ), p, k )
		u
	}
    
    factors.sparse <- function( p, k, prob=0.10, ... ) {
        if( prob < 0 || prob > 1)
            stop( "prob must be between 0 and 1" )
            
        pzero <- 1 - prob
        pone  <- 0.5*prob
        pnone <- pone
        
        u <- matrix( sample( c(-1, 0, 1), p*k, replace=TRUE, 
                             prob=c(pnone, pzero, pone) ), p, k )

        n <- apply( abs( u ), 2, sum )
        
        if( any( n == 0 ) ) {
            u <- factors.sparse( p, k, prob, ... )
        }

        u <- u / sqrt( p*prob )
        u
    }
    
    type <- match.arg( type )
    u <- switch( type,
             basis=factors.basis( p, k, ... ),
           uniform=factors.uniform( p, k, ... ),
	      gaussian=factors.gaussian( p, k, ... ),
	        sparse=factors.sparse( p, k, ... )
         )
    u
}

noise <- function( n, p, type=c("white", "heavy", "colored"), ... ) {
    noise.white <- function( n, p, var=1, ... ) {
        if( var <= 0 )
            stop( "var must be positive" )

        e <- matrix( rnorm( n*p, sd=sqrt( var ) ), n, p )
        e
    }
    
    noise.heavy <- function( n, p, var=1, df=5, ... ) {
        if( var <= 0 )
            stop( "var must be positive" )
        if( df <= 2 )
            stop( "df must be greater than 2")

        t.var <- df/( df-2 )
        scale <- sqrt( var/t.var )
        
        e <- matrix( scale * rt( n*p, df=df ), n, p )
        e
    }
    
    noise.colored <- function( n, p, var=1, rowdf=5, coldf=5, ... ) {
        if( var <= 0 )
            stop( "var must be positive" )
        if( rowdf <= 2 )
            stop( "rowdf must be greater than 2")
        if( coldf <= 2 )
            stop( "coldf must be greater than 2")
        
        row.var <- 1/( rowdf-2 )
        col.var <- 1/( coldf-2 )
        elt.var <- row.var + col.var
        scale   <- sqrt( var/elt.var )
        
        sigma.row <- sqrt( 1/rchisq( n, df=rowdf ) )
        sigma.col <- sqrt( 1/rchisq( p, df=coldf ) )
        e.row   <- matrix( rnorm( n*p, sd=sigma.row ), n, p, byrow=FALSE )
        e.col   <- matrix( rnorm( n*p, sd=sigma.col ), n, p, byrow=TRUE )
        e       <- scale*( e.row+e.col )
        e
    }
    
    type <- match.arg( type )
    e <- switch( type,
             white=noise.white( n, p, ... ),
             heavy=noise.heavy( n, p, ... ),
           colored=noise.colored( n, p, ... ),
         )
    e
}

spiked.data <- function( spike, n, p, 
                         left=c("basis", "uniform", "gaussian", "sparse"), 
                        right=c("basis", "uniform", "gaussian", "sparse"),
                        noise=c("white", "heavy", "colored"),
                    var.noise=1,
                  compute.svd=TRUE,
                         ... ) {
    k <- length( spike )
    u <- factors( n, k, left, ... )
    v <- factors( p, k, right, ... )
    e <- noise( n, p, noise, var=var.noise, ... )

    if( k > 0 ) {
        mu <- u %*% diag( sqrt( n*spike ), k, k ) %*% t( v )
        x  <- mu + e
    } else {
        mu <- matrix( 0, n, p )
        x  <- e
    }
    
    res <- list( n=n, p=p, k=k, u=u, spike=spike, v=v, 
                 x=x, signal=mu, noise=e, var.noise=var.noise )
    
    if( compute.svd ) {
        x.svd <- svd( x )
        res <- c( res, list( u.est=x.svd$u, d.est=x.svd$d, spike.est=(x.svd$d)^2/n, v.est=( x.svd$v ) ) )
    }
    
    res
}

msep.spiked <- function( data ) {
    with( data, {
        if( !( exists( "u.est" ) 
               && exists( "d.est" ) 
               && exists( "v.est" ) ) ) {
            x.svd <- svd( x )
            u.est <- x.svd$u
            d.est <- x.svd$d
            v.est <- x.svd$v 
        }
        
        n     <- nrow( x )
        p     <- ncol( x )
        np    <- min( n, p )
        resid <- signal
        msep  <- rep( NA, np+1 )
        
        msep[ 1 ] <- mean( signal^2 )
        for( k in seq_len( np ) ) {
            resid <- ( resid 
                       - 
                       ( ( d.est[k] * u.est[,k,drop=FALSE] ) 
                         %*% t( v.est[,k,drop=FALSE] ) ) )
            msep[ k+1 ] <- mean( resid^2 )
        }
        
        msep
    })
}

msem.spiked <- function( data ) {
    msep.spiked( data ) + data$var.noise
}
