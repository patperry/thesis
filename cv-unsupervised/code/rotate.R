
rotate <- function( x, rows=TRUE, cols=TRUE ) {
    n <- nrow( x )
    p <- ncol( x )
    
    if( rows ) {
        z1 <- matrix( rnorm( n*n ), n, n )
        s1 <- sample( c(-1,1), n, replace=TRUE )
        px <- qr.qy( qr( z1 ), diag( s1, n, n ) %*% x)
    } else{
        px <- x
    }
    
    if( cols ) {
        z2   <- matrix( rnorm( p*p ), p, p )
        s2   <- sample( c(-1,1), p, replace=TRUE )
        pxqt <- t( qr.qy( qr( z2 ), diag( s2, p, p ) %*% t( px ) ) )
    } else {
        pxqt <- px
    }
    
    pxqt
}
