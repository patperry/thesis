
require( "ggplot2" )
require( "plyr" )

m <- function( x, y, gamma ) {
    my.sqrt <- function( z, a, b ) {
        res.1    <-  sqrt( ( z-a )*( z-b ) )
        ifelse( sign( Im( z ) ) == sign( Im( res.1 ) ),
                res.1, -res.1 )
    }

    z   <- complex( real=x, imag=y )
    a   <- ( 1 - 1/sqrt( gamma ) )^2
    b   <- ( 1 + 1/sqrt( gamma ) )^2
    res <- ( -( z + 1 - 1/gamma ) + my.sqrt( z, a, b ) )/( 2*z )

    list( x=Re( res ), y=Im( res ) )
}

support <- function( gamma ) {
    a    <- ( 1 - 1/sqrt( gamma ) )^2
    b    <- ( 1 + 1/sqrt( gamma ) )^2
    n.0  <- 5000

    x <- seq( from=a, to=b, len=n.0 )
    y <- rep( 0, n.0 )

    m.xy <- m( x, y, gamma )
    res  <- data.frame( real=m.xy$x, imag=m.xy$y )
    res
}

contour <- function( gamma ) {
    a    <- ( 1 - 1/sqrt( gamma ) )^2
    b    <- ( 1 + 1/sqrt( gamma ) )^2
    xmax <- b + 0.01
    xmin <- a - 0.01
    ymin <- -0.1
    ymax <-  0.1

    n.1  <-   50
    n.2  <-  500
    n.3  <- 1000
    n.4  <-  n.2

    x <- c( rep( xmax, n.1 ),
            seq( from=xmax, to=xmin, len=n.2 ),
            rep( xmin, n.3 ),
            seq( from=xmin, to=xmax, len=n.4 ) )
    y <- c( seq( from=ymin, to=ymax, len=n.1 ),
            rep( ymax, n.2 ),
            seq( from=ymax, to=ymin, len=n.3 ),
            rep( ymin, n.4 ) )
    part <- c( rep( 1, n.1 ),
               rep( 2, n.2 ),
               rep( 3, n.3 ),
               rep( 4, n.4 ) )

    m.xy <- m( x, y, gamma )
    res  <- data.frame( real=m.xy$x, imag=m.xy$y, part=factor( part ) )
    res
}

make.paths <- function( type, gamma=c( 0.2, 5 ) ) {
    param <- data.frame( gamma=c( 0.2, 5 ) )
    within( ddply( param, .(gamma), splat( type ) ), {
                 gamma <- factor( gamma, levels=sort( gamma, dec=TRUE ) )
    } )
}
supp <- make.paths( support )
path <- make.paths( contour ) 

                 
p <- ( ggplot( supp, aes( real, imag ) )
       + facet_grid( gamma ~ . )
       + geom_path( aes( colour=I( "red" ) ) )
       + geom_path( aes( colour=part ), data=path )
       + scale_colour_hue( "Contour", 
             breaks=c("red", 1, 2, 3, 4),
             labels=c("Support", "I", "II", "III", "IV") )
       + labs( x="Real Part", y="Imaginary Part" )
       + theme_bw() 
       )
p

pdf( file="../plots/wishart-clt-contour.pdf", width=6, height=4 )
print( p )
dev.off()
