
source( "../../common/code/spiked-data" )



plot_scree <- function( d.est, frob2, spec2, elbow=NA ) {
    np     <- length( d.est )
    kmax.f <- length( frob2 ) - 1
    kmax.2 <- length( spec2 ) - 1
    

    frob2.min <- which.min( frob2 ) - 1
    spec2.min <- which.min( spec2 ) - 1
    frob2 <- frob2 / min( frob2 )
    spec2 <- spec2 / min( spec2 )
    
    p <- 
      ( ggplot( data.frame( x=c(0:np), y=0, type=factor("f1", "f2", "f3") ), 
                aes( x, y, colour=type ) ) 
      + facet_grid( type ~ ., scales="free_y") 
      + geom_vline( aes( xintercept=xint, colour=I("f1") ), linetype="dashed",
                    data=data.frame( xint=elbow,
                                     type=factor( c("f1", "f2", "f3") ) ) )
      + geom_vline( aes( xintercept=xint, colour=I("f2") ), linetype="dotdash",
                    data=data.frame( xint=frob2.min,
                                     type=factor( c("f1", "f2", "f3") ) ) )
      + geom_vline( aes( xintercept=xint, colour=I("f3") ), linetype="twodash",
                    data=data.frame( xint=spec2.min,
                                     type=factor( c("f1", "f2", "f3") ) ) )
      + layer( data=data.frame( x=1:np, y=(d.est^2 / sum( d.est^2 ) ), type="f1" ),
               geom=c("point") )
      + layer( data=data.frame( x=0:kmax.f, y=frob2, type="f2" ),
               geom="point" )
      + layer( data=data.frame( x=0:kmax.2, y=spec2, type="f3" ),
               geom="point" )
      + theme_bw()
      + xlab( "Rank" )
      + ylab( paste( "      Resid. Spec. Sq.         ",
                     "Resid. Frob. Sq.      ",
                     "Singular Value Sq.")
            )
      + opts( strip.background=theme_blank(),
              strip.text.y=theme_blank(),
              legend.position="none" )
      )
    p
}


scree_sim <- function( spike, n, p, ... ) {
    np    <- min( n, p )
    sim   <- spiked.data( spike, n, p, 
                          left="uniform", 
                         right="uniform", 
                         noise="white" )

    resid <- sim$signal
    spec2 <- rep( NA, np+1 )
    frob2 <- rep( NA, np+1 )

    if( np > 0 ) {
        spec2[1] <- svd( resid, nu=0, nv=0 )$d[1]^2
        frob2[1] <- sum( resid^2 )
    }

    for( i in seq_len( np ) ) {
        u <- sim$u.est[,i,drop=FALSE]
        v <- sim$v.est[,i,drop=FALSE]
        d <- sim$d.est[i]
    
        resid <- resid - (d * u) %*% t(v)
    
        spec2[i+1] <- svd( resid, nu=0, nv=0 )$d[1]^2
        frob2[i+1] <- sum( resid^2 )
    }

    frob2 <- frob2
    spec2 <- spec2

    plot_scree( sim$d.est, frob2, spec2, ... )
}


n <- 100
p <- n
spike.right <- seq(5, 0.25, length=20)
spike.left  <- c(20, 15, 10, spike.right)
elbow.right <- 13
elbow.left  <- 4

set.seed( 0, "Mersenne-Twister" )
pdf( "../plots/scree-elbow-right.pdf", width=3, heigh=5.75 )
print( scree_sim( spike.right, n, p, elbow=elbow.right ) )
dev.off()

set.seed( 2, "Mersenne-Twister" )
pdf( "../plots/scree-elbow-left.pdf", width=3, heigh=5.75 )
print( scree_sim( spike.left, n, p, elbow=elbow.left ) )
dev.off()
