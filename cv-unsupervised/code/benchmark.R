
source( "spikesim.R" )
require( "plyr" )

reps <- 100

benchmark <- function( reps=100, n=100, p=50, maxrank=20, 
                       spike=5:10, factors="gaussian", noise="white",
                       factor.seed.start=1, noise.seed.start=10001 ) {
    
    raw <- ldply( seq_len( reps ), function( i ) {
               factor.seed <- factor.seed.start + i - 1
               noise.seed  <- noise.seed.start  + i - 1
           
               sim <- spikesim( spike, n, p, maxrank, factors, factors, noise,
                                factor.seed=factor.seed, noise.seed=noise.seed )
               est <- summary( sim )$est
               
               res <- data.frame( rep=i, 
                                  factor.seed=factor.seed, 
                                  noise.seed=noise.seed,
                                  est )
               res
            })
    
    class( raw ) <- c("benchmark", "data.frame")
    raw
}

summary.benchmark <- function( object, ... ) {
    raw   <- object
    ranks <- raw[, -match( c("rep", "factor.seed", "noise.seed", "true"), 
                       colnames( raw ) ) ] 
    true  <- raw$true

    nbins.pos <- max(1, max( ranks - true ))
    pos <- adply( ranks - true, 2, tabulate, nbins=nbins.pos )
    colnames( pos ) <- c("method", 
                         laply( seq_len( nbins.pos ), function( x) 
                             paste( "+", x, sep='' ) ) )
    pos <- pos[,-1, drop=FALSE ]

    nbins.neg <- max(1, max( true  - ranks ))
    neg <- adply( true - ranks, 2, tabulate, nbins=nbins.neg )
    colnames( neg ) <- c("method", 
                         laply( seq_len( nbins.neg ), function( x) 
                             paste( "-", x, sep='' ) ) )
    neg <- neg[ , c(1, nbins.neg:1 + 1), drop=FALSE ]

    zero <- adply( true - ranks, 2, function( x ) sum( x == 0 ) )
    colnames( zero ) <- c("method", "0")
    zero <- zero[,-1, drop=FALSE]

    res <- cbind( neg, zero, pos )
    res
}

n            <- 100
p            <- 50
maxrank      <- 20
spike.weak   <- 5:10
spike.strong <- 10*spike.weak
reps         <- 100

wgw <- benchmark( reps, n, p, maxrank, spike.weak, "gaussian", "white" )
save( wgw, file="../objects/wgw.rda" )

wgc <- benchmark( reps, n, p, maxrank, spike.weak, "gaussian", "colored" )
save( wgc, file="../objects/wgc.rda" )

wgh <- benchmark( reps, n, p, maxrank, spike.weak, "gaussian", "heavy" )
save( wgh, file="../objects/wgh.rda" )

sgw <- benchmark( reps, n, p, maxrank, spike.strong, "gaussian", "white" )
save( sgw, file="../objects/sgw.rda" )

sgc <- benchmark( reps, n, p, maxrank, spike.strong, "gaussian", "colored" )
save( sgc, file="../objects/sgc.rda" )

sgh <- benchmark( reps, n, p, maxrank, spike.strong, "gaussian", "heavy" )
save( sgh, file="../objects/sgh.rda" )

wsw <- benchmark( reps, n, p, maxrank, spike.weak, "sparse", "white" )
save( wsw, file="../objects/wsw.rda" )

wsc <- benchmark( reps, n, p, maxrank, spike.weak, "sparse", "colored" )
save( wsc, file="../objects/wsc.rda" )

wsh <- benchmark( reps, n, p, maxrank, spike.weak, "sparse", "heavy" )
save( wsh, file="../objects/wsh.rda" )

ssw <- benchmark( reps, n, p, maxrank, spike.strong, "sparse", "white" )
save( ssw, file="../objects/ssw.rda" )

ssc <- benchmark( reps, n, p, maxrank, spike.strong, "sparse", "colored" )
save( ssc, file="../objects/ssc.rda" )

ssh <- benchmark( reps, n, p, maxrank, spike.strong, "sparse", "heavy" )
save( ssh, file="../objects/ssh.rda" )


latex <- function( object, ... ) UseMethod( "latex" )

latex.benchmark <- function( object, ... ) {
    print.counts <- function( row, min=-7, max=7, 
                                   less.bin=FALSE, greater.bin=TRUE ) {
        row.ub <- ncol( row )
        row.lb <- -row.ub
        
        counts  <- rep( 0, max-min+1 )
        greater <- 0
        less    <- 0
        
        for( i in seq( from=row.lb, to=row.ub, by=1 ) ) {
            if( i <= 0 ) {
                col <- match( i, colnames( row ) )
            } else {
                col <- match( paste( "+", i, sep='' ), colnames( row ) )
            }
            
            if( !is.na( col ) ) {
                if( i < min ) {
                    less <- less + row[[ col ]]
                } else if ( i > max ) {
                    greater <- greater + row[[ col ]]
                } else {
                    counts[ i-min+1 ] <- counts[ i-min+1 ] + row[[ col ]]
                }
            }
        }
        
        if( less.bin ) {
            if( less > 0 ) {
                cat( " & ", less )
            } else{
                cat( " & " )
            }
        }
        
        for( c in counts ) {
            if( c > 0 ) {
                cat( " & ", c )
            } else {
                cat( " & " )
            }
            
        }
        
        if( greater.bin ) {
            if( greater > 0 ) {
                cat( " & ", greater )
            } else {
                cat( " & " )
            }
        }
        
        cat( "\\\\ \n")
    }
    
    sum <- summary( object )
    method <- sum$method
    
    cat( " CV-Gabriel");  print.counts( sum[method == "gabriel.best",] )
    cat( " GCV-Gabriel"); print.counts( sum[method == "gabriel.rot.best",] )
    cat( " CV-Wold");     print.counts( sum[method == "wold.best",] )
    cat( " GCV-Wold");    print.counts( sum[method == "wold.rot.best",] )
    cat( " AIC");         print.counts( sum[method == "aic",] )
    cat( " BIC$_1$");     print.counts( sum[method == "bic1",] )
    cat( " BIC$_2$");     print.counts( sum[method == "bic2",] )
    cat( " BIC$_3$");     print.counts( sum[method == "bic3",] )
    cat( " $F$");         print.counts( sum[method == "fk",] )
    cat( " MDL");         print.counts( sum[method == "mdl",] )
    cat( " UIP");         print.counts( sum[method == "kn",] )
}

