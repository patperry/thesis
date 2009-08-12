
require( "bcv", "~/lib/R" )
source( "../../common/code/spiked-data.R" )
source( "rotate.R" )

source( "bai-ng.R" )
source( "faber-kowalski.R" )
source( "kritchman-nadler.R" )
source( "malinowski.R" )
source( "perry-wolfe.R" )
source( "rao-edelman.R" )
source( "wax-kailath.R" )

spikesim <- function( spike, n, p, 
                      maxrank=min(n,p),
                         left=c("basis", "uniform", "gaussian", "sparse"), 
                        right=c("basis", "uniform", "gaussian", "sparse"),
                        noise=c("white", "heavy", "colored"),
                    var.noise=1,
                  factor.seed=NA,
                   noise.seed=NA ) {
    data <- spiked.data( spike, n, p, left, right, noise, var.noise,
                         factor.seed=factor.seed, noise.seed=noise.seed )
    rank  <- 0:maxrank
    msem  <- msem.spiked( data )[ rank+1 ]
    x     <- data$x
    x.rot <- rotate( x )
    
    cvw     <- cv.svd.wold( x,     maxrank=maxrank )
    cvw.rot <- cv.svd.wold( x.rot, maxrank=maxrank )
    cvg     <- cv.svd.gabriel( x,     maxrank=maxrank )
    cvg.rot <- cv.svd.gabriel( x.rot, maxrank=maxrank )

    true    <- rank[ which.min( msem ) ]
    
    gabriel.best     <- summary( cvg )$rank.best
    gabriel.1se      <- summary( cvg )$rank.1se
    gabriel.rot.best <- summary( cvg.rot )$rank.best
    gabriel.rot.1se  <- summary( cvg.rot )$rank.1se

    wold.best     <- summary( cvw )$rank.best
    wold.1se      <- summary( cvw )$rank.1se
    wold.rot.best <- summary( cvw.rot )$rank.best
    wold.rot.1se  <- summary( cvw.rot )$rank.1se

    aic <- rank.est.aic( x, maxrank=maxrank, center=FALSE )
    
    bic1 <- rank.est.bic1( x, maxrank=maxrank, center=FALSE )
    bic2 <- rank.est.bic2( x, maxrank=maxrank, center=FALSE )
    bic3 <- rank.est.bic3( x, maxrank=maxrank, center=FALSE )
    
    fk <- rank.est.fk( x, maxrank=maxrank, center=FALSE )
    kn <- rank.est.kn( x, maxrank=maxrank, center=FALSE )
    mF <- rank.est.mF( x, maxrank=maxrank, center=FALSE )
    
    pw.naive <- rank.est.pw.naive( x, maxrank=maxrank, center=FALSE )
    pw.fbk   <- rank.est.pw.fbk( x, maxrank=maxrank, center=FALSE )
    pw.kn    <- rank.est.pw.kn( x, maxrank=maxrank, center=FALSE )
    
    mdl <- rank.est.mdl( x, maxrank=maxrank, center=FALSE )
    
    res <- list(  msem=msem,
                  rank=rank,
                  
                   cvg=cvg,
               cvg.rot=cvg.rot,
                   cvw=cvw,
               cvw.rot=cvw.rot,
         
                  true=true, 
        
          gabriel.best=gabriel.best, 
           gabriel.1se=gabriel.1se, 
      gabriel.rot.best=gabriel.rot.best, 
       gabriel.rot.1se=gabriel.rot.1se, 
       
             wold.best=wold.best, 
              wold.1se=wold.1se, 
         wold.rot.best=wold.rot.best, 
          wold.rot.1se=wold.rot.1se, 

                   aic=aic,
                   
                  bic1=bic1,
                  bic2=bic2,
                  bic3=bic3,
      
                    fk=fk,
                    kn=kn,
                    mF=mF,
                    
              pw.naive=pw.naive, 
                pw.fbk=pw.fbk, 
                 pw.kn=pw.kn, 
                 
                   mdl=mdl
        )
      
    class( res ) <- "spikesim"
    res
}



summary.spikesim <- function( object, ... ) {
    rank     <- c()
    type     <- c()
    rotation <- c()
    msep     <- c()
    se       <- c()
    
    nrank    <- length( object$rank )
    
    rank     <- c(rank, object$rank)
    type     <- c(type,     rep( "gabriel", nrank ))
    rotation <- c(rotation, rep( "none", nrank ))
    msep     <- c(msep, summary( object$cvg )$msep.mean)
    se       <- c(se,   summary( object$cvg )$msep.se)
    
    rank     <- c(rank, object$rank)
    type     <- c(type,     rep( "gabriel", nrank ))
    rotation <- c(rotation, rep( "rotated", nrank ))
    msep     <- c(msep, summary( object$cvg.rot )$msep.mean)
    se       <- c(se,   summary( object$cvg.rot )$msep.se)

    rank     <- c(rank, object$rank)
    type     <- c(type,     rep( "wold", nrank ))
    rotation <- c(rotation, rep( "none", nrank ))
    msep     <- c(msep, summary( object$cvw )$msep.mean)
    se       <- c(se,   summary( object$cvw )$msep.se)
    
    rank     <- c(rank, object$rank)
    type     <- c(type,     rep( "wold", nrank ))
    rotation <- c(rotation, rep( "rotated", nrank ))
    msep     <- c(msep, summary( object$cvw.rot )$msep.mean)
    se       <- c(se,   summary( object$cvw.rot )$msep.se)
    
    cv <- data.frame( rank=rank, 
                      type=factor( type ),
                  rotation=factor( rotation ),
                      msep=msep,
                        se=se )
    
    truth <- data.frame( rank=object$rank,
                         msem=object$msem )
    
    est <- data.frame( 
                  true=object$true, 
        
          gabriel.best=object$gabriel.best, 
           gabriel.1se=object$gabriel.1se, 
      gabriel.rot.best=object$gabriel.rot.best, 
       gabriel.rot.1se=object$gabriel.rot.1se, 
       
             wold.best=object$wold.best, 
              wold.1se=object$wold.1se, 
         wold.rot.best=object$wold.rot.best, 
          wold.rot.1se=object$wold.rot.1se, 

                   aic=object$aic,
                   
                  bic1=object$bic1,
                  bic2=object$bic2,
                  bic3=object$bic3,
      
                    fk=object$fk,
                    kn=object$kn,
                    mF=object$mF,
                    
              pw.naive=object$pw.naive, 
                pw.fbk=object$pw.fbk, 
                 pw.kn=object$pw.kn, 
                 
                   mdl=object$mdl )
                         
    list( cv=cv, truth=truth, est=est )
}

plot.spikesim <- function( x, type=c("gabriel", "wold"), ... ) {
    type <- match.arg( type )
    with( x, {
        ylim <- range( c(as.numeric( cvg$msep ),
                         as.numeric( cvg.rot$msep ),
                         as.numeric( cvw$msep ), 
                         as.numeric( cvw.rot$msep ), 
                         as.numeric( msem ) ) )
        
        plot( rank, msem, ylim=ylim, t='l', lty=1, col='red' )
        
        if( type == "gabriel" ) {
            plot( cvg, add=TRUE, lty=2 )
            plot( cvg.rot, add=TRUE, lty=5 )
        } else {
            plot( cvw, add=TRUE, col="green", lty=2 )
            plot( cvw.rot, add=TRUE, col="green", lty=4 )
        }
    } )
}
