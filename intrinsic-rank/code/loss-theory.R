
require( "plyr" )


loss.theory <- function( kmax, spike, n=NA, p=NA, concentration=n/p, var=1, ... ) {
    loss.theory0 <- function( k, spike, n=NA, p=NA, concentration=n/p, var=1 ) {
        spike   <- sort( spike, decreasing=TRUE )
        gamma   <- concentration
        n.spike <- length( spike )
    
        loss.theory.block <- function( i ) {
            mu <- ifelse( i <= n.spike, spike[i], 0 )
            if( i <= k ) {
                if( mu > 1/sqrt( gamma ) ) {
                    tr  <- ( ( var/( gamma*mu ) ) 
                             * ( 3*var + ( gamma+1 )*mu ) )
                    det <- ( ( var/( gamma*mu ) )^2
                             * ( mu + var )
                             * ( gamma*mu + var ) )
                } else {
                    tr  <- mu + var * ( 1 + 1/sqrt( gamma ) )^2
                    det <- mu * var * ( 1 + 1/sqrt( gamma ) )^2
                }
            } else {
                tr  <- mu
                det <- 0
            }
        
            frob2 <- tr
            spec2 <- ( tr + sqrt( tr^2 - 4*det ) )/2
        
            data.frame( frob2=frob2, spec2=spec2 )
        }

        block <- ldply( seq_len( max( k, n.spike ) ), loss.theory.block )
        spec2 <- max( block$spec2 )
        frob2 <- sum( block$frob2 )

        data.frame( spec2=spec2, frob2=frob2 )
    }
    
    k <- seq_len( kmax+1 ) - 1
    ldply( k, function( k ) 
        cbind( k=k, 
               loss.theory0( k, spike=spike, concentration=concentration, 
                             var=var ) ) )
}
