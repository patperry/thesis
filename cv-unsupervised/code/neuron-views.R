
require( "bcv", "~/lib/R" )

motor <- read.csv( "../data/motor.csv" )

num_neuron    <- 49
num_condition <- 27
num_time      <- 2105/5 + 1

A <- matrix( NA, num_condition*num_time, num_neuron )

for( r in 1:nrow(motor) ) {
    neuron    <- motor$neuron[r]
    condition <- motor$condition[r]
    time      <- motor$time[r]
    response  <- motor$response[r]
    
    i      <- ( condition-1 )*num_time + time/5 + 1
    j      <- neuron
    A[i,j] <- response
}


T <- matrix( NA, num_condition, num_neuron*num_time )

for( r in 1:nrow(motor) ) {
    neuron    <- motor$neuron[r]
    condition <- motor$condition[r]
    time      <- motor$time[r]
    response  <- motor$response[r]
    
    i      <- condition
    j      <- ( neuron-1 )*num_time + time/5 + 1 
    T[i,j] <- response
}

cv.A <- cv.svd( scale(A, center=TRUE, scale=FALSE), 2, 2 )