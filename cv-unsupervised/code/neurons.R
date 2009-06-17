
gs      <- "/opt/local/bin/gs"
mogrify <- "/opt/local/bin/mogrify"

require( "plyr" )
require( "ggplot2" )


motor <- read.csv( "../data/motor.csv" )

plot_neuron <- function( neuron ) {
    null <- function( x ) { "" }

    ( ggplot( neuron, aes( time, response, colour=factor( condition ) ) )
      + geom_line()
      + labs( x=NULL, y=NULL )
      + scale_x_continuous( limits=range( motor$time ),     formatter=null )
      + scale_y_continuous( limits=range( motor$response ), formatter=null )
      + theme_bw() 
      + opts( axis.ticks=theme_blank() ) # theme_segment( size=0 ) )
      + opts( axis.text.x=theme_blank() )
      + opts( axis.text.y=theme_blank() )
      + opts( axis.title.x=theme_blank() )
      + opts( axis.title.y=theme_blank() )
      + opts( panel.margin=unit( rep( 0,4 ), "lines" ) )
      + opts( legend.position="none" ) 
      + opts( plot.margin=unit( rep( 0,4 ), "lines") )
      )
 }
  

pdf( "../plots/neurons.pdf", width=8, height=8 )
grid.newpage()
pushViewport( viewport( layout=grid.layout( 7,7 ) ) )
d_ply( motor, .(neuron), 
    function( neuron, ... ) {
        n <- neuron$neuron[ 1 ]
        c <- ( n-1 ) %%  7 + 1
        r <- ( n-1 ) %/% 7 + 1
        p <- plot_neuron( neuron )
        print( p, vp=viewport( layout.pos.row=r, layout.pos.col=c ) )
    }
)
dev.off()


neuron1 <- motor[ motor$neuron == 1, -1 ]
pdf( "../plots/neuron1.pdf", width=2.75, height=2.75 )
( ggplot( neuron1, aes( time, response, colour=factor( condition ) ) )
      + geom_line()
      + labs( x="Time (ms)", y="Response" )
      + scale_x_continuous( limits=range( motor$time ) )
      + scale_y_continuous( limits=range( motor$response ) )
      + theme_bw() 
      + opts( legend.position="none" ) 
      )
dev.off()


system( paste( gs,
			   "   -dSAFTER -dBATCH -dNOPAUSE -sDEVICE=png16m",
			   "   -dGraphicsAlphaBits=4 -dTextAlphaBits=4 -r600",
   			   "   -dBackgroundColor='16#ffffff'",
   			   "   -sOutputFile=../plots/neurons.png > /dev/null",
   			   "   ../plots/neurons.pdf &&",
			   mogrify,
   			   "   -resize 1800 ../plots/neurons.png", sep=" \\\n" ) )
