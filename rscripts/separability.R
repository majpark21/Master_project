# from:
# https://stats.stackexchange.com/a/78855/43610

separability.measures <- function ( Vector.1 , Vector.2 ) {
  # convert vectors to matrices in case they are not
  Matrix.1 <- as.matrix (Vector.1)
  Matrix.2 <- as.matrix (Vector.2)
  # define means
  mean.Matrix.1 <- mean ( Matrix.1 )
  mean.Matrix.2 <- mean ( Matrix.2 )
  # define difference of means
  mean.difference <- mean.Matrix.1 - mean.Matrix.2
  # define covariances for supplied matrices
  cv.Matrix.1 <- cov ( Matrix.1 )
  cv.Matrix.2 <- cov ( Matrix.2 )
  # define the halfsum of cv's as "p"
  p <- ( cv.Matrix.1 + cv.Matrix.2 ) / 2
  # --%<------------------------------------------------------------------------
  # calculate the Bhattacharryya index
  bh.distance <- 0.125 *t ( mean.difference ) * p^ ( -1 ) * mean.difference +
    0.5 * log (det ( p ) / sqrt (det ( cv.Matrix.1 ) * det ( cv.Matrix.2 )
    )
    )
  # --%<------------------------------------------------------------------------
  # calculate Jeffries-Matusita
  # following formula is bound between 0 and 2.0
  jm.distance <- 2 * ( 1 - exp ( -bh.distance ) )
  # also found in the bibliography:
  # jm.distance <- 1000 * sqrt (   2 * ( 1 - exp ( -bh.distance ) )   )
  # the latter formula is bound between 0 and 1414.0
  # --%<------------------------------------------------------------------------
  # calculate the divergence
  # trace (is the sum of the diagonal elements) of a square matrix
  trace.of.matrix <- function ( SquareMatrix ) {
    sum ( diag ( SquareMatrix ) ) }
  # term 1
  divergence.term.1 <- 1/2 * trace.of.matrix (( cv.Matrix.1 - cv.Matrix.2 ) * 
                                                ( cv.Matrix.2^ (-1) - cv.Matrix.1^ (-1) )
  )
  # term 2
  divergence.term.2 <- 1/2 * trace.of.matrix (( cv.Matrix.1^ (-1) + cv.Matrix.2^ (-1) ) *
                                                ( mean.Matrix.1 - mean.Matrix.2 ) *
                                                t ( mean.Matrix.1 - mean.Matrix.2 )
  )
  # divergence
  divergence <- divergence.term.1 + divergence.term.2
  # --%<------------------------------------------------------------------------
  # and the transformed divergence
  transformed.divergence  <- 2 * ( 1 - exp ( - ( divergence / 8 ) ) )
  indices <- data.frame(
    jm=jm.distance,bh=bh.distance,div=divergence,tdiv=transformed.divergence)
  return(indices)
}