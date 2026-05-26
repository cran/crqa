## written by Moreno I. Coco (moreno.cocoi@gmail.com)
## License: GPL > 2.

## Simulate a Roessler attractor with user-defined parameters
##
## Arguments:
##   numsteps: number of simulated points
##   dt      : integration time step (Euler)
##   a, b, c : Roessler system parameters
##
## System:
##   dx/dt = -y - z
##   dy/dt =  x + a*y
##   dz/dt =  b + z*(x - c)
##
## Standard parameters (Marwan & Kraemer 2023, Eur. Phys. J. Spec. Top.,
## benchmark protocol in Appendix A.3):
##   a = 0.25, b = 0.25, c = 4, dt = 0.05
## with the first 1000 points discarded as transient.

rosslerattractor <- function(numsteps, dt = 0.05, a = 0.25, b = 0.25, c = 4){

  x <- matrix(0, nrow = numsteps, ncol = 3)

  ## small random initial state (mirrors lorenzattractor() convention)
  for (gn in 1:3) x[1, gn] <- rnorm(1)

  x0 <- x[1, ]
  for (i in 2:numsteps) {
    x[i, ]    <- x0
    x[i, 1]   <- x[i, 1] + dt * (-x0[2] - x0[3])
    x[i, 2]   <- x[i, 2] + dt * ( x0[1] + a * x0[2])
    x[i, 3]   <- x[i, 3] + dt * ( b + x0[3] * (x0[1] - c))
    x0        <- x[i, ]
  }

  return(x)
}
