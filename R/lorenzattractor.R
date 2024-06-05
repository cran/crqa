## written by Moreno I. Coco (moreno.cocoi@gmail.com)
## License: GPL > 2.

## Simulate a lorenz attractor with user-defined parameters

## Arguments:

## numsteps: the number of simulated points 
## dt, sigma, b, r: the parameters shaping the motion

## possible initial values (standard parameters, Lorenz, 1963)
# numsteps = 2 ^ 11; dt = .01; sigma = 10; r = 28; b = 8/3;

lorenzattractor <- function(numsteps, dt, sigma, r, b){
    
#    require(plot3D)
    
    ## initialize the matrix
    x = matrix(0, nrow = numsteps, ncol = 3)
    
    for (gn in 1:3) x[1, gn] = 10*rnorm(1) 

    x0 = x[1,]
    for (i in 2:numsteps){ ## fill in the matrix
        # print(x0)
        
        x[i, ] = x0       
        x[i,1] = x[i,1] + dt*sigma*(x0[2] - x0[1])
        x[i,2] = x[i,2] + dt*(-x0[2] + x0[1] * (r - x0[3]))
        x[i,3] = x[i,3] + dt*(x0[1]*x0[2]- b*x0[3])
        x0 = x[i,]
    }

    return(x)
    
}
