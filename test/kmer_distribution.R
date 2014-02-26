# this is a test script to prove see how xor changes the bit pattern distribution when using 2 bits to represent DNA,

# experiment holds because because 
#    1. complement of a letter is the bit inverse, 
#    2. 2 bits, 
#    3. mutually exclusive letters.

# if probability for p(A) = p(T) = p, then p(C) = p(G) = 0.5 - p
# xor of a character from seq and a char of its reverse complement has patterns 00, 01, 10, 11, each is a aggregate of specific
#   combinations of A,T,C,G.  e.g. 00 = AA, TT, CC, GG. 
#   p(00) = p(11) = 0.5 - 2x + 4x^2,  p(01) = p(10) = 2x - 4x^2
#   quadratic, with min/max at x = 0.25
# compared to y=x, the 2 functions cross at x = 0 and x = 0.25

# to evaluate distribution, 2 sets of plots are analyzed:
#   1. p(s) vs x, x in [0, 0.5].  s is a sequence of 1 .. 6 characters.  choosing 6 because 12 bits cover 4096 processors.
#         this set of experiments show how probably a subsequence with a specific proportion of {00, 11} vs {01, 10} changes as p changes
#         p from real genome, e.g. 0.18 for p(A) = p(T) in human
#   2. p(s) vs proc_id, for number of procs from {4, 16,64, 256, 1024,4096}.  s is set to the processor id (key for distribution)
#         
# Results:
#   1. xor produces an even function for the probability distribution, and is more flat as p changes.  therefore, xor is more robust
#         w.r.t. p changes.  Also effect of p being big (close to 0.5) is less pronounced
#   2. when p = 0.25, as expected, even distribution.  when p <> 0.25, xor version is a lot smoother, but not completely flat.
#
# Conclusion:
#   xor is a good operator to promote load balancing when distributing kmers, under the assumption that 
#      we use 2bit DNA alphabet, and complement can be obtained using bit inverse.


# NOT YET ANALYZED:  using 3 bits (DNA5) or 4 bits (IUPAC)
#       DNA5 issue:  condition 1 fails.
#       IUPAC issue:  condition 3 fails, so some letter's probability are additive. (or alternatively, certain NOT some letters.)


plotfunc <- function(p, y, z, prange) {
  op <- par(mar = par("mar")/2)
  plot(p, y, ylim=prange, type='l', col='#FF0000')
  par(new=T)
  plot(p, z, ylim=prange, axes=F, type='l', col='#0000FF', xlab='', ylab='')  
  par(op)
}

plotloads <- function(p, i) {
  q = 0.5 - 2 * p + 4 * p * p
  
  nprocs <- 1
  for (x in 1:i) {
    nprocs <- nprocs * 4
  }
  
  # allocate
  procs <- array(0:(nprocs - 1))
  
  # now for each proc, determine it's values
  y <- rep(1, nprocs)
  z <- rep(1, nprocs)
  div <- nprocs / 4
  for (x in 1:i) {
    pt <- (procs %/% div) %% 4
    
    y[pt == 0] = y[pt == 0] * p
    y[pt == 1] = y[pt == 1] * (0.5 - p)
    y[pt == 2] = y[pt == 2] * (0.5 - p)
    y[pt == 3] = y[pt == 3] * p

    z[pt == 0] = z[pt == 0] * q
    z[pt == 1] = z[pt == 1] * (0.5 - q)
    z[pt == 2] = z[pt == 2] * (0.5 - q)
    z[pt == 3] = z[pt == 3] * q
        
    div <- div / 4
  }

  r <- range(c(y,z))
  
  
  op <- par(mar = par("mar")/2)
  plot(procs, y, ylim=r, type='l', col='red')
  par(new=T)
  plot(procs, z, ylim=r, axes=F, type='l', col='blue', xlab='', ylab='')  
  par(op)
  
  
}

plotempty <- function() {
  op <- par(mar = par("mar")/2)
  plot.new()
  par(op)
}

p <- array(0:100) * 0.005
pp <- 0.5 - p

q <- 2*p - 4 * p * p
qq <- 0.5 - q



#get( getOption( "device" ) )()
par(mfrow=c(6,6))


# i = 1
r <- c(0, 0.5)
y <- p
z <- q
plotfunc(p, y, z, r)

plotempty()
plotempty()
plotempty()

plotloads(0.25, 1)
plotloads(0.32, 1)

# i = 2
r <- c(0, 0.25)
y <- p * p
z <- q * q
plotfunc(p, y, z, r)

y <- p * pp
z <- q * qq
plotfunc(p, y, z, r)

plotempty()
plotempty()

plotloads(0.25, 2)
plotloads(0.32, 2)


# i = 3
r <- c(0, 0.25)
y <- p * p * p
z <- q * q * q
plotfunc(p, y, z, r)

y <- p * p * pp
z <- q * q * qq
plotfunc(p, y, z, r)

plotempty()
plotempty()

plotloads(0.25, 3)
plotloads(0.32, 3)


# i = 4
r <- c(0, 0.125)
y <- p * p * p * p
z <- q * q * q * q
plotfunc(p, y, z, r)

y <- p * p * p * pp
z <- q * q * q * qq
plotfunc(p, y, z, r)

y <- p * p * pp * pp
z <- q * q * qq * qq
plotfunc(p, y, z, r)

plotempty()

plotloads(0.25, 4)
plotloads(0.32, 4)

# i = 5
r <- c(0, 0.0625)
y <- p * p * p * p * p
z <- q * q * q * q * q
plotfunc(p, y, z, r)

y <- p * p * p * p * pp
z <- q * q * q * q * qq
plotfunc(p, y, z, r)

y <- p * p * p * pp * pp
z <- q * q * q * qq * qq
plotfunc(p, y, z, r)

plotempty()

plotloads(0.25, 5)
plotloads(0.32, 5)

# i = 6
r <- c(0, 0.03125)
y <- p * p * p * p * p * p
z <- q * q * q * q * q * q
plotfunc(p, y, z, r)

y <- p * p * p * p * p * pp
z <- q * q * q * q * q * qq
plotfunc(p, y, z, r)

y <- p * p * p * p * pp * pp
z <- q * q * q * q * qq * qq
plotfunc(p, y, z, r)

y <- p * p * p * pp * pp * pp
z <- q * q * q * qq * qq * qq
plotfunc(p, y, z, r)

plotloads(0.25, 6)
plotloads(0.32, 6)
