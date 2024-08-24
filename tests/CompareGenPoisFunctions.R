#==================================================================================#
# Purpose: compare the GenPois we are using with the genpois from vgam package
#
# Notes: Our implementation produces NA in the inverse function for positive alpha,
#        that VGAM computes. Also vgam doesnt allow for negative alpha, giving
#        NA, but we compute something that looks wrong.
#==================================================================================#

library(VGAM)

# a = 0
p = 0.9
m = 0.1:120
a = 0
qgenpois2(p,m,a)
qGpois(p,a,m)

x = 0:4
m = 3
dgenpois2(x,m,a)
dGpois(x,a,m)

q = seq(-1,21,by=0.2)
pgenpois2(q,m,a)
pGpois(q,a,m)

# positive a
p = 0.9
m = 0.1:120
a = 0.2
qgenpois2(p,m,a)
qGpois(p,a,m)

x = 0:4
m = 3
dgenpois2(x,m,a)
dGpois(x,a,m)

q = seq(-1,21,by=0.2)
pgenpois2(q,m,a)
pGpois(q,a,m)

# negative a
p = 0.9
m = 0.1:120
a = -0.5
qgenpois2(p,m,a)
qGpois(p,a,m)

x = 0:4
m = 3
dgenpois2(x,m,a)
dGpois(x,a,m)

q = seq(-1,21,by=0.2)
pgenpois2(q,m,a)
pGpois(q,a,m)
