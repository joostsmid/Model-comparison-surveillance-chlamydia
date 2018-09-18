
theta.f.inf <- theta[,((9+(nyears+2)*6)+1):(9+(nyears+2)*6*2)]
plot3d(1:14,1:6,colMeans(theta.f.inf))

theta.m.inf <- theta[,(9+1):(9+(nyears+2)*6)]
plot3d(1:14,1:6,colMeans(theta.m.inf))
