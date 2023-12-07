library(rfracture)
library(rgl)

source("lib.R")

pres = function(x,sel) {
  p = matrix(0,n,n)
  p[sel] = x
  p
}

deform = function(p) {
  Re(fft(kern*fft(p),inverse = TRUE)/(n*n))
}

fun = function(d0) {
  w = f1 - f2 - d0
  sel = w < 0
  if (any(sel)) {
    for (iter in 1:100) {
      A = function(x) {
        deform(pres(x,sel))[sel]
      }
      ret = cg(A, w[sel])
      p = pres(ret$x,sel)
      d = deform(p)
      wd = w-d
      wd[sel] = 0
      print(c(mean(sel),range(wd),range(p)))
      nsel = p < 0 | wd < 0
      if (all(sel == nsel)) break
      sel = nsel
    }
  } else {
    p = matrix(0,n,n)
    wd = w
  }
  p = -p*pi
  list(p=p, w=w, wd=wd, sump=mean(p),frac=mean(sel))
}


n = 256
L = 1

k = expand.grid(k1=seq_circ(n),k2=seq_circ(n))
kern = 1/sqrt(k$k1^2+k$k2^2)
kern[1] = 0
kern = matrix(kern, n, n)

x = seq(0,L,len=n+1)[-(n+1)]
points = expand.grid( x = x, y = x )

if (TRUE) {
  froll = 1
  alpha = 4
  frac = fracture_matrix(dims = c(n,n), power.iso=function(f) 0.0001*ifelse(f<froll, 1,f/froll)^{-alpha}, corr.profile = function(l) 0, seed=124)
  frac_s = sqrt(frac$var.diff)
  s = 0.2
  f1 = frac$f1/frac_s*s
  f2 = frac$f2/frac_s*s
}
if (FALSE) {
  f1 = matrix(0,n,n)
  f2 = f1
  spheres = data.frame(x=c(0.25,0.75),y=c(0.25,0.75),r=0.25*sqrt(2))
  spheres = do.call(rbind, apply(expand.grid(-1:1,-1:1),1,function(x) data.frame(x=spheres$x+x[1],y=spheres$y+x[2],r=spheres$r) ))
  for (i in seq_len(nrow(spheres))) {
    h = spheres$r[i]^2 - (points[,1]-spheres$x[i])^2 - (points[,2]-spheres$y[i])^2
    h = pmax(h,0)
    f2[] = pmax(f2[], sqrt(h))
  }
  f2[is.nan(f2)] = 0
}

if (FALSE) {
  f1 = matrix(0,n,n)
  f2 = f1
  h = 0.3^2 - (points[,1]-0.5)^2 - (points[,2]-0.5)^2
  h = ifelse(h>0, 1, 0)
  f2[] = h
}


clear3d()
surface3d(x,x,f1,asp="iso", col="white")
surface3d(x,x,f2,asp="iso", col="white")

K = 30
w = f1-f2
d = seq(min(w),mean(w),len=K)
d = mean(w) + (min(w)-mean(w))*10^seq(0,-3,len=K)
ret = lapply(d, fun)

sump = sapply(ret, function(x)x$sump)
vol = -d
plot(vol,sump,log="x")
lines(vol, s*log((1.5*(vol/s))^-2+1))


plot(vol,sump,log="x")
# sump = s*log((1.5*(vol/s))^-2+1)
# s/1.5/sqrt(exp(sump/s) - 1) = vol
lines(s/(1.5*sqrt(exp(sump/s) - 1)), sump)




perc = sapply(ret, function(x)x$frac)
plot(d,perc,log="y")

vol = sapply(ret, function(x) mean(x$wd))
plot(-d,vol,log="")
abline(0,1)


plot(log((d-d[1])/(-d)),-log(sump),log="")
abline(4,-1.2)
abline(1.5,-0.25)

lnd = seq(-2,7,len=100)
plot(log((d-d[1])/(-d)),-log(sump),log="")
lines(lnd, 3.5+lnd*-1.2)
lines(lnd, 1.5+lnd*-0.25)

lines(lnd, log((exp(3.5)*exp(lnd)^-1 + exp(1.5))*exp(lnd)^-0.25))
lines(lnd, log((30*exp(lnd)^-1 + 5)*exp(lnd)^-0.25))

plot(-d,1/sump,log="xy")
#lines(-d, ((-25*d-d[1]*5)/(d-d[1])^1.25)*(-d)^0.25)
lines(-d, ((25*(-d)^1.25-d[1]*5*(-d)^0.25)/(d-d[1])^1.25))

plot(vol,sump)
plot(diff(sump)/diff(vol))
plot(vol,sump,log="x")
abline(-0.5227,-0.4418)

plot(log10(vol),sump)
lines(log10(vol), -0.5227-0.4418*log10(vol))
smax = function(x,y,k) log(exp(x*k)+exp(y*k))/k

#lines(log10(vol), smax(-0.5227-0.4418*log10(vol),0,10))

plot(log10(vol),sump)
lines(log10(vol), log(0.00536*vol^(-1.9187)+1)/10)

plot(vol,sump)
lines(vol, 1/10*log(0.00536*vol^(-1.9187)+1))

plot(vol,sump)
lines(vol, 0.1*log((15.26*vol)^-2+1))



sel = 10:30
lm(sump[sel]~log10(vol[sel]))


plot(d-d[1],sump,log="xy")
abline(0,2)
abline(-0.9,-1/4)
plot(vol,sump,log="")
lines(vol,1/(vol^{1/4}*10^{0.9} + (d-d[1])^2))

plot(sump,perc,log="xy",asp=1)
abline(0,1)

perc2 = sapply(seq(min(w),max(w),len=length(d)), function(d) mean(w-d<0))
plot(d,perc2)

plot(d,perc)

clear3d()
val = ret[[20]]
w = f1 - f2
wd = val$wd
p = -val$p
surface3d(x,x,f1+(wd-w)/2,asp="iso", col=colrange(p,c("white","red","blue")))
surface3d(x,x,f2-(wd-w)/2,asp="iso", col=colrange(p,c("white","red","blue")))


to_n_plus_1 = function(f) as.vector(rbind(cbind(f,f[,1]),cbind(f[1,,drop=FALSE],f[1,1])))

f = f2
points = expand.grid(
      x = seq(0,L,len=n+1),
      y = seq(0,L,len=n+1))

I = matrix(seq_len((n+1)*(n+1)),n+1,n+1)
triangles = rbind(cbind(
  as.vector(I[1:n,1:n]),
  as.vector(I[1:n+1,1:n]),
  as.vector(I[1:n,1:n+1])
),cbind(
  as.vector(I[1:n+1,1:n]),
  as.vector(I[1:n+1,1:n+1]),
  as.vector(I[1:n,1:n+1])
))

w = f1 - f2
bias = 0.5
for (i in seq_along(ret)) {
  print(i)
  val = ret[[i]]
  wd = val$wd
  p = -val$p
  points$z = to_n_plus_1(f1+(wd-w)*bias)
  write.vtk.tri(
    points,
    triangles,
    point_data = list(p=to_n_plus_1(p)),
    filename = sprintf("vtk_output/frac_A_%04d.vtk",i)
  )
  points$z = to_n_plus_1(f2-(wd-w)*(1-bias))
  write.vtk.tri(
    points,
    triangles,
    point_data = list(p=to_n_plus_1(p)),
    filename = sprintf("vtk_output/frac_B_%04d.vtk",i)
  )
}



n = 256
L = 1

k = expand.grid(k1=seq_circ(n),k2=seq_circ(n))
kern = 1/sqrt(k$k1^2+k$k2^2)
kern = 1/sqrt(k$k1^2+k$k2^2)^3
kern[1] = 0
kern = matrix(kern, n, n)

kern1 = fft(kern, inverse = TRUE)/(n*n)

plot(1:n-1,Re(kern1)[,1],log="")
abline(2.5,-1)
abline(1,-2)
abline(1,-3)

x = (1:n-1)/n
y = Re(kern1)[,1]
y = (max(y)-y)/(max(y)-min(y))
plot(x,y,log="")
lines(x, (1-x)*x*4,col=2)

y = Re(kern1)
y = (max(y)-y)/(max(y)-min(y[1,]))
x = (1:n-1)/n
x = expand.grid(x,x)
ny = 4*(x[,1]*(1-x[,1])+x[,2]*(1-x[,2]))
ny = matrix(ny,n,n)
range((y-ny))
plot(y[1,])
lines(ny[1,])

plot((y-ny)[128,])

plot(as.vector(y)[seq(1,n*n,n+1)])
abline(0,4/n)
range(Re(kern1))
image(y)


a = seq(0,2*pi,len=100)
plot(a,(cos(a)-1)^2+sin(a)^2)
plot(y[,1])
