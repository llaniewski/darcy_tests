library(rfracture)

froll = 1
alpha = 3.5
n = 256
frac = fracture_matrix(dims = c(n,n), power.iso=function(f) 0.0001*ifelse(f<froll, 1,f/froll)^{-alpha},gap=0.1, corr.profile = function(l) 0)

f1 = frac$f1
f2 = frac$f2
clear3d()
surface3d(x,x,f1,asp="iso", col="white")
surface3d(x,x,f2,asp="iso", col="white")


cg = function(A,b,x0=rep(0,length(b)),itmax=1000,reslim=1e-9) {
  x = x0
  r = b - A(x)
  q = r
  for (i in 1:itmax) {
    res = sqrt(sum(r*r))
    if (res < reslim) break
    Aq = A(q)
    alpha = sum(r * q)/sum(q * Aq)
    x = x + alpha*q
    r = r - alpha*Aq
    beta = sum(r * Aq)/sum(q * Aq)
    q = r - beta * q
  }
  list(x=x, b=b, res=res, iter=i)
}

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
      print(c(sum(sel),range(wd),range(p)))
      nsel = p < 0 | wd < 0
      if (all(sel == nsel)) break
      sel = nsel
    }
  } else {
    p = matrix(0,n,n)
    wd = w
  }
  list(p=p, w=w, wd=wd, sump=sum(p),frac=sum(sel)/length(p))
}


k = expand.grid(k1=seq_circ(n),k2=seq_circ(n))
kern = 1/sqrt(k$k1^2+k$k2^2)
kern[1] = 0

f1 = frac$f1
f2 = frac$f2
#f1[] = 0
#f2[] = sqrt(0.3^2 - (frac$points[,1]-0.5)^2 - (frac$points[,2]-0.5)^2)
#f2[is.nan(f2)] = 0

clear3d()
surface3d(x,x,f1,asp="iso", col="white")
surface3d(x,x,f2,asp="iso", col="white")


x  = seq(0,frac$span[1,1],len=n+1)[-(n+1)]

K = 100
w = f1-f2
d = seq(min(w),mean(w),len=K)
d = mean(w) + (min(w)-mean(w))*10^seq(0,-3,len=K)
ret = lapply(d, fun)
sump = sapply(ret, function(x)x$sump)
plot(d,-sump,log="y")

perc = sapply(ret, function(x)x$frac)
plot(d,perc,log="y")

plot(-sump/(n*n)*pi,perc,log="xy")
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
      x = seq(0,frac$span[1,1],len=n+1),
      y = seq(0,frac$span[2,2],len=n+1))

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
for (i in seq_along(ret)) {
  print(i)
  val = ret[[i]]
  wd = val$wd
  p = -val$p
  points$z = to_n_plus_1(f1+(wd-w)/2)
  write.vtk.tri(
    points,
    triangles,
    point_data = list(p=to_n_plus_1(p)),
    filename = sprintf("vtk_output/frac_A_%04d.vtk",i)
  )
  points$z = to_n_plus_1(f2-(wd-w)/2)
  write.vtk.tri(
    points,
    triangles,
    point_data = list(p=to_n_plus_1(p)),
    filename = sprintf("vtk_output/frac_B_%04d.vtk",i)
  )
}

