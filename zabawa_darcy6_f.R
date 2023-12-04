library(rfracture)
library(rgl)
library(RcppArmadillo)

Sys.setenv("PKG_LIBS"="-lsuperlu")
Rcpp::sourceCpp("spsolve.cpp")

write.vtk.tri = function(points, triangles, point_data=list(), cell_data=list(), filename) {
  f = file(filename,"w")
  wt = function(x) write.table(x, file=f, row.names = FALSE, col.names = FALSE)
  wl = function(...) writeLines(c(...), con=f)
  wl(
    "# vtk DataFile Version 5.1",
    "vtk output",
    "ASCII",
    "DATASET POLYDATA",
    paste("POINTS", nrow(points),"float")
  )
  wt(points[,1:3])
  wl("")
  wl(
    paste("POLYGONS",nrow(triangles)+1,nrow(triangles)*3),
    "OFFSETS vtktypeint64"
  )
  wt((seq_len(nrow(triangles)+1)-1)*3)
  wl(
    "CONNECTIVITY vtktypeint64"
  )
  wt(triangles[,1:3]-1)
  wl("")
  if (length(cell_data) > 0) {
    wl(
      paste("CELL_DATA",nrow(triangles))
    )
    for (n in names(cell_data)) {
      v = cell_data[[n]]
      if (is.null(nrow(v))) v = data.frame(x=v)
      if (nrow(v) != nrow(triangles)) stop("wrong length of cell_data")
      if (ncol(v) == 1) {
        wl(
          paste("SCALARS", n, "double"),
          "LOOKUP_TABLE default"
        )
        wt(v)
      } else if (ncol(v) == 3) {
        wl(
          paste("VECTORS", n, "double")
        )
        wt(v)
      } else {
        stop("wrong ncol to cell_data")
      }
      wl("")
    }
  }
  if (length(point_data) > 0) {
    wl(
      paste("POINT_DATA",nrow(points))
    )
    for (n in names(point_data)) {
      v = point_data[[n]]
      if (is.null(nrow(v))) v = data.frame(x=v)
      if (nrow(v) != nrow(points)) stop("wrong length of point_data")
      if (ncol(v) == 1) {
        wl(
          paste("SCALARS", n, "double"),
          "LOOKUP_TABLE default"
        )
        wt(v)
      } else if (ncol(v) == 3) {
        wl(
          paste("VECTORS", n, "double")
        )
        wt(v)
      } else {
        stop("wrong ncol to point_data")
      }
      wl("")
    }
  }
  close(f)
}

vecprod = function(x,y) cbind(x[,2]*y[,3]-x[,3]*y[,2],x[,3]*y[,1]-x[,1]*y[,3],x[,1]*y[,2]-x[,2]*y[,1])
vecnorm = function(x) sqrt(rowSums(x^2))
inv = function(M) cbind(M[,4],-M[,2],-M[,3],M[,1])/(M[,1]*M[,4]-M[,2]*M[,3])
mul = function(A,B) cbind(
  A[,1]*B[,1]+A[,3]*B[,2],
  A[,2]*B[,1]+A[,4]*B[,2],
  A[,1]*B[,3]+A[,3]*B[,4],
  A[,2]*B[,3]+A[,4]*B[,4]
)
grad = function(p1,p2,p3) {
  V2 = p1 - p2
  V1 = p3 - p2
  M11 = rowSums(V1*V1)
  M12 = rowSums(V1*V2)
  M22 = rowSums(V2*V2)
  DET = M11*M22 - M12*M12
  (-M12*V1 + M11*V2) / DET
}
colrange = function(v, colors=c("black","red","yellow")) {
  col = v
  col = (col-min(col))/(max(col)-min(col))
  col = colorRamp(colors=colors)(col)
  col = rgb(col,max=255)
  col
}

equate.lists = function(x) {
  lapply(x, function(x) {
    l = max(sapply(x, length))
    lapply(x, function(x) if (length(x) == 1) rep(x,l) else x)
  })
}

glue.lists = function(x) {
  nm = names(x[[1]])
  nm = setNames(as.list(nm),nm)
  lapply(nm, function(nm) do.call(c, lapply(x, "[[", nm)))
}
myappend = function(x,y) append(x,list(y))




frac = fracture_geom(refine=1, power.iso=function(f) 0.00001*exp(-(f-1)^2),gap=0.1)
frac = cut(frac)
#write.stl(frac,"~/vtk_test/test.vtk")

frac = fracture_geom(refine=6, power.iso=function(f) 0.00001*exp(-(f-1)^2),gap=0.1)
frac = cut(frac)

  frac$points$fm = 0

obj = frac
if (FALSE) {
a = obj$points$x
b = obj$points$y
obj$points$x = a*2
obj$points$y = cos(pi*b)
obj$points$z = -sin(pi*b)
} else {
  obj$points$z = 0  
}

p = data.frame(x1=obj$points$x,
               x2=obj$points$y,
               x3=obj$points$z)


f = data.frame(i1=obj$triangles[,1],
               i2=obj$triangles[,2],
               i3=obj$triangles[,3])

plot3d(p, aspect = "iso")
i = as.vector(t(as.matrix(f)))
triangles3d(p[i,], col=3)


v1 = p[f$i2,] - p[f$i1,]
v2 = p[f$i3,] - p[f$i1,]
n = vecprod(v1,v2)
f$area = vecnorm(n)/2
n = n/vecnorm(n)
g1 = grad(p[f$i1,],p[f$i2,],p[f$i3,])
g2 = grad(p[f$i2,],p[f$i3,],p[f$i1,])
g3 = grad(p[f$i3,],p[f$i1,],p[f$i2,])


tmp = cbind(rowSums(v1*v1),rowSums(v1*v2),rowSums(v2*v1),rowSums(v2*v2))
tmp = inv(tmp)
tmp = mul(tmp,cbind(rowSums(v1*v1),rowSums(v2*v2),0,0))[,1:2]
tmp = (tmp[,1]*v1 + tmp[,1]*v2)/2
fp = p[f$i1,]+tmp

fp = (p[f$i1,]+p[f$i2,]+p[f$i3,])/3

e = data.frame(i1=c(f[,1],f[,2],f[,3]),
               i2=c(f[,2],f[,3],f[,1]),
               t1=rep(seq_len(nrow(f)),times=3))
v3 = p[e$i2,] - p[e$i1,]
e$len = vecnorm(v3)
en =  vecprod(v3,n[e$t1,]) / e$len
sel = e$i1 > e$i2
e[sel,1:2] = e[sel,2:1]
sel = order(e$i1,e$i2,e$t1)
en = en[sel,]
e = e[sel,]

t1n = en
t2n = en
t2n[] = NA
sel = diff(e$i1) == 0 & diff(e$i2) == 0
e$t2[c(sel,FALSE)] = e$t1[c(FALSE,sel)]
t2n[c(sel,FALSE),] = t1n[c(FALSE,sel)]
sel = !c(FALSE,sel)
e = e[sel,]
t1n = t1n[sel,]
t2n = t2n[sel,]

b = data.frame(i1=obj$edge[,1],i2=obj$edge[,2],border=obj$border)
sel = b$i1 > b$i2
b[sel,1:2] = b[sel,2:1]
eb = merge(b,cbind(e,eid=seq_len(nrow(e))))
e$border[eb$eid] = eb$border

e[is.na(e$t2),]

background_pressure = function(x) { rho*(as.matrix(x) %*% g) }

leakoff_coef = 0.01

e$zones = "interior"
e$zones[is.na(e$t2)] = "wall"
sel = p[,1] == min(p[,1]) & p[,2] >= 0.6
sel = p[,1] == min(p[,1]) & p[,2] >= 0.35 & p[,2] <= 0.65
e$zones[sel[e$i1] & sel[e$i2] & (!is.na(e$border)) & e$border] = "inlet"
sel = p[,1] == max(p[,1]) & p[,2] >= 0.6 
#e$zones[sel[e$i1] & sel[e$i2] & (!is.na(e$border)) & e$border] = "outlet"


d1 = fp[e$t1,] - p[e$i1,]
d2 = fp[e$t2,] - p[e$i1,]

e$d1 = -rowSums(t1n*d1)/vecnorm(t1n)
e$d2 = -rowSums(t2n*d2)/vecnorm(t1n)
e$d2[is.na(e$d2)] = 1e-6

ep = (p[e$i1,] + p[e$i2,])/2


g = c(0,-3,0)

tab = data.frame(phi = rep(0,nrow(f)))
tab$h = 0.075
rho = 1
rho_s = 2.5
tab$phi_rcp = 0.5
tab$nu = 0.01
tab$d = 0.015
tab$h = 0.075

#sel = fp[,1] > 0.30 & fp[,1] < 0.7 & fp[,2] > 0.30 & fp[,2] < 0.7
#tab$phi[sel] = 0.4

#tab$phi = 0.3

D = 0.0005
dt = 0.05

tab$rho_m = rho*(1-tab$phi) + rho_s*tab$phi


o = cumsum(c(0,rep(nrow(f),2),rep(nrow(e),6),nrow(f)))
N = o[length(o)]
e$idx = seq_len(nrow(e))
f$idx = seq_len(nrow(f))

only_wall = !any(e$zones %in% c("inlet","outlet"))

iter = -1

start_iter = iter+1
X11()
for (iter in start_iter + 1:10000-1) {

  
  tab$F = tab$phi/pmax(tab$phi_rcp-tab$phi,0.00001)*1
  #tab$F = 0
  tab$G = (18*tab$nu)/(tab$d^2)*tab$phi*rho
  tab$H = (12*tab$nu)/(tab$h^2)#*(1-tab$phi)*rho
  
  
perA = cbind(tab$H+tab$G,-tab$G,-tab$G,tab$F+tab$G)
perP = cbind((1-tab$phi)*rho, 0, 0, tab$phi*rho_s)
perW = cbind(1,1,1/rho, 1/rho_s)
perM = mul(perP, mul(inv(perA), perP))
sel = tab$phi < 1e-5
perM[sel,] = mul(perP, mul(cbind(1/perA[,1],0,0,0), perP))[sel,]
perC = mul(perM,perW)
perB = mul(perW,perC)

# M' = perB[,1] * g + perB[,3] * dp
# V' = perB[,2] * g + perB[,4] * dp


g1 = perB[e$t1,2] * as.vector(t1n %*% g)
g2 = perB[e$t2,2] * as.vector(t2n %*% g)
trans1 = -perB[e$t1,4] / e$d1
trans2 = -perB[e$t2,4] / e$d2
trans1m = -perB[e$t1,3] / e$d1
trans2m = -perB[e$t2,3] / e$d2

trans1D = -D / e$d1
trans2D = -D / e$d2


RHS = rep(0,N)

RHS[o[2]+f$idx] = tab$rho_m*f$area
RHS[o[5]+e$idx] = perB[e$t1,2] * as.vector(t1n %*% g)
RHS[o[6]+e$idx] = perB[e$t2,2] * as.vector(t2n %*% g)
RHS[o[7]+e$idx] = perB[e$t1,1] * as.vector(t1n %*% g)
RHS[o[8]+e$idx] = perB[e$t2,1] * as.vector(t2n %*% g)

sel = !is.na(e$t2)
ret = list(
  list( i=o[1]+e$t1, j=o[4]+e$idx, x=e$len*dt   ),
  list( i=o[1]+e$t2, j=o[5]+e$idx, x=e$len*dt ),
  list( i=o[2]+e$t1, j=o[6]+e$idx, x=e$len*dt ),
  list( i=o[2]+e$t2, j=o[7]+e$idx, x=e$len*dt ),
  list( i=o[2]+f$idx, j=o[2]+f$idx, x=f$area ),
  list( i=o[3]+e$idx, j=o[4]+e$idx, x=1 ),
  list( i=o[3]+e$idx, j=o[5]+e$idx, x=1 ),
  list( i=o[4]+e$idx, j=o[6]+e$idx, x=1 ),
  list( i=o[4]+e$idx, j=o[7]+e$idx, x=1 ),
  list( i=o[5]+e$idx, j=o[4]+e$idx, x=1 ),
  list( i=o[5]+e$idx, j=o[1]+e$t1,  x= trans1 ),
  list( i=o[5]+e$idx, j=o[3]+e$idx, x=-trans1 ),
  list( i=o[6]+e$idx, j=o[5]+e$idx, x=1 ),
  list( i=o[6]+e$idx, j=o[1]+e$t2,  x= trans2 ),
  list( i=o[6]+e$idx, j=o[3]+e$idx, x=-trans2 ),
  list( i=o[7]+e$idx, j=o[6]+e$idx, x=1 ),
  list( i=o[7]+e$idx, j=o[1]+e$t1,  x= trans1m ),
  list( i=o[7]+e$idx, j=o[3]+e$idx, x=-trans1m ),
  list( i=o[8]+e$idx, j=o[7]+e$idx, x=1 ),
  list( i=o[8]+e$idx, j=o[1]+e$t2,  x= trans2m ),
  list( i=o[8]+e$idx, j=o[3]+e$idx, x=-trans2m ),
  list( i=o[7]+e$idx, j=o[2]+e$t1,  x= trans1D ),
  list( i=o[7]+e$idx, j=o[8]+e$idx, x=-trans1D ),
  list( i=o[8]+e$idx, j=o[2]+e$t2,  x= trans2D ),
  list( i=o[8]+e$idx, j=o[8]+e$idx, x=-trans2D ),
  list( i=o[9]+f$idx, j=o[1]+f$idx, x=-leakoff_coef),
  list( i=o[9]+f$idx, j=o[9]+f$idx, x=1),
  list( i=o[1]+f$idx, j=o[9]+f$idx, x=1*f$area*dt ),
  list( i=o[2]+f$idx, j=o[9]+f$idx, x=1*f$area*dt )
)

RHS[o[9]+f$idx] = - leakoff_coef * background_pressure(fp)

ret = equate.lists(ret)
ret = glue.lists(ret)

sel = !is.na(ret$i)
ret = lapply(ret, "[", sel)

eqidx = e$idx[is.na(e$t2)]
eqidx = c(o[6] + eqidx, o[8] + eqidx)
sel = rep(TRUE, N)
sel[eqidx] = FALSE
sel = sel[ret$i]
ret = lapply(ret, "[", sel)

ret_add = list()

eidx = which(e$zones == "inlet")
if (length(eidx) > 0) {
dofidx = o[3] + eidx
eqidx  = o[6] + eidx
val    = background_pressure(ep[eidx,]) + 2
ret_add = myappend(ret_add, list( i=eqidx, j=dofidx, x=1 ))
RHS[eqidx] = val
dofidx = o[7] + eidx
eqidx  = o[8] + eidx
val    = 0
ret_add = myappend(ret_add, list( i=eqidx, j=dofidx, x=1))
ret_add = myappend(ret_add, list( i=eqidx, j=o[5]+eidx, x=-1.15))
RHS[eqidx] = val
}
eidx = which(e$zones == "outlet")
if (length(eidx) > 0) {
dofidx = o[3] + eidx
eqidx  = o[6] + eidx
val    = background_pressure(ep[eidx,])
ret_add = myappend(ret_add, list( i=eqidx, j=dofidx, x=1))
RHS[eqidx] = val
dofidx = o[8] + eidx
eqidx  = o[8] + eidx
val    = 0
ret_add = myappend(ret_add, list( i=eqidx, j=dofidx, x=1))
ret_add = myappend(ret_add, list( i=eqidx, j=o[2]+e$t1[eidx], x=-1))
RHS[eqidx] = val
}
eidx = which(e$zones == "wall")
if (length(eidx) > 0) {
dofidx = o[5] + eidx
eqidx  = o[6] + eidx
val    = 0
ret_add = myappend(ret_add, list( i=eqidx, j=dofidx, x=1))
RHS[eqidx] = val
dofidx = o[7] + eidx
eqidx  = o[8] + eidx
val    = 0
ret_add = myappend(ret_add, list( i=eqidx, j=dofidx, x=1))
RHS[eqidx] = val
}
if (FALSE) {
  plot(M$j,N-M$i,asp=1,xlim=c(0,N),ylim=c(0,N),col=ifelse(M$x == 1, 2,3),pch=16,cex=0.5)
  abline(h=N-o); abline(v=o)
}

if (only_wall) {
  ret_add = myappend(ret_add, list( i=o[1]+f$idx, j=N+1, x=f$area ))
  ret_add = myappend(ret_add, list( i=N+1, j=o[1]+f$idx, x=f$area ))
  RHS = c(RHS,0)
}


ret = glue.lists(c(equate.lists(ret_add),list(ret)))

#M = sparseMatrix(i=ret$i,j=ret$j,x=ret$x,dims=rep(length(RHS),2))
#x = solve(M,RHS)

x = arma_spsolve(rbind(ret$i,ret$j)-1, ret$x, RHS)

#x = .solve.dgC.lu(M,RHS)

#print(range(x[o[8]+e$idx]))



if (FALSE) {
  plot(fp[,1:2],asp=1,col=colrange(x[o[1]+f$idx]))
  segments(p[e$i1,1],p[e$i1,2],p[e$i2,1],p[e$i2,2],col=colrange(x[o[3]+e$idx]))
  
  plot(fp[,1:2],asp=1,col=colrange(x[o[2]+f$idx]))
  s = x[o[4]+e$idx]
  segments(ep[,1],ep[,2],ep[,1]+t1n[,1]*s,ep[,2]+t1n[,2]*s,colrange(abs(s)))
}


#print(range(x[o[2]+f$idx]))

nrhom = x[o[2]+f$idx]
drhom = (tab$rho_m - nrhom)/dt
#print(range(drhom))
tab$rho_m = nrhom
tab$rho_m = pmax(0,tab$rho_m)
#plot(fp[,1:2],asp=1,col=colrange(tab$rho_m))

tab$phi = (tab$rho_m - rho)/(rho_s - rho)

if (iter %% 50 == 0) {
  k = nrow(f)
  S = data.frame(
    i=c(e$t1,e$t1+k,e$t1+k*2,e$t2,e$t2+k,e$t2+k*2),
    j=c(e$idx,e$idx,e$idx,e$idx,e$idx,e$idx),
    x=c(t1n[,1],t1n[,2],t1n[,3],t2n[,1],t2n[,2],t2n[,3])
  )
  S = S[!is.na(S$i),]
  S = sparseMatrix(i=S$i,j=S$j,x=S$x,dims=c(nrow(f)*3,nrow(e)))
  dp = S %*% (x[o[3]+e$idx] * e$len)
  dp = as.vector(dp) / f$area
  dp = matrix(dp,k,3)
  #range(S %*% (e$len * x[o[4]+e$idx]))
  Qf = perC[,1]*outer(rep(1,k),g) - perC[,3]*dp
  Qs = perC[,2]*outer(rep(1,k),g) - perC[,4]*dp
  
  plot(fp[,1:2],asp=1,col=colrange(tab$phi),pch=16)
  #plot3d(p, aspect = "iso")
  #i = as.vector(t(as.matrix(f[,1:3])))
  #triangles3d(p[i,], col=rep(colrange(tab$phi),each=3))
  
  write.vtk.tri(p, f, 
    cell_data=list(
      phi=tab$phi,
      pres=x[o[1]+f$idx],
      pres_diff=x[o[1]+f$idx] - background_pressure(fp),
      feakoff=x[o[9]+f$idx],
      drhom=drhom,
      dp=dp,
      Qf=Qf,
      Qs=Qs
    ),filename=sprintf("vtk_output/test_%08d.vtk",iter))
}
toprint = numeric(0)
for (k in c("inlet","outlet","wall")) {
  sel = e$zones == k
  toprint[[paste0(k,"_flux")]] = sum(x[o[5]+e$idx[sel]]*e$len[sel])
}
toprint["total_phi"] = sum(x[o[2]+f$idx]*f$area)
toprint["leakoff"] = sum(x[o[9]+f$idx]*f$area)
toprint["mean_pressure_diff"] = sum((x[o[1]+f$idx] - background_pressure(fp))*f$area) / sum(f$area)
print(toprint)
}

for (i in dev.list()) dev.off(i)
