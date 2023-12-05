library(rfracture)
library(rgl)
library(RcppArmadillo)
library(tibble)
rm(list=ls())

Sys.setenv("PKG_LIBS"="-lsuperlu")
Rcpp::sourceCpp("spsolve.cpp")

write.vtk.tri = function(points, triangles, point_data=list(), cell_data=list(), filename) {
  face = file(filename,"w")
  wt = function(x) write.table(x, file=face, row.names = FALSE, col.names = FALSE)
  wl = function(...) writeLines(c(...), con=face)
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
    paste("POLYGONS",nrow(triangles)+1L,nrow(triangles)*3L),
    "OFFSETS vtktypeint64"
  )
  wt((seq_len(nrow(triangles)+1L)-1L)*3L)
  wl(
    "CONNECTIVITY vtktypeint64"
  )
  wt(triangles[,1:3]-1L)
  wl("")
  if (length(cell_data) > 0) {
    wl(
      paste("CELL_DATA",nrow(triangles))
    )
    for (field_name in names(cell_data)) {
      v = cell_data[[field_name]]
      if (is.null(nrow(v))) v = data.frame(x=v)
      if (nrow(v) != nrow(triangles)) stop("wrong length of cell_data")
      if (ncol(v) == 1) {
        wl(
          paste("SCALARS", field_name, "double"),
          "LOOKUP_TABLE default"
        )
        wt(v)
      } else if (ncol(v) == 3) {
        wl(
          paste("VECTORS", field_name, "double")
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
    for (field_name in names(point_data)) {
      v = point_data[[field_name]]
      if (is.null(nrow(v))) v = data.frame(x=v)
      if (nrow(v) != nrow(points)) stop("wrong length of point_data")
      if (ncol(v) == 1) {
        wl(
          paste("SCALARS", field_name, "double"),
          "LOOKUP_TABLE default"
        )
        wt(v)
      } else if (ncol(v) == 3) {
        wl(
          paste("VECTORS", field_name, "double")
        )
        wt(v)
      } else {
        stop("wrong ncol to point_data")
      }
      wl("")
    }
  }
  close(face)
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

frac = fracture_geom(refine=6, power.iso=function(face) 0.00001*exp(-(face-1)^2),gap=0.1)
frac = cut(frac)

frac$points$fm = 0


p = cbind(frac$points$x,frac$points$y,0)
p = cbind(p[,1],p[,2]-0.5,p[,3])
triangles = as.vector(t(frac$triangles))
p = p[triangles,]

sel = p[,1] >= 0.5
sel = rep(colSums(matrix(sel,3)) == 3,each=3)
np = p[sel,]
np = p
np = cbind(0.5,np[,2],np[,1]-0.5)
p = rbind(p,np)

a = pi/4
p = p %*% matrix(c(1,0,0,0,cos(a),sin(a),0,-sin(a),cos(a)), 3,3)

i = factor(sprintf("%3f %3f %3f",p[,1],p[,2],p[,3]))
i = as.integer(i)
sel = !duplicated(i)
p = p[sel,]
p[i[sel],] = p
f = matrix(i, ncol=3,byrow=TRUE)

obj = list(
  points=data.frame(x=p[,1],y=p[,2],z=p[,3]),
  triangles=data.frame(i1=f[,1],i2=f[,2],i3=f[,3])
)


vertex = tibble(point=cbind(
    obj$points$x,
    obj$points$y,
    obj$points$z
))

face = tibble(i1=obj$triangles[,1],
               i2=obj$triangles[,2],
               i3=obj$triangles[,3])
face$idx = seq_len(nrow(face))

v1 = vertex$point[face$i2,] - vertex$point[face$i1,]
v2 = vertex$point[face$i3,] - vertex$point[face$i1,]
face$normal = vecprod(v1,v2)
face$area = vecnorm(face$normal)/2
face$normal = face$normal/vecnorm(face$normal)
face$center = (vertex$point[face$i1,]+vertex$point[face$i2,]+vertex$point[face$i3,])/3

flux = tibble(i1=c(face$i1,face$i2,face$i3),
              i2=c(face$i2,face$i3,face$i1),
              faceidx=rep(face$idx,times=3),
              edgeidx=as.integer(NA))
flux$idx = seq_len(nrow(flux))
v3 = vertex$point[flux$i2,] - vertex$point[flux$i1,]
flux$len = vecnorm(v3)
flux$normal =  vecprod(v3,face$normal[flux$faceidx,]) / flux$len
rm(v3)

edge = tibble(i1=flux$i1, i2=flux$i2, fluxidx=flux$idx, len=flux$len)
flux$i1 = NULL
flux$i2 = NULL
flux$len = NULL
sel = edge$i1 > edge$i2
edge[sel,1:2] = edge[sel,2:1]
edge = edge[order(edge$i1, edge$i2),]
sel = !duplicated(edge[,c("i1","i2")])
edge$idx = cumsum(sel)
flux$edgeidx[edge$fluxidx] = edge$idx
edge = edge[sel,]
edge$fluxidx = NULL

edge$center = (vertex$point[edge$i1,] + vertex$point[edge$i2,])/2

print(table(table(flux$edgeidx)))
print(table(table(flux$faceidx)))

background_pressure = function(x) { rho*(as.matrix(x) %*% g) }
leakoff_coef = 0.01

edge$order = table(flux$edgeidx)
edge$border = edge$order == 1

flux$zone = as.factor("interior")

sel = vertex$point[,1] == min(vertex$point[,1]) & sqrt(vertex$point[,2]^2+vertex$point[,3]^2) <= 0.15
sel = sel[edge$i1] & sel[edge$i2]
flux = rbind(flux,tibble(faceidx=NA, edgeidx=edge$idx[sel], idx=NA, normal=NA, zone="inlet"))


flux = flux[order(flux$edgeidx),]
flux$idx = seq_len(nrow(flux))
sel = !duplicated(flux$edgeidx)
i = flux$idx[sel][flux$edgeidx]
i = flux$idx-i
flux = flux[order(i,flux$edgeidx),]
plot(flux$edgeidx)

flux$idx = seq_len(nrow(flux))

v1 = face$center[flux$faceidx,] - edge$center[flux$edgeidx,]
flux$d = -rowSums(flux$normal*v1)/vecnorm(flux$normal)
rm(v1)

g = c(0,-3,0)

phi = rep(0,nrow(face))
h = 0.075
rho = 1
rho_s = 2.5
phi_rcp = 0.5
nu = 0.01
d = 0.015
h = 0.075

D = 0.0005
dt = 0.05

rho_m = rho*(1-phi) + rho_s*phi

o = cumsum(c(0,rep(nrow(face),2),rep(nrow(edge),2),rep(nrow(flux),2),nrow(face)))
N = o[length(o)]

only_wall = !any(flux$zone %in% c("inlet","outlet"))

iter = -1

start_iter = iter+1
X11()
for (iter in start_iter + 1:100-1) {

  F = phi/pmax(phi_rcp-phi,0.00001)*1
  #F = 0
  G = (18*nu)/(d^2)*phi*rho
  H = (12*nu)/(h^2)#*(1-phi)*rho
  
  
  perA = cbind(H+G,-G,-G,F+G)
  perP = cbind((1-phi)*rho, 0, 0, phi*rho_s)
  perW = cbind(1,1,1/rho, 1/rho_s)
  perM = mul(perP, mul(inv(perA), perP))
  sel = phi < 1e-5
  perM[sel,] = mul(perP, mul(cbind(1/perA[,1],0,0,0), perP))[sel,]
  perC = mul(perM,perW)
  perB = mul(perW,perC)
  
  # M' = perB[,1] * g + perB[,3] * dp
  # V' = perB[,2] * g + perB[,4] * dp
  
  flux_g_vol = perB[flux$faceidx,2] * as.vector(flux$normal %*% g)
  flux_g_mass = perB[flux$faceidx,1] * as.vector(flux$normal %*% g)
  flux_trans_vol = -perB[flux$faceidx,4] / flux$d
  flux_trans_mas = -perB[flux$faceidx,3] / flux$d
  flux_trans_mas_diff = -D / flux$d

  RHS = rep(0,N)
  
  RHS[o[2]+face$idx] = rho_m * face$area
  RHS[o[5]+flux$idx] = flux_g_vol
  RHS[o[6]+flux$idx] = flux_g_mass

  ret = list(
    list( i=o[1]+flux$faceidx, j=o[5]+flux$idx,      x=edge$len[flux$edgeidx]*dt   ),
    list( i=o[1]+face$idx,     j=o[7]+face$idx,      x=1*face$area*dt ),
    list( i=o[2]+face$idx,     j=o[2]+face$idx,      x=face$area ),
    list( i=o[2]+flux$faceidx, j=o[6]+flux$idx,      x=edge$len[flux$edgeidx]*dt ),
    list( i=o[2]+face$idx,     j=o[7]+face$idx,      x=1*face$area*dt ),
    list( i=o[3]+flux$edgeidx, j=o[5]+flux$idx,      x=1 ),
    list( i=o[4]+flux$edgeidx, j=o[6]+flux$idx,      x=1 ),
    list( i=o[5]+flux$idx,     j=o[1]+flux$faceidx,  x= flux_trans_vol ),
    list( i=o[5]+flux$idx,     j=o[3]+flux$edgeidx,  x=-flux_trans_vol ),
    list( i=o[5]+flux$idx,     j=o[5]+flux$idx,      x=1 ),
    list( i=o[6]+flux$idx,     j=o[1]+flux$faceidx,  x= flux_trans_mas ),
    list( i=o[6]+flux$idx,     j=o[2]+flux$faceidx,  x= flux_trans_mas_diff ),
    list( i=o[6]+flux$idx,     j=o[3]+flux$edgeidx,  x=-flux_trans_mas ),
    list( i=o[6]+flux$idx,     j=o[4]+flux$edgeidx,  x=-flux_trans_mas_diff ),
    list( i=o[6]+flux$idx,     j=o[6]+flux$idx,      x=1 ),
    list( i=o[7]+face$idx,     j=o[1]+face$idx,      x=-leakoff_coef),
    list( i=o[7]+face$idx,     j=o[7]+face$idx,      x=1)
  )

  RHS[o[7]+face$idx] = - leakoff_coef * background_pressure(face$center)
  
  ret = equate.lists(ret)
  ret = glue.lists(ret)
  
  eqidx = unique(ret$i[is.na(ret$j)])
  sel = is.na(ret$i) | ret$i %in% eqidx
  ret = lapply(ret, "[", !sel)

  ret_add = list()
  
  bc = flux[flux$zone == "inlet",]
  if (nrow(bc) > 0) {
    eqidx  = o[5] + bc$idx
    ret_add = myappend(ret_add, list( i=eqidx, j=o[3] + bc$edgeidx, x=1 ))
    RHS[eqidx] = background_pressure(edge$center[bc$edgeidx,]) + 2
    eqidx  = o[6] + bc$idx
    ret_add = myappend(ret_add, list( i=eqidx, j=o[6] + bc$idx, x=1 ))
    ret_add = myappend(ret_add, list( i=eqidx, j=o[5] + bc$idx, x=-1.15 ))
    RHS[eqidx] = 0
  }
  bc = flux[flux$zone == "outlet",]
  if (nrow(bc) > 0) {
    stop("can't do outlet yet")
    eqidx  = o[5] + bc$idx
    ret_add = myappend(ret_add, list( i=eqidx, j=o[3] + bc$edgeidx, x=1 ))
    RHS[eqidx] = background_pressure(edge_center[eidx,]) + 2
    eqidx  = o[6] + bc$idx
    ret_add = myappend(ret_add, list( i=eqidx, j=o[6] + bc$idx, x=1 ))
    ret_add = myappend(ret_add, list( i=eqidx, j=o[5] + bc$idx, x=-1.15 ))
    RHS[eqidx] = 0
  }
  
  if (only_wall) {
    ret_add = myappend(ret_add, list( i=o[1]+face$idx, j=N+1, x=face$area ))
    ret_add = myappend(ret_add, list( i=N+1, j=o[1]+face$idx, x=face$area ))
    RHS = c(RHS,0)
  }

  ret = glue.lists(c(equate.lists(ret_add),list(ret)))
  if (FALSE) {
    plot(ret$j,N-ret$i,asp=1,xlim=c(0,N),ylim=c(0,N),col=ifelse(ret$x == 1, 2,3),pch=16,cex=0.5)
    abline(h=N-o); abline(v=o)
  }
  
  #M = sparseMatrix(i=ret$i,j=ret$j,x=ret$x,dims=rep(length(RHS),2))
  #x = solve(M,RHS)
  
  x = arma_spsolve(rbind(ret$i,ret$j)-1, ret$x, RHS)
  
  #x = .solve.dgC.lu(M,RHS)
  
  #print(range(x[o[8]+edge$idx]))
  
  
  
  if (FALSE) {
    plot(face$center[,1:2],asp=1,col=colrange(x[o[1]+face$idx]))
    segments(vertex[edge$i1,1],vertex[edge$i1,2],vertex[edge$i2,1],vertex[edge$i2,2],col=colrange(x[o[3]+edge$idx]))
    
    plot(face$center[,1:2],asp=1,col=colrange(x[o[2]+face$idx]))
    s = x[o[4]+edge$idx]
    segments(edge_center[,1],edge_center[,2],edge_center[,1]+t1n[,1]*s,edge_center[,2]+t1n[,2]*s,colrange(abs(s)))
  }
  
  
  #print(range(x[o[2]+face$idx]))
  
  n_rho_m = x[o[2]+face$idx]
  d_rho_m = (rho_m - n_rho_m)/dt
  rho_m   = n_rho_m
  rho_m = pmin(pmax(rho_m,rho),rho_s)
  
  #plot(face$center[,1:2],asp=1,col=colrange(rho_m))
  
  phi = (rho_m - rho)/(rho_s - rho)
  
  if (iter %% 50 == 0) {
    k = nrow(face)
    S = data.frame(
      i=c(flux$faceidx,flux$faceidx+k,flux$faceidx+k*2),
      j=rep(flux$idx,times=3),
      x=c(flux$normal[,1],flux$normal[,2],flux$normal[,3])
    )
    S = S[!is.na(S$i),]
    S = sparseMatrix(i=S$i,j=S$j,x=S$x,dims=c(nrow(face)*3,nrow(flux)))
    dp = S %*% (x[o[3]+flux$edgeidx] * edge$len[flux$edgeidx])
    dp = as.vector(dp) / face$area
    dp = matrix(dp,k,3)
    #range(S %*% (edge$len * x[o[4]+edge$idx]))
    Qf = perC[,1]*outer(rep(1,k),g) - perC[,3]*dp
    Qs = perC[,2]*outer(rep(1,k),g) - perC[,4]*dp
    
    plot(face$center[,1:2],asp=1,col=colrange(phi),pch=16)
    #plot3d(vertex, aspect = "iso")
    #i = as.vector(t(as.matrix(face[,1:3])))
    #triangles3d(vertex[i,], col=rep(colrange(phi),each=3))
    
    write.vtk.tri(vertex$point, face, 
                  cell_data=list(
                    phi=phi,
                    pres=x[o[1]+face$idx],
                    pres_diff=x[o[1]+face$idx] - background_pressure(face$center),
                    feakoff=x[o[7]+face$idx],
                    drhom=d_rho_m,
                    dp=dp,
                    Qf=Qf,
                    Qs=Qs
                  ),filename=sprintf("vtk_output/darcy7b_%08d.vtk",iter))
  }
  toprint = numeric(0)
  for (k in c("inlet","outlet","wall")) {
    sel = flux$zone == k
    toprint[[paste0(k,"_flux")]] = sum(x[o[5]+flux$idx[sel]]*edge$len[flux$edgeidx[sel]])
  }
  toprint["total_phi"] = sum(x[o[2]+face$idx]*face$area)
  toprint["leakoff"] = sum(x[o[7]+face$idx]*face$area)
  toprint["mean_pressure_diff"] = sum((x[o[1]+face$idx] - background_pressure(face$center))*face$area) / sum(face$area)
  toprint["clamped_mass"] = sum((rho_m - n_rho_m)*face$area)
  print(toprint)
}

for (i in dev.list()) dev.off(i)
