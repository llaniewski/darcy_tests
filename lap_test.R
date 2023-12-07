
o = cumsum(c(0,nrow(face),nrow(edge),nrow(flux)))
N = o[length(o)]

flux_trans = 1 / flux$d


ret = list(
  list( i=o[1]+flux$faceidx, j=o[3]+flux$idx,      x = edge$len[flux$edgeidx]   ),
  list( i=o[2]+flux$edgeidx, j=o[3]+flux$idx,      x = 1 ),
  list( i=o[3]+flux$idx,     j=o[3]+flux$idx,      x = 1 ),
  list( i=o[3]+flux$idx,     j=o[1]+flux$faceidx,  x = flux_trans ),
  list( i=o[3]+flux$idx,     j=o[2]+flux$edgeidx,  x =-flux_trans )
)

ret = equate.lists(ret)
ret = glue.lists(ret)

eqidx = unique(ret$i[is.na(ret$j)])
sel = is.na(ret$i) | ret$i %in% eqidx
ret = lapply(ret, "[", !sel)


plot(ret$j, ret$i, col=colrange(ret$x),pch=16,asp=1,ylim=c(N,0),xlim=c(0,N))
abline(v=o); abline(h=o)

M = sparseMatrix(i=ret$i, j=ret$j, x=ret$x)

A = M[(o[1]+1):o[2],(o[3]+1):o[4]]
B = M[(o[2]+1):o[3],(o[3]+1):o[4]]
C = M[(o[3]+1):o[4],(o[2]+1):o[3]]
D = M[(o[3]+1):o[4],(o[1]+1):o[2]]


G = A %*% D - A %*% C %*% solve(B %*% C) %*% B %*% D



Gd = as.matrix(G)
#Gdi = solve(Gd)
Gde = eigen(Gd)

n = nrow(f)

v = Gde$values
sel = v<max(v)*1e-5
v = 1/sqrt(v)
v[sel] = 0
plot(v)

Gp = Gde$vectors %*% diag(v) %*% t(Gde$vectors)

write.vtk.tri(vertex$point, face, 
              cell_data=list(
                v1 = Gde$vectors[,n-11],
                v2 = Gde$vectors[,n-12],
                v3 = Gde$vectors[,n-13],
                v4 = Gde$vectors[,n-14]
              ),filename="lap.vtk")

write.vtk.tri(vertex$point, face, 
              cell_data=list(
                v1 = Gp[,1],
                v2 = Gp[,2],
                v3 = Gp[,3],
                v4 = Gp[,4]
              ),filename="lap.vtk")

