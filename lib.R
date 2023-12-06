
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
