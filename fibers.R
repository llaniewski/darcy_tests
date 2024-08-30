n = 128
A = array(0,c(n,n,n))

x = seq(0,1,len=n)
p = as.matrix(expand.grid(x=x,y=x,z=x))

sel = FALSE
r1 = 0.02
r2 = 0.5

for (it in 1:100) {
  M = diag(c(20,1,1)) %*% matrix(rnorm(9),3,3)
  M = qr.Q(qr(M))
  C = runif(3)
  pm = p
  pm = cbind(pm,-1) %*% rbind(diag(3),C)
  pm = pm %*% M
  sel = sel | sqrt((pm[,1]/r1)^2+(pm[,2]/r1)^2+(pm[,3]/r2)^2) < 1
}

A[sel] = 1
mean(A)

library(reticulate)
vtk = import("vtk")
numpy_support = import("vtk.util.numpy_support")

im = vtk$vtkImageData()
im$SetDimensions(dim(A)+1L)
vtk_sel = numpy_support$numpy_to_vtk(as.vector(A)*1.0)
vtk_sel$SetName("A")
im$GetCellData()$AddArray(vtk_sel)
writer = vtk$vtkXMLImageDataWriter()
writer$SetFileName("out.vti")
writer$SetInputData(im)
writer$Write()


