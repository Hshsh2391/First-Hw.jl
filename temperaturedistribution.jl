 begin
 using SparseArrays
 using IterativeSolvers
 using Plots
 L=1#m
 L_c=0.1#m
n=11#nodes
 d=L/(n-1)#m
 h=10#w/m^2.c
 k=100#w/m.c
 Bi=h*L_c/k
c=(d^2*L^2/(L_c)^2)*Bi+2
    a=spzeros(n,n)
     for i=2:n-1
        a[i,i]=-c
        a[i,i-1]=a[i,i+1]=1
     end
        #drichlet boundary cind.
        a[1,1]=1;
        for i=2:n
        a[1,i]=a[i,1]=0
        #neumann B.C.
        a[n,n]=c+(2*d*Bi*L/L_c)
        a[n,n-1]=-2
        end
 #rhs
     b=spzeros(n,1)
     for i=2:n
        b[i,1]=0-a[i-1,1]
     end
     #Dirichlet B.C.
     b[1,1]=1
     c1=-(L/L_c*(Bi)^(1/2)*sinh(L/L_c*(Bi)^(1/2))+h*L/k*cosh(L/L_c*(Bi)^(1/2)))/(L/L_c*(Bi)^(1/2)*cosh(L/L_c*(Bi)^(1/2))+h*L/k*sinh(L/L_c*(Bi)^(1/2)))
As(x)=c1*sinh.(L/L_c*(Bi)^(1/2)*x)+cosh.(L/L_c*(Bi)^(1/2)*x)
g=IterativeSolvers.gmres!(b,a,1e-6)
#cgg=cg!(spzeros(n,1), a, b,maxiter=150)
x=LinRange(0,L,n)
plt1=plot(xlabel="\$x^*\$", linewidth=2,legend=:topleft)
plot!(plt1 ,x,g,label="Gauss-siedle solution",color=:blue)
#plot!(plt1,LinRange(0,L,n),cgg,label="Cg solution",color=:green)
 plot!(plt1,x,As(x),label="Analtical solution",color=:red,linestyle=:dash)
end

