begin
    using SparseArrays
    using LinearAlgebra
    using Plots
    using IterativeSolvers
    L=1#m
    L_c=0.1#m
    h=10#w/m^2.c
     k=10#w/m.c
     Bi=h*L_c/k
    d1=spzeros(1,200)
    CG=spzeros(1,200)
    G=spzeros(1,200)
    trecg=spzeros(1,200)
    treg=spzeros(1,200)
    c1=-(L/L_c*(Bi)^(1/2)*sinh(L/L_c*(Bi)^(1/2))+h*L/k*cosh(L/L_c*(Bi)^(1/2)))/(L/L_c*(Bi)^(1/2)*cosh(L/L_c*(Bi)^(1/2))+h*L/k*sinh(L/L_c*(Bi)^(1/2)))
As(x)=c1*sinh.(L/L_c*(Bi)^(1/2)*x)+cosh.(L/L_c*(Bi)^(1/2)*x)
    for j=1:200
        n=1*j+1
        d1[1,j]=L/(n-1)
      c=(d1[1,j]^2*L^2/(L_c)^2)*Bi+2
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
           a[n,n]=c+(2*d1[1,j]*Bi*L/L_c)
           a[n,n-1]=2
           end
    #rhs
        b=spzeros(n,1)
        for i=2:n
           b[i,1]=0-a[i-1,1]
        end
        #Dirichlet B.C.
        b[1,1]=1
    cg=cg!(spzeros(n,1),a,b,maxiter=100)
    CG[1,j]=cg[n,1]
    g=gauss_seidel!(spzeros(n,1),a,b,maxiter=2500)
    G[1,j]=g[n,1]
    end
    scecg=spzeros(1,199)
    sce=spzeros(1,199)
    for z=1:199
        scecg[1,z]=@.abs(CG[1,z+1]-CG[1,z])
        trecg[1,z]=@.abs(CG[1,z]-As(1))
        sce[1,z]=@.abs(G[1,z+1]-G[1,z])
        treg[1,z]=@.abs(G[1,z]-As(1))
    end
    
    end
    plt3=plot(xaxis=:log10,title="Errors",xlabel="\$d\$",legend=:topleft,linewidth=2)
    plot!(plt3,d1[1,1:199],scecg[1,:],label=" Cg successive error",color=:blue)
     scatter!(plt3,d1[1,1:199],scecg[1,:],label=" Cg successive error")  
    plot!(plt3,d1[1,:],trecg[1,:],label=" Cg true error",linestyle=:dash,color=:red)
     scatter!(plt3,d1[1,:],trecg[1,:],label=" Cg true error")
     plot!(plt3,d1[1,1:199],sce[1,:],label=" Gauss-Siedel successive error")
 scatter!(plt3,d1[1,1:199],sce[1,:],label=" Gauss-Siedel successive error")    
 plot!(plt3,d1[1,:],treg[1,:],label=" Gauss-siedel true error",linestyle=:dash,color=:purple)
     scatter!(plt3,d1[1,:],treg[1,:],label=" Gauss-siedel true error")