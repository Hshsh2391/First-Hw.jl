begin
    c1=-(L/L_c*(Bi)^(1/2)*sinh(L/L_c*(Bi)^(1/2))+h*L/k*cosh(L/L_c*(Bi)^(1/2)))/(L/L_c*(Bi)^(1/2)*cosh(L/L_c*(Bi)^(1/2))+h*L/k*sinh(L/L_c*(Bi)^(1/2)))
As(x)=c1*sinh.(L/L_c*(Bi)^(1/2)*x)+cosh.(L/L_c*(Bi)^(1/2)*x)
    using SparseArrays
    using LinearAlgebra
    using Plots
    using IterativeSolvers
    L16=(1)#m
    L_c16=Float16(0.1)#m
    h16=Float16(1)#w/m^2.c
     k16=Float16(10)#w/m.c
    Bi16=Float16(h16*L_c16/k16)
    d16=spzeros(1,100)
    G=spzeros(1,100)
    for j=1:100
        n=1*j+1
        d16[1,j]=Float16(L/(n-1))
        c16=Float16((d16[1,j]^2*L16^2/(L_c16)^2)*Bi16+2)
        a=spzeros(n,n)
         for i=2:n-1
            a[i,i]=-c16
            a[i,i-1]=a[i,i+1]=1
         end
            #drichlet boundary cind.
            a[1,1]=1;
            for i=2:n
            a[1,i]=a[i,1]=0
            #neumann B.C.
            a[n,n]=c16+(2*d16[1,j]*Bi16*L16/L_c16)
            a[n,n-1]=-2
            end
     #rhs
         b=spzeros(n,1)
         for i=2:n
            b[i,1]=0-a[i-1,1]
         end
         #Dirichlet B.C.
         b[1,1]=1
    g=gauss_seidel!(zeros(n,1), a, b,maxiter=20000);
    G[1,j]=g[n,1]
    end
    Ger=spzeros(1,100)
    for i=1:100
        Ger[1,i]=@.abs(G[1,i]-As(1));
    end
end
    plt4=plot(xaxis=:log10,title="Error",xlabel="\$d\$",legend=:topright)
     scatter!(plt4,d16[1,1:100],Ger[1,:],label=" Half precision  error")
    begin
     L_c32=Float32(0.1)#m
    h32=Float32(1)#w/m^2.c
     k32=Float32(10)#w/m.c
    Bi32=Float32(h32*L_c32/k32)
    d32=spzeros(1,100)
    G32=spzeros(1,100)
    for j=1:100
        n=1*j+1
        d32[1,j]=Float32(L/(n-1))
        c32=Float32((d32[1,j]^2*L^2/(L_c32)^2)*Bi32+2)
        a=spzeros(n,n)
         for i=2:n-1
            a[i,i]=-c32
            a[i,i-1]=a[i,i+1]=1
         end
            #drichlet boundary cind.
            a[1,1]=1;
            for i=2:n
            a[1,i]=a[i,1]=0
            #neumann B.C.
            a[n,n]=c32+(2*d32[1,j]*Bi32*L/L_c32)
            a[n,n-1]=-2
            end
     #rhs
         b=spzeros(n,1)
         for i=2:n
            b[i,1]=0-a[i-1,1]
         end
         #Dirichlet B.C.
         b[1,1]=1

    g32=gauss_seidel!(zeros(n,1), a, b,maxiter=20000);
    G32[1,j]=g32[n,1]
    end
    Ger32=spzeros(1,100)
    for i=1:100
        Ger32[1,i]=@.abs(G[1,i]-As(1));
    end
    end
     scatter!(plt4,d32[1,1:100],Ger32[1,:],label=" single precision error")

     begin
        L_c64=Float64(0.1)#m
       h64=Float64(1)#w/m^2.c
        k64=Float64(10)#w/m.c
       Bi64=Float64(h64*L_c64/k64)
       d64=spzeros(1,100)
       G64=spzeros(1,100)
       for j=1:100
           n=1*j+1
           d64[1,j]=Float64(L/(n-1))
           c64=Float64((d64[1,j]^2*L^2/(L_c64)^2)*Bi64+2)
           a=spzeros(n,n)
            for i=2:n-1
               a[i,i]=-c64
               a[i,i-1]=a[i,i+1]=1
            end
               #drichlet boundary cind.
               a[1,1]=1;
               for i=2:n
               a[1,i]=a[i,1]=0
               #neumann B.C.
               a[n,n]=c64+(2*d64[1,j]*Bi64*L/L_c64)
               a[n,n-1]=-2
               end
        #rhs
            b=spzeros(n,1)
            for i=2:n
               b[i,1]=0-a[i-1,1]
            end
            #Dirichlet B.C.
            b[1,1]=1
       g64=gauss_seidel!(zeros(n,1), a, b,maxiter=20000);
       G64[1,j]=g64[n,1]
       end
       Ger64=spzeros(1,100)
       for i=1:100
           Ger64[1,i]=@.abs(G64[1,i]-As(1));
       end
   
       end
       scatter!(plt4,d64[1,1:100],Ger64[1,:],label=" double precision error")

       begin
        L_cBig=BigFloat(0.1)#m
       hBig=BigFloat(1)#w/m^2.c
        kBig=BigFloat(10)#w/m.c
       BiBig=BigFloat(hBig*L_cBig/kBig)
       dBig=spzeros(1,100)
       GBig=spzeros(1,100)
       for j=1:100
           n=1*j+1
           dBig[1,j]=BigFloat(L/(n-1))
           cBig=BigFloat((dBig[1,j]^2*L^2/(L_cBig)^2)*BiBig+2)
           a=spzeros(n,n)
            for i=2:n-1
               a[i,i]=-cBig
               a[i,i-1]=a[i,i+1]=1
            end
               #drichlet boundary cind.
               a[1,1]=1;
               for i=2:n
               a[1,i]=a[i,1]=0
               #neumann B.C.
               a[n,n]=cBig+(2*dBig[1,j]*BiBig*L/L_cBig)
               a[n,n-1]=-2
               end
        #rhs
            b=spzeros(n,1)
            for i=2:n
               b[i,1]=0-a[i-1,1]
            end
            #Dirichlet B.C.
            b[1,1]=1
       gBig=gauss_seidel!(zeros(n,1), a, b,maxiter=20000);
       GBig[1,j]=gBig[n,1]
       end
       GerBig=spzeros(1,100)
       for i=1:100
           GerBig[1,i]=@.abs(GBig[1,i]-As(1))
       end
       end
       
       scatter!(plt4,dBig[1,1:100],GerBig[1,:],label=" arbitrary precision error")