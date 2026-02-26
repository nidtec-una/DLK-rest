function [xf, logres, VectResNorm]= lgmres_v2(A,b,mL,lL,tol,itermax,x0)
maxit=itermax;
m=mL;
d=lL;
n=size (A,2);
flag=0;
r=b-A*x0;
normR0 = max(norm(b), eps);
logres(1,:)= norm(r)/normR0;
miteracion(1,1)=m;

%Preallocating for speed
w=zeros(n,m+d);
z=zeros(n,1);
ij=1;

% Residuos normalizados
VectResNorm=r/max(norm(r), eps);

while flag==0
    beta=norm(r);
    v(:,1)=r/beta;
    h=zeros(m+1,m);
    if size(logres,1)==1
        for j=1:m                       %modified gram schmidt--Arnoldi
            w(:,j)=A*v(:,j);
            for i=1:j
                h(i,j)=w(:,j)'*v(:,i);
                w(:,j)=w(:,j)-h(i,j)*v(:,i);
            end
            h(j+1,j)=norm(w(:,j));
            if h(j+1,j)==0
                m=j;
                h2=zeros(m+1,m);
                for k=1:m
                    h2(:,k)=h(:,k);
                end
                h=h2;
            else
                v(:,j+1)=w(:,j)/h(j+1,j);
            end
        end
        g=zeros(m+1,1);
        g(1,1)=beta;
        for j=1:m                       %plane rotations (QR decompostion)
            P=eye(m+1);
            sin=h(j+1,j)/(sqrt(h(j+1,j)^2 + h(j,j)^2));
            cos=h(j,j)/(sqrt(h(j+1,j)^2 + h(j,j)^2));
            P(j,j)=cos;
            P(j+1,j+1)=cos;
            P(j,j+1)=sin;
            P(j+1,j)=-sin;
            h=P*h;
            g=P*g;
        end
        R=zeros(m,m);
        G=zeros(m,1);
        V=zeros(n,m);
        for k=1:m
            G(k)=g(k);
            V(:,k)=v(:,k);
            for i=1:m
                R(k,i)=h(k,i);
            end
        end
        minimizer=R\G;
        Z=V*minimizer;
        xm=x0 + Z;
        r=b-A*xm;
        VectResNorm(:,size(VectResNorm,2)+1)=r/max(norm(r), eps);    %Residuo normalizado
        miteracion(size(miteracion,1)+1,1)=m;
        logres(size(logres,1)+1,:)=norm(r)/normR0;
        
        if logres(size(logres,1)) <tol
            flag=1;
        else
            x0=xm;                        %update and restart
        end
        
        %Calculo de z(k)
        z(:,ij)= Z;
        %ij=ij+1;
    else
        if ij<=lL
            d=ij;
            ij=ij+1;
        end
        s=m+d;
        %Modified gram schmidt--Arnoldi
        for j=1:s           
            if j<=m
                w(:,j)=A*v(:,j);
            else
                w(:,j)=A*z(:,d-(j-m-1));
            end
            for i=1:j
                h(i,j)=w(:,j)'*v(:,i);
                w(:,j)=w(:,j)-h(i,j)*v(:,i);
            end
            h(j+1,j)=norm(w(:,j));
            if h(j+1,j)==0
                s=j;
                h2=zeros(s+1,s);    %VERIFICAR!!!...
                for k=1:s
                    h2(:,k)=h(:,k);
                end
                h=h2;
            else
                v(:,j+1)=w(:,j)/h(j+1,j);
            end
        end
        g=zeros(s+1,1);
        g(1,1)=beta;

        %Plane rotations (QR decompostion)
        for j=1:s                       
            P=eye(s+1);
            sin=h(j+1,j)/(sqrt(h(j+1,j)^2 + h(j,j)^2));
            cos=h(j,j)/(sqrt(h(j+1,j)^2 + h(j,j)^2));
            P(j,j)=cos;
            P(j+1,j+1)=cos;
            P(j,j+1)=sin;
            P(j+1,j)=-sin;
            h=P*h;
            g=P*g;
        end
        R=zeros(s,s);
        G=zeros(s,1);
        V=zeros(n,s);
        for k=1:s
            G(k)=g(k);
            V(:,k)=v(:,k);
            for i=1:s
                R(k,i)=h(k,i);
            end
        end
        for k=m+1:s
            V(:,k)=z(:,d-(k-m-1));
        end
        minimizer=R\G;
        xm=x0+V*minimizer;
        r=b-A*xm;
        VectResNorm(:,size(VectResNorm,2)+1)=r/max(norm(r), eps);  % Residuo Normalizad
        miteracion(size(miteracion,1)+1,1)=m;
        logres(size(logres,1)+1,:)=norm(r)/normR0;
        aux=V*minimizer;
        Z=z;
        if size(z,2)<lL
            z(:,size(z,2)+1)=aux;
        else
            for j=2:lL
                z(:,j-1)=Z(:,j);
            end
            z(:,lL)=aux;
        end
        
        if logres(size(logres,1),1) <tol || size(logres,1)==maxit
            flag=1;
        else
            x0=xm;                        %update and restart
        end
        
        
    end
end 
xf = xm;
