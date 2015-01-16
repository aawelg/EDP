
%% 1) Matrices de differentiation

ecart_4=zeros(100,1);
ecart_2=ecart_4;
k=1;

for n=10:1e3:1e5

    h=2*pi/n;
 
    % D_2
    i=[1:(n-1),2:n]; j=[2:n,1:(n-1)];
    D_2=sparse(i,j,[ ones(1,n-1), -ones(1,n-1)] ,n,n);
    D_2(1,n)=-1; D_2(n,1)=1;
    D_2=1/(2*h)*D_2;

    %full(D_2)

    % D_4
    i=[1:(n-1),2:n]; j=[2:n,1:(n-1)];
    S1=sparse(i,j,2/3*[ ones(1,n-1), -ones(1,n-1)],n,n);

    i=[1:(n-2),3:n]; j=[3:n,1:(n-2)];
    S2=sparse(i,j,-1/12*[ ones(1,n-2), -ones(1,n-2)],n,n);

    D_4=(S1+S2);
    D_4(1,n)=-2/3; D_4(n,1)=2/3; D_4(1,n-1)=1/12; D_4(n-1,1)=-1/12;
    D_4(2,n)=1/12; D_4(n,2)=-1/12;
    D_4=1/h*D_4;

    %full(D_4)

    x=-pi+h:h:pi; x=x';

    ux=exp(sin(x));
    dux=cos(x).*exp(sin(x));

    ecart_2(k)=norm(D_2*ux-dux);
    ecart_4(k)=norm(D_4*ux-dux);

    k=k+1;
    
end


% affichage des ecarts 
 n=10:1e3:1e5;
 figure(1);
 hold on;
 plot(log10(n),log10(ecart_2),log10(n),-2*log10(n))
 plot(log10(n),log10(ecart_4),log10(n),-4*log10(n))
 hold off;
 
 
 %%  2) DFT
 
 nbf=4;
 N=52:2:100; err=zeros(nbf,length(N));
 
for j=1:length(N)
    
    n=N(j);
    h=2*pi/n;
    k=[0:n/2-1 0 -n/2+1:-1];
    x=-pi+h:h:pi; x=x';
    
    y=zeros(n,nbf);
    y(:,1)=exp(sin(x)); y(:,2)=abs(sin(x)).^3;  %y(:,5)=exp(-1./(sin(x)).^2);
    y(:,3)=1./(sin(x).^2+1); y(:,4)=sin(x).^10;
    
    I=sqrt(-1);
    ck=fft(y);
    
    du=zeros(n,nbf);
        for l=1:nbf   
            du(:,l)=ifft( I* k'.*ck(:,l));
        end
    
    yu=zeros(n,nbf);
    yu(:,1)=cos(x).*exp(sin(x)); yu(:,2)=3*abs(cos(x).*sin(x).^2);
    yu(:,3)=-2*cos(x).*sin(x)./(sin(x).^2+1).^2;
    yu(:,4)=10*cos(x).*sin(x).^9;
      %yu(:,5)=(2*cos(x)./sin(x).^3).*exp(-1./(sin(x).^2));
    
    dyu=zeros(1,nbf);
        for l=1:nbf   
            dyu(l)=norm(du(:,l)-yu(:,l));
        end    
    err(:,j)=dyu'; 
    
end

% affichage

figure(1)

hold on
subplot(2,2,1)
plot(err(1,:),'-b');
grid 

subplot(2,2,2)
plot(err(2,:),'-r');
grid 

subplot(2,2,3)
plot(err(3,:),'-m');
grid

subplot(2,2,4)
plot(err(4,:),'-g');
grid 
hold off


 %%  3) 

 




