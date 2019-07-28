clear all
clc
disp('VIBRATION PROJECT')
m=input('Enter the number of masses: ');
    for i=1:m
    b(i)= input(['Enter value of mass ' num2str(i) ': ']);
    M=[diag(b)];
    end
    disp('The mass matrix is: ');
    disp(M);
n = input('Enter number of springs: ');
for c=1:n
k(c) = input(['Enter Stiffness of spring ' num2str(c) ': ']);
end
for i=1:1:m
    for j=1:1:m
        if (i==j&&i<n)
            a(i,j) = k(i)+k(i+1);
        elseif(i==j&&i==n)
            a(i,j) = k(n);
        elseif(j==i+1)
            a(i,j) = -k(j);
        elseif(i==j+1)
            a(i,j) = -k(i);
        end
    end
end
disp('The stiffness matrix is: ');
disp(a);
    [X,D]=eig(a,M);
     for i=1:m
    freq(i)=sqrt(D(i,i));
    modes(:,i)=X(:,i)/X(1,i);
     end
disp('The natural frequencies are: ');
disp(freq);
disp('The mode shapes are: ');
disp(modes);
plot(modes,'LineWidth',2);
grid on
title('Eigen Vector');
Legend=cell(m,1);
 for iter=1:m
   Legend{iter}=strcat('Mode shapes', num2str(iter));
 end
 legend(Legend);
for i=1:m
   u(:,i)=modes(:,i);
end
disp('The modal matrix is: ');
disp(u);
% free vibration response %
for i=1:m
    x(i)=input(['Enter initial displacement of mass ' num2str(i) ': ']);
    v(i)=input(['Enter initial velocity of mass ' num2str(i) ': ']);
end
syms J phi
A=u;
B=transpose(x);
C=linsolve(A,B);
for i=1:m
D(:,i)=[u(:,i)*transpose(-freq(i))];
end
E=transpose(v);
G=linsolve(D,E);
for i=1:m
    J(i)=sqrt(C(i)^2+G(i)^2);
    phi(i)=acos(C(i)/J(i));
end
syms t
uu=transpose(u);
for i=1:m
   one(:,i)=J(i)*u(:,i)*cos((freq(i)*t)+phi(i));
end
final=sum(one,2);
disp('The free vibration response of system is: ');
disp(vpa(final));
startf=input('Enter the start of time interval: ');
endingf=input('Enter the end of time interval: ');
rangef=input('Enter the time gap: ');
aa=startf:rangef:endingf;
gg1=subs(final,t,aa);
plot(aa,gg1,'Linewidth',2);
grid on
title('Free Vibration Response');
xlabel('time');
ylabel('displacement');
Legend=cell(m,1);
 for iter=1:m
   Legend{iter}=strcat('Mass', num2str(iter));
 end
 legend(Legend);
syms t
for i=1:m
    mg(i)=transpose(u(:,i))*M*u(:,i);
    kg(i)=transpose(u(:,i))*a*u(:,i);
    U(:,i)=u(:,i)/(mg(i)^0.5);
end
% forced viration response %
num=input('Enter the number of mass on which force is acting: ');
L=input('Enter the value of external force as given i.e A*cos(w*t): ');
o=input('Enter the value of frequency of force: ');
mat=zeros(1,m);
f=subs(mat,mat(num),L);
for i=1:m
    if i<num
        f(i)=0;
    end
    if i>num
        f(i)=0;
    end
    end
F=U'*transpose(f);
disp('The generalised mass matrix is: ');
disp(mg);
disp('The generalised stiffness matrix is: ');
disp(kg);
disp('The orthonormal mode matrix is: ');
disp(U);
disp('The generalised force matrix is: ');
disp(vpa(F));
for i=1:m
    q(i)=((F(i)/kg(i))*cos(o*t))/(1-((o/freq(i))^2));
end
for i=1:m
  d(:,i)=U(:,i)*q(i);
end
S=sum(d,2);
disp('The forced vibration response of the system is: ');
SS=vpa(S);
disp(SS);
start=input('Enter the start of time interval: ');
ending=input('Enter the end of time interval: ');
a=start:0.049:ending;
gg=subs(SS,t,a);
plot(a,gg,'Linewidth',2);
grid on
title('Forced Vibration Response');
xlabel('time');
ylabel('displacement');
Legend=cell(m,1);
 for iter=1:m
   Legend{iter}=strcat('Mass', num2str(iter));
 end
 legend(Legend);