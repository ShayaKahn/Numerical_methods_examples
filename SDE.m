%% Wiener Process
clear all
clc
T=2;

% parameters
N=10000;
dt=T/N;
t=0:dt:T;
reals=1000;
x_f=nan(1,reals);
x_mid=nan(1,reals);

for i=1:reals
    dw=randn(1,N)*sqrt(dt);
    w=[0,cumsum(dw)];
    x_f(i)=w(end);
    x_mid(i)=w(round(length(w)/2));
    if mod(i,20)==0
        plot(t,w)
        hold on
    end    
end 

% fit distribution to the final time
pd = fitdist(x_f','Normal');
x=-5:0.1:5;
y=pdf(pd,x);
plot(2-y,x,'k-',LineWidth=1.5)
hold on

% fit distribution to the middle time
pd = fitdist(x_mid','Normal');
x=-5:0.1:5;
y=pdf(pd,x);
plot(1-y,x,'k-',LineWidth=1.5)
hold on
syms t
fplot(t,sqrt(t),'b-',LineWidth=2)
xlabel('t', 'Interpreter','latex')
ylabel('x', 'Interpreter','latex')
legend('$$\sqrt{t}$$','Interpreter','latex',Location='northwest')
title(['Weiner process - ',num2str(reals),' realizations'], ...
    'Interpreter','latex')
xlim([0 2])
%% Original Langevin
clear
clc
close all

% parameters
gamma = 1/2;
m = 1;
ksi = 1;
v0 = 10;
t0 = 0;
tf = 20;
n = 1000;
dt = (tf-t0)/n;
t = t0:dt:tf;

% initialization
v = zeros(size(t));
v(1) = v0;
dW = sqrt(dt)*randn(1,n);

% solution using Euler Maruyama for one iteration
for i=1:n
    v(i+1) = v(i) - dt*gamma/m*v(i) + ksi/m*dW(i); 
end
figure; hold on
plot(t,v)
plot(xlim,[0 0],'--k')
xlabel('t','Interpreter','latex');
ylabel('v(t)','Interpreter','latex')
%% solution using Euler Maruyama for 5000 iterations
real = 5000; 
v = zeros(real,n+1);
v(:,1) = v0;
for i=1:real
    dW = sqrt(dt)*randn(1,n);
    for j=1:n
        v(i,j+1) = v(i,j) - dt*gamma/m*v(i,j) + ksi/m*dW(j); 
    end
end
figure;
plot(t,mean(v))
ylabel('$\overline{v}$','Interpreter','latex');
xlabel('t','Interpreter','latex')
title('mean velocity vs time','Interpreter','latex')

%% velocity distributions for different times
figure; hold on
legendinfo = {};
for j=1:7
    idxj = 40*(j-1)+1;
    [counts, centers] = hist(v(:,idxj),31) ;
    legendinfo{j} = ['t = ', num2str(t(idxj))];
    plot(centers, counts/trapz(centers, counts))
end
xlim([-6,12])
xlabel('v','Interpreter','latex');
ylabel('P(v)','Interpreter','latex');
legend(legendinfo, 'Location', 'northwest');
title('Distribution of velocity in different times','Interpreter','latex')
%% Euler Maruyama
clear all
clc

% t vector
N=1000;
t0=0;
tf=5;
dt=tf/N;
t=t0:dt:tf;

% initial condition
x0=0;

% number of realizations
reals=10000;

% parameters
sigma=0.1;

% initialization
x=nan(reals, N+1);
x(:,1) = x0;

% Create dw
dw = sqrt(dt)*randn(reals, N+1);

% Force
F=@(x)x-x.^3;

% Solution
for i=1:reals
    for j=1:N
        F_x=F(x(i,j));
        x(i,j+1) = x(i,j) + F_x*dt + sigma*dw(i,j);
    end
end

% Plot of one iteration
figure;
plot(t, x(end,:));
xlabel('t',Interpreter='latex');
ylabel('X(t)',Interpreter='latex');
title('$F=x-x^3$',Interpreter='latex');
hold off

% Compute P(x)
figure;
numBins = 100; 
[counts, edges] = histcounts(x(:,end), numBins, ...
    'Normalization', 'probability');
% Plot the probability distribution
bar(edges(1:end-1), counts, 'hist');
xlabel('x',Interpreter='latex');
ylabel('P(x)',Interpreter='latex');
hold off
%% Focker Plank

% x vector
M=1000;
x0=-4;
xf=4;
dx=(2*xf)/M;
x=x0:dx:xf;

% t vector
N=1000;
t0=0;
tf=5;
dt=tf/N;
t=t0:dt:tf;

% parameters
sigma=0.1;

% initialization
P=nan(M+1,N+1);
P(:,1)=zeros(1,M+1);
P(M/2+1,1)=1;

% define the F
F=@(x)x-x.^3;

% diagonals
d=(-dt/(2*dx^2)*sigma^2).*ones(1,M+1);
d_up=dt/(4*dx^2)*sigma^2-dt/(4*dx)*F(x(2:end));
d_down=dt/(4*dx^2)*sigma^2+dt/(4*dx)*F(x(1:end-1));

% matrix S
S=diag(d)+diag(d_down,-1)+diag(d_up,1);
S(1,2)=2.*S(1,2);
S(end,end-1)=2.*S(end,end-1);

% Crank Nicolson matrix
CN=sparse((eye(M+1)-S)\(eye(M+1)+S));

% solution
for i=1:N
    P(:,i+1)=CN*P(:,i);
end

% plots
figure;
plot(x,P(:,end),"Color",'blue','LineWidth',2)
xlabel('x',Interpreter='latex')
ylabel('P(x)',Interpreter='latex')
figure;
imagesc(t,x,P)
xlabel('t',Interpreter='latex')
xlabel('x',Interpreter='latex')
title('P(x,t)',Interpreter='latex')
colorbar
set(gca,'ColorScale','log')