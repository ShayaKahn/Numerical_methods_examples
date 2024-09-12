%% PDE
clear all
clc
%% Crank Nicholson Example
%% Parameters
D=0.1;
dx=0.01;
dt=0.005;
alpha=D*dt/dx^2;
tmax=2;
L=1;
t=0:dt:tmax;
x=dx:dx:L-dx;
N=length(x);
T=length(t);
%% Set initial conditions
T0=zeros(N,1);
T0(1:round(N/2))=1;
T_D=nan(N,T);
T_D(:,1)=T0;
T_N=nan(N,T);
T_N(:,1)=T0;
T_P=nan(N,T);
T_P(:,1)=T0;
%% Create matrices
% Main diagonal for S
d=-2*ones(N,1);

% diagonals for dirichlet and periodic Boundary conditions a
d_DP_up=ones(N,1);
d_DP_down=ones(N,1);

% diagonals for Neumann Boundary conditions
d_N_up=[0;2;ones(N-2,1)];
d_N_down=[ones(N-2,1);2;0];

% Specify the S matrices
S_D=spdiags([d d_DP_up d_DP_down],[0 1 -1],N,N);
S_N=spdiags([d d_N_up d_N_down],[0 1 -1],N,N);
S_P=spdiags([d d_DP_up d_DP_down],[0 1 -1],N,N);

% Set periodic Boundary conditions
S_P(end,1)=1;
S_P(1,end)=1;

% Specify the CN matrices
CN_D=sparse((eye(N)-0.5*alpha*S_D)\(eye(N)+0.5*alpha*S_D)); 
CN_N=sparse((eye(N)-0.5*alpha*S_N)\(eye(N)+0.5*alpha*S_N));
CN_P=sparse((eye(N)-0.5*alpha*S_P)\(eye(N)+0.5*alpha*S_P));
%% Solution
for i = 2:length(t)
    T_D(:,i)=CN_D*T_D(:,i-1);
    T_N(:,i)=CN_N*T_N(:,i-1);
    T_P(:,i)=CN_P*T_P(:,i-1);
end
%% Plots
close all
T_cell={T_D,T_N,T_P};
title_cell={'dirichlet BC','Neumann BC','periodic BC'};

% Plots
for i=1:numel(T_cell)
    figure;
    %imagesc(0:dx:1,t,T_cell{i}')
    imagesc(x,t,T_cell{i}')
    axis('xy')
    xlabel('x','FontSize',20);ylabel('t','FontSize',20);
    title(title_cell{i})
    colorbar
    ax=axes('Position',[.6,.7,.2,.2],'Units','normalized');
    plot(x,T_cell{i}(:,end),'-',LineWidth=2)
    hold on
    plot(x,T_cell{i}(:,1),'-',LineWidth=2)
    xlabel('x','Color','w');ylabel('T','Color','w')
    legend('f','i')
    hold off
end
%% Operator splitting example

% parameters
tf=3;
l=10;
dx=0.1;
dt=0.01;
t=0:dt:tf;
x=0:dx:l;
L=length(x);
D=0.1;% diffusion coefficient
v=1; % wave velocity
alpha=D*dt/dx^2;
r=v*dt/dx;
u0=exp(-4*(x-7).^2); % initial condition

% CN diagonals
d=-2*ones(L,1);
d_up=ones(L,1);
d_down=ones(L,1);

% CN matrix
S=spdiags([d d_up d_down],[0 1 -1],L,L);
CN=sparse((eye(L)-0.5*alpha*S)\(eye(L)+0.5*alpha*S));

% Euler diagonals
d=(1-r)*ones(L,1);
d_up=r*ones(L,1);

% Euler matrix
Euler=spdiags([d d_up],[0 1],L,L);

% Solution
u=nan(L,length(t));
u(:,1)=u0;
for i=1:length(t)-1
    u_temp=CN*u(:,i);
    u(:,i+1)=Euler*u_temp;
end
% Plots
figure;
p=plot(x,u0);
axis([0 l 0 1])
xlabel('x');ylabel('u')
T=text(l/2,.5,'t=','HorizontalAlignment','center');
for i=1:length(t)
    p.YData=u(:,i);
    pause(dt);
    if ~mod(t(i),.1)
    T.String=['t=' num2str(t(i))];
    end
end
hold off
figure;
imagesc(x,t,u');axis('xy')
xlabel('x','FontSize',20);
ylabel('t','FontSize',20);