%% Solution to one body motion around the sun - Euler.
% Parameters
clear all
m=1;
f=@(u)[u(2)
    -m*u(1)./((u(1).^2+u(3).^2).^(3/2))
    u(4)
    -m*u(3)./((u(1).^2+u(3).^2).^(3/2))];
u0=[2;0;0;0.5];
T=50*pi;
h=1/1000;
t = 0:h:T; % time 
% Euler initiation.
u_e= nan(4,length(t));
u_e(:,1)= u0;
% Euler method solution.
for i=1:length(t)-1
    u_e(:,i+1)=u_e(:,i) + h*f(u_e(:,i));
end
figure;
plot(t,u_e,LineWidth=1.5)
xlabel('t',Interpreter='latex')
legend('x','$v_x$','y','$v_y$',Interpreter='latex')
figure;
plot(u_e(1,:),u_e(3,:),'-',0,0,'ro')
xlabel('x',Interpreter='latex');ylabel('y',Interpreter='latex')
title('Solution by Euler',Interpreter='latex')
%% One body motion around the sun - Stormer - Verlet.
% Parameters
m=1;
% Set initial conditions
x0=2; 
y0=0;
vx0=0;
vy0=0.5;
% time.
T = 50*pi;
h = 1/1000;
t=0:h:T;
% Initiation.
Lx = zeros(1,length(t));
Ly = zeros(1,length(t));
Lx(1) = x0;
Ly(1) = y0;
% First step outside the loop.
Lx(2) = Lx(1)-(h^2/2)*m*x0/(x0^(3/2));
Ly(2) = Ly(1)+h*vy0;
%  Stormer - Verlet solution.
for i=2:length(t)
    D=((Lx(i))^2+(Ly(i))^2)^(3/2);
    ax = -m*Lx(i)/D;
    ay= -m*Ly(i)/D;
    Lx(i+1) = 2*Lx(i)-Lx(i-1)+h^2*ax;
    Ly(i+1) = 2*Ly(i)-Ly(i-1)+h^2*ay;
end
figure; 
plot(Lx,Ly,'-',0,0,'ro')
xlabel('x',Interpreter='latex');ylabel('y',Interpreter='latex')
title('Solution by Stormer - Verlet',Interpreter='latex')
%% pendulum solution by Euler.
clear all
f=@(u) [u(2) ; -sin(u(1))]; % ODE function
theta0 = 0.1;% initial amplitude
h = 0.001; % step size
t=0:h:3*pi*50; % time span
% initiation for Euler method
u_e= nan(2,length(t));
u_e(:,1)= [theta0;0];
% Euler method solution
for i=1:length(t)-1
    u_e(:,i+1)=u_e(:,i) + h*f(u_e(:,i));
end
figure % Phase space
plot(u_e(1,:),u_e(2,:))
xlabel('\theta','FontSize',20,'Rotation',0)
legend('Euler','Location','eastoutside')
ylabel('$\dot{\theta}$','Interpreter','latex','FontSize',20,'Rotation',0)
title('Phase space',Interpreter='latex')
%% pendulum solution by Velocity - Verlet.
theta0 = 0.1;% initial amplitude
h = 0.001; % step size
t=0:h:3*pi*50; % time span
% initiation for Euler method
u_vv= nan(2,length(t));
u_vv(:,1)= [theta0;0];
% Leap-Frog method solution
for i=1:length(t)-1
    u_vv(1,i+1)=u_vv(1,i) + h*u_vv(2,i)+h^2/2*(-sin(u_vv(1,i)));
    u_vv(2,i+1)=u_vv(2,i) + h/2*(-sin(u_vv(1,i))-sin(u_vv(1,i+1)));
end
figure % Phase space
plot(u_vv(1,:),u_vv(2,:))
xlim([-0.15,0.15]);ylim([-0.15,0.15])
xlabel('\theta','FontSize',20,'Rotation',0)
legend('Velocity Verlet','Location','eastoutside')
ylabel('$\dot{\theta}$','Interpreter','latex','FontSize',20,'Rotation',0)
title('Phase space',Interpreter='latex')