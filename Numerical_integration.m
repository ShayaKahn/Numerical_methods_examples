%% Baton Thrown into Air.
clear all
% Parameters.
m1 = 0.1;m2 = 0.1;L = 1;g = 9.81;
% Initial conditions.
u0 = [0; 4; L; 20; -pi/2; 2];
% Equations.
f = @(q) [1 zeros(1,5);0 m1+m2 zeros(1,3 ...
    ) -m2*L*sin(q(5));zeros(1,2) 1 zeros(1,3);zeros( ...
    1,3) m1+m2 0 m2*L*cos(q(5)); zeros(1,4) 1 0;0 -L*sin( ...
    q(5)) 0 L*cos(q(5)) 0 L^2]\[q(2);m2*L*q(6)^2*cos(q(5));q( ...
    4);m2*L*q(6)^2*sin(q(5))-(m1+m2)*g;q(6);-g*L*cos(q(5))];
h = 0.1; % Step size
t = 0:h:4; % time 
% Euler initiation.
u_e= nan(6,length(t));
u_e(:,1)= u0;
% Euler method solution.
for i=1:length(t)-1
    u_e(:,i+1)=u_e(:,i) + h*f(u_e(:,i));
end
% RK initiation.
u_rk= nan(6,length(t));
u_rk(:,1)= u0;
% RK method solution.
for i=1:length(t)-1
    u=u_rk(:,i);
    k1 = h*f(u);
    k2 = h*f(u+k1/2);
    k3= h*f(u+k2/2);
    k4=h*f(u+k3);
    u_rk(:,i+1)=u+1/6*(k1+2*k2+2*k3+k4);
end
figure;
title('Motion of a Thrown Baton, Solved by Euler',Interpreter='latex');
axis([0 22 0 22])
hold on
% Here we plot the position of m1, m2 and the rod between them for both
% methods.
for j = 1:length(t)
   theta = u_e(5,j);
   X = u_e(1,j);
   Y = u_e(3,j);
   xvals = [X X+L*cos(theta)];
   yvals = [Y Y+L*sin(theta)];
   plot(xvals,yvals,xvals(1),yvals(1),'r.',xvals(2),yvals(2),'g.')
end
xlabel('x',Interpreter='latex');ylabel('y',Interpreter='latex')
hold off
figure;
title('Motion of a Thrown Baton, Solved by RK',Interpreter='latex');
axis([0 22 0 22])
hold on
for j = 1:length(t)
   theta = u_rk(5,j);
   X = u_rk(1,j);
   Y = u_rk(3,j);
   xvals = [X X+L*cos(theta)];
   yvals = [Y Y+L*sin(theta)];
   plot(xvals,yvals,xvals(1),yvals(1),'r.',xvals(2),yvals(2),'g.')
end
xlabel('x',Interpreter='latex');ylabel('y',Interpreter='latex')
hold off