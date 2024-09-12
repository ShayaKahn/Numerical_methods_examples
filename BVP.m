%% Shooting method - Example 1 - estimate the initial guess
% variables 
Equations = @(x,y) [y(2); 3/2*y(1)^2];
time_span= [0 1];
y1_0 = 4; % Value of y at x=0.
y1_1 = 1; % Value of y at x=1.
z = -100:0; % Initial guesses for y'(0).
final = zeros(size(z));

% solve for all the guesses.
for i=1:length(z)
[x, y] = ode45(Equations , time_span, [y1_0 z(i)]');
final(i) = y(end,1);
end
figure; hold on
plot(z, final)
plot(xlim,[1 1],'--k')
title('buandary value for each guess', Interpreter='latex')
ylabel('y(1)', Interpreter='latex')
xlabel('y''(0)', Interpreter='latex')
legend('y(1)','1', Interpreter='latex')
%% Shooting method - Example 1 - Apply shooting method
figure; hold on
legendinfo = {};
% function that calculates the error of our guess from the actual value
% at the end of the integration process.
Error = @(guess)F(guess, y1_0, y1_1, Equations, time_span);
% Try two guesses and plot the integration results.
for guess = [-40 0]
y2_0 = fzero(Error, guess);
[x,y]= ode45(Equations, time_span, [y1_0 y2_0]');
plot(x,y(:,1))
legendinfo = [legendinfo , ['y''(0)=', num2str(y2_0)]];
end
scatter(time_span, [y1_0 y1_1],'r')
title('right solutions', Interpreter='latex');
xlabel('x', Interpreter='latex'); ylabel('y', Interpreter='latex')
legend(legendinfo)
grid
%% Example 2 - Two coupled nonlinear BVP, Shooting two variables

% Define equations and boundray conditions.
Equations = @(t,y)[y(3);y(4);-y(4)^2-y(1)*y(2)+2+exp(t)*(exp(t)+t^2)...
    ;-y(3)+y(2)+2*t];
y1_0=0;
y2_0=0;
y1_1=1;
y2_1=exp(1);
time_span = [0 1];

% Two initial guesses.
s1 = 1;
s2 = 1;
S = [s1; s2];

% Define h for numerical differentiation and tolerance for NR.
h=0.001;
Tol=10^-5;

% The jaccobian and the function f.
f=@(S)[F1(S(1), S(2)); F2(S(1), S(2))];
J=@(S)[(F1(S(1)+h, S(2))-F1(S(1), S(2)))/h, (F1(S(1), S(2)+h)- ...
    F1(S(1), S(2)))/h;
    (F2(S(1)+h, S(2))-F2(S(1), S(2)))/h, (F2(S(1), S(2)+h)- ...
    F2(S(1), S(2)))/h];

% Apply NR and solve for the correct initial guesses by ode45.
[result, ~] = NR(f, J, S, Tol);
[t,y]= ode45(Equations, time_span, [y1_0 ;y2_0; result(1); result(2)]);
figure;
plot(t,y)
xlabel('t',Interpreter='latex')
hold on
scatter([0 1],[0 1],'r')
scatter([0 1],[0 exp(1)],'r')
legend('$y_1$', '$y_2$', '$y_3$', '$y_4$', Interpreter='latex', ...
    Location='best')
%% Example 3 - schrodinger equation.
clear all
clc
figure; hold on
Tol=0.0001;
legendinfo = {};
E_vector=linspace(0,2,100); % The initial guesses for the Energy.
En_vector=zeros();
for i=1:length(E_vector)
    En = fzero(@schroidinger,E_vector(i));
    En_vector(i)=En;
end
% Find the unique elements in En_vector within tolerance.
unique_En=uniquetol(En_vector,Tol); 
for i=1:length(unique_En)
    % Plot the Eginstates.
    [x,y]= ode45(@(x,y)[y(2);-2*unique_En(i)*y(1)],[0 10],[0;1]);
    plot(x,y(:,1))
    % Plot the Energies as legend.
    legendinfo = [legendinfo , ['E_',num2str(i),'=',...
        num2str(unique_En(i))]];
end    
scatter([0 10],[0 0],'r')
title('Egin states and Energies', Interpreter='latex');
xlabel('x', Interpreter='latex'); ylabel('\psi(x)')
legend(legendinfo)
grid
%% Collocation method Example.
h=.01;
theta0=pi/4;
tf=7;
t=h:h:tf-h;
n=length(t);
% coupled equations:
f=@(theta)[(theta(2)+theta0-2*theta(1))/h^2+sin(theta(1))
 (theta(3:n)+theta(1:n-2)-2*theta(2:n-1))/h^2+sin(theta(2:n-1))
 (theta(n-1)+0-2*theta(n))/h^2+sin(theta(n))];
% Jacobian:
J=@(theta) diag((-2/h^2+cos(theta)).*ones(n,1),0)...
 +diag(1/h^2*ones(n-1,1),1)+diag(1/h^2*ones(n-1,1),-1);
% initial guess:
theta=theta0*cos(t'*(3*pi/2)/tf);
% solve using NR:
theta=NRs(f,J,theta,1e-5);
theta_dot=([theta; 0]-[theta0;theta])/h;
figure;
p=polarplot(0,1,'o',[0 0],[0 1]);
ax=gca;
ax.ThetaAxisUnits='radians';
ax.ThetaZeroLocation = 'bottom';
T=text(0,0.1,'t=0','HorizontalAlignment','center');
for i=1:n
    thetai=theta(i);
    p(1).ThetaData=thetai;
    p(2).ThetaData=[0 thetai];
    T.String=['t=' num2str(t(i))];
    pause(h/2)
end
clf
%figure;
plot([0 t tf],[theta0 ;theta ;0])
hold on
plot([0 t],theta_dot)
plot(0,theta0,'bo',tf,0,'bo')
xlabel('t','FontSize',20);%ylabel('\theta ','FontSize',20,'Rotation',0)
title(['h=' num2str(h) ', \theta_0=' num2str(theta0)])
legend({'$\theta$','$\dot{\theta}$'},'Interpreter','latex','Location', ...
    'north','FontSize',20)
%% Functions.
function Error = F(z, y1_0, y1_1, Equations, time_span)
[~, y] = ode45(Equations, time_span, [y1_0 z]');
Error = y(end,1)-y1_1;
end
% The function of E.
function final = schroidinger(E)
    [~, y] = ode45(@(x,y)[y(2);-2*E*y(1)],[0 10],[0; 1]);
    final = y(end,1);
end
function final = F1(y3_0, y4_0)
Equations = @(t,y)[y(3);y(4);-y(4)^2-y(1)*y(2)+2+exp(t)*(exp(t)+t^2)...
    ;-y(3)+y(2)+2*t];
time_span = [0 1];
y1_0=0;
y2_0=0;
y1_1=1;
[~, y] = ode45(Equations, time_span, [y1_0; y2_0; y3_0; y4_0]);
final = y(end,1)-y1_1;
end
function final = F2(y3_0, y4_0)
Equations = @(t,y)[y(3);y(4);-y(4)^2-y(1)*y(2)+2+exp(t)*(exp(t)+t^2);...
    -y(3)+y(2)+2*t];
time_span = [0 1];
y1_0=0;
y2_0=0;
y2_1=exp(1);
[~, y] = ode45(Equations, time_span, [y1_0; y2_0; y3_0; y4_0]);
final = y(end,2)-y2_1;
end
function [x0,n] = NR(f,J,x0,h)
    n=0;
    while norm(f(x0))>h
        x0 = x0-J(x0)\f(x0);
        n=n+1;
    end
end
function [x0,n] = NRs(f,J,x0,h)
    n=0;
    while norm(f(x0))>h
        x0=x0-sparse(J(x0))\f(x0);
        n=n+1;
    end
end