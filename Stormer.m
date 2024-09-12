%% Stormer solution - one planet and sun.

% Parameters.
m=200;
% Initial conditions.
x0=3; 
y0=3;
vx0=0;
vy0=1;

% Time.
T = 100;
h = 0.00001;
t=0:h:T;

% Initiation.
x = zeros(1,length(t));
y = zeros(1,length(t));

% Set initial conditions.
x(1) = x0;
y(1) = y0;

% Step 1.
x(2) = x(1)+h*vx0+(h^2/2)*(-m*x0/(y0^2+x0^2)^(3/2));
y(2) = y(1)+h*vy0+(h^2/2)*(-m*y0/(y0^2+x0^2)^(3/2));

for i=2:length(t)
    D=((x(i))^2+(y(i))^2)^(3/2);
    ax = -m*x(i)/D;
    ay= -m*y(i)/D;
    x(i+1) = 2*x(i)-x(i-1)+h^2*ax;
    y(i+1) = 2*y(i)-y(i-1)+h^2*ay;
end

figure; 
plot(x,y)
%% Euler solution - one planet and sun.
clear all
T = 100;
h = 0.00001;
t=0:h:T;
m=200;
r_c=zeros(length(t),1);
x=zeros(length(t),1);
y=zeros(length(t),1);
v_x=zeros(length(t),1);
v_y=zeros(length(t),1);
x(1)=3;
y(1)=3; 
v_x(1)=0;  
v_y(1)=1; 
for i= 1:length(t)-1 
    r_c(i)=sqrt(x(i)^2+y(i)^2);
    v_x(i+1)=v_x(i)-(m*x(i)*h)/(r_c(i)^3);
    v_y(i+1)=v_y(i)-(m*y(i)*h)/(r_c(i)^3);
    x(i+1)=x(i)+v_x(i+1)*h;
    y(i+1)=y(i)+v_y(i+1)*h;
end
figure;plot(x,y)