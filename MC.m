%% Moment of inertia - Torus.

clear all
clc
% Define the box.
rho=@(x,y,z) 1;
x_high=4;
x_low=1;
y_high=4;
y_low=-3;
z_high=1;
z_low=-1;
V=(x_high-x_low)*(y_high-y_low)*(z_high-z_low);
N=1000;
n=1000;
Ix=nan(1,N);
f=0;
for i=1:N
    for j=1:n
        % Generate random point inside the box.
        x=x_low+(x_high-x_low)*rand;
        y=y_low+(y_high-y_low)*rand;
        z=z_low+(z_high-z_low)*rand;
         if z^2+(sqrt(x^2+y^2)-3)^2<=1
             % Calculate dI around x for the point and add
             % to the previous calculations.
             f=f+(y^2+z^2)*rho(x,y,z); 
         end
    end
    % Calculate the moment of inertia around x.
    Ix(i)=V*f/(i*n);
end
% Plots.
plot(1:n:n*N,Ix,'.-')
hold on
xlabel('N','FontSize',20)
ylabel('I_x','FontSize',20)

%% Particles on sphere - Metropolis.

clear all
clc
% Number of paticles.
m=50;
% Temprature.
T=0.01;
% Energy function.
J=.1/m;
E=@(X) J*sum(1./pdist(X)); % X is x y z coordinates.
% Initial state.
phi = 2*pi*rand(m,1);
theta = pi*rand(m,1);
R=[theta phi]; % R is theta phi coordinates (R=1)
RtoX = @(R)[cos(R(:,1)).*sin(R(:,2)),sin(R(:,1)).*sin(R(:,2)),cos(R(:,2))];
% # of moves:
N=1e7;
En= nan(N,1); % Enrergy n
Pn=nan(N,1); % prob n
En(1)=E(RtoX(R));
Pn(1)=exp(-En(1)/T);
for i=1:N
    % Choose particle randomly:
    j = randi(m);
    % create an alternative state -> move:
    Rtemp=R;
    Rtemp(j,1)=mod(Rtemp(j,1)+0.1*rand,pi);
    Rtemp(j,2)=mod(Rtemp(j,2)+0.1*rand,2*pi);
    % Calculate the energy of the alternative state.
    Etemp=E(RtoX(Rtemp));
    % Apply metropolis algorithm:.
    if Etemp<En(i) || rand<exp(-(Etemp-En(i))/T)
        En(i+1)=Etemp;
        Pn(i+1)=exp(-Etemp/T);
        R=Rtemp;
    else
        En(i+1)=En(i);
        Pn(i+1)=Pn(i);
    end
end
% Plots.
subplot(1,2,1)
x=0:.05:pi+.05;
y=0:.05:2*pi;
[th,p]=meshgrid(x,y);hold on
x=cos(th).*sin(p);
y=sin(th).*sin(p);
z=cos(p);
mesh(x,y,z);
r=RtoX(R);
scatter3(r(:,1),r(:,2),r(:,3),30,'filled')
title('final configuration')
view(20,30)
subplot(1,2,2)
plot(En,'.')
xlabel('n','FontSize',20)
ylabel('E_n','FontSize',20,'Rotation',0,'HorizontalAlignment','right');
% axis([0 N 0 J*m^2])
Emean=sum(En.*Pn)/sum(Pn);
emean=mean(En);
title(['\langleE\rangle= ' num2str(Emean)])