%% Gillespie example
clear variables
% Parameters.
nA=150;
nX=10;
nY=10;
k0=0.01;
k1=0.1;
k2=1;
t=0;
tf=15;
% Array initiation
r=[k0*nA*nX,k1*nX*nY,k2*nY];
nA_vector=[nA];
nX_vector=[nX];
nY_vector=[nY];
time=[0];
% Gillespie
while t<tf
    R=sum(r);
    num=rand;
    if num<=r(1)/R
        nX=nX+1;
        nA=nA-1;
    elseif (r(1)/R<num) && (num<=(r(2)/R+r(1)/R))
        nY=nY+1;
        nX=nX-1;
    elseif (r(1)/R+r(2)/R<num) && (num<=1)
        nY=nY-1;
    end
    t=t+exprnd(1/R); % Progress time with exponential distribution
    time(end+1)=t;
    % Updates
    r=[k0*nA*nX,k1*nX*nY,k2*nY];
    nX_vector(end+1)=nX;
    nY_vector(end+1)=nY;
    nA_vector(end+1)=nA;
end
plot(time,nX_vector,'r-',time,nY_vector,'g-',time,nA_vector,'b-')
xlabel('x');ylabel('Abundances')
legend('Rabbits','Foxes','Grass')
title('Gillespie example')
