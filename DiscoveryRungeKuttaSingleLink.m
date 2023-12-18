%% Research code by Agus Hasan

clear;
clc;

%% time horizon
tf  = 15;
dt  = 0.001;
t   = dt:dt:tf;

%% number of variables and coefficients
n = 4;
r = 60;

%% system description
A = eye(n);
C = eye(n);

%% noise
R = 0;

%% state initialization
x        = zeros(n,1);
xbar     = x;
xhat     = x;
y        = x;
thetabar = zeros(r,1);
thetahat = zeros(r,1);
 
%% true parameters
Ks      = 1.61;         % tahu
Jh      = 0.0021;       % tahu
M       = 0.403;        % tahu
g       = -0.98;        % tahu
h       = 0.06;         % tahu
Km      = 0.00767;
Kg      = 70;
Jl      = 0.0059;       % tahu
Rm      = 2.6;          % tahu

alpha   = Ks/Jh;
beta    = (Km^2*Kg^2)/(Rm*Jh);
gamma   = Ks/Jl;
delta   = M*g*h/Jl;
omega   = (Km*Kg)/(Rm*Jh);

%% initial control inputs
u     = 0;

%% for plotting
uArray          = [];
xArray          = [];
xbarArray       = [];
xhatArray       = [];
yArray          = [];
thetabarArray   = [];
thetahatArray   = [];

%% Initialization for estimator

% lambdav = 0.999;
% lambdat = 0.99999;
% Rx      = 1*eye(n);
% Rt      = 1*eye(n);
% Px      = 100000*eye(n);
% Pt      = 10*eye(r);
% Gamma   = 1*zeros(n,r);
% 
% Pplus       = 1000000*eye(n);
% QF          = 1*eye(n);
% RF          = 1*eye(n);
% a           = 0.999;
% UpsilonPlus = 0*zeros(n,r);
% S           = 10000*eye(r);
% lambda      = 0.99999;

lambdav = 0.999;
lambdat = 0.999;
Rx      = 0.1*eye(n);
Rt      = 1*eye(n);
Px      = 100*eye(n);
Pt      = 10*eye(r);
Gamma   = 0*ones(n,r);

Pplus       = 100000000*eye(n);
QF          = 0.01*eye(n);
RF          = 1*eye(n);
a           = 0.999;
UpsilonPlus = 0*zeros(n,r);
S           = 10000000*eye(r);
lambda      = 0.99999;

%% simulation
for i=1:(tf/dt)

    u = cos(0.5*i*dt);    
%    u = 0.4*cos(0.1*i*dt);
    uArray         = [uArray u];
    xArray         = [xArray x];
    xbarArray      = [xbarArray xbar];    
    xhatArray      = [xhatArray xhat];    
    yArray         = [yArray y];
    thetabarArray  = [thetabarArray thetabar]; 
    thetahatArray  = [thetahatArray thetahat]; 

    % Simulate the system using Runge-Kutta
    k1 = x(3);
    l1 = x(4);
    m1 = alpha*x(2)-beta*x(3)+omega*u;
    n1 = -(alpha+gamma)*x(2)+delta*sin(x(1)+x(2))+beta*x(3)-omega*u;
    k2 = x(3)+0.5*dt*m1;
    l2 = x(4)+0.5*dt*n1;
    m2 = alpha*(x(2)+0.5*dt*l1)-beta*(x(3)+0.5*dt*m1)+omega*u;
    n2 = -(alpha+gamma)*(x(2)+0.5*dt*l1)+delta*sin((x(1)+0.5*dt*k1)+(x(2)+0.5*dt*l1))+beta*(x(3)+0.5*dt*m1)-omega*u;
    k3 = x(3)+0.5*dt*m2;
    l3 = x(4)+0.5*dt*n2;
    m3 = alpha*(x(2)+0.5*dt*l2)-beta*(x(3)+0.5*dt*m2)+omega*u;
    n3 = -(alpha+gamma)*(x(2)+0.5*dt*l2)+delta*sin((x(1)+0.5*dt*k2)+(x(2)+0.5*dt*l2))+beta*(x(3)+0.5*dt*m2)-omega*u;
    k4 = x(3)+dt*m3;
    l4 = x(4)+dt*n3;
    m4 = alpha*(x(2)+dt*l3)-beta*(x(3)+dt*m3)+omega*u;
    n4 = -(alpha+gamma)*(x(2)+dt*l3)+delta*sin((x(1)+dt*k3)+(x(2)+dt*l3))+beta*(x(3)+dt*m3)-omega*u;
    x(1) = x(1) + (dt/6)*(k1+2*k2+2*k3+k4);
    x(2) = x(2) + (dt/6)*(l1+2*l2+2*l3+l4);
    x(3) = x(3) + (dt/6)*(m1+2*m2+2*m3+m4);
    x(4) = x(4) + (dt/6)*(n1+2*n2+2*n3+n4);    
    y = C*x+dt*R^2*randn(n,1);

%    Phi = [y(3) 0 0 0 0 0 0 0 0;
%           0 y(4) 0 0 0 0 0 0 0;
%           0 0 y(2) y(3) u 0 0 0 0;
%           0 0 0 0 0 y(2) sin(y(1)+y(2)) y(3) u];
    
    Phi = [y(1) y(2) y(3) y(4) y(1)^2 y(2)^2 y(3)^2 y(4)^2 sin(y(1)+y(2)) sin(y(1)+y(3)) sin(y(1)+y(4)) sin(y(2)+y(3)) sin(y(2)+y(4)) sin(y(3)+y(4)) u zeros(45,1)';
           zeros(15,1)' y(1) y(2) y(3) y(4) y(1)^2 y(2)^2 y(3)^2 y(4)^2 sin(y(1)+y(2)) sin(y(1)+y(3)) sin(y(1)+y(4)) sin(y(2)+y(3)) sin(y(2)+y(4)) sin(y(3)+y(4)) u zeros(30,1)';
           zeros(30,1)' y(1) y(2) y(3) y(4) y(1)^2 y(2)^2 y(3)^2 y(4)^2 sin(y(1)+y(2)) sin(y(1)+y(3)) sin(y(1)+y(4)) sin(y(2)+y(3)) sin(y(2)+y(4)) sin(y(3)+y(4)) u zeros(15,1)';
           zeros(45,1)' y(1) y(2) y(3) y(4) y(1)^2 y(2)^2 y(3)^2 y(4)^2 sin(y(1)+y(2)) sin(y(1)+y(3)) sin(y(1)+y(4)) sin(y(2)+y(3)) sin(y(2)+y(4)) sin(y(3)+y(4)) u];
    
    % Estimation using adaptive observer
    Kx = Px*C'*inv(C*Px*C'+Rx);
    Kt = Pt*Gamma'*C'*inv(C*Gamma*Pt*Gamma'*C'+Rt);
    Gamma = (eye(n)-Kx*C)*Gamma;

    xbar = xbar+(Kx+Gamma*Kt)*(y-C*xbar);
    thetabar = thetabar-Kt*(y-C*xbar);

    xbar = A*xbar+Phi*thetabar;

    thetabar = thetabar;
    Px = (1/lambdav)*eye(n)*(eye(n)-Kx*C)*Px*eye(n);
    Pt = (1/lambdat)*(eye(r)-Kt*C*Gamma)*Pt;
    Gamma = eye(n)*Gamma-Phi;

    % Estimation using adaptive KF
    Pmin  = eye(n)*Pplus*eye(n)'+QF;
    Sigma = C*Pmin*C'+RF;
    KF    = Pmin*C'*inv(Sigma);
    Pplus = (eye(n)-KF*C)*Pmin;
     
    ytilde = y-C*xhat;
    QF    = a*QF + (1-a)*(KF*(ytilde*ytilde')*KF');    
    RF    = a*RF + (1-a)*(ytilde*ytilde'+C*Pmin*C');
 
    Upsilon = (eye(n)-KF*C)*(eye(n))*UpsilonPlus+(eye(n)-KF*C)*Phi;
    Omega   = C*eye(n)*UpsilonPlus+C*Phi;
    Lambda  = inv(lambda*Sigma+Omega*S*Omega');
    Gamma1  = S*Omega'*Lambda;
    S       = (1/lambda)*S-(1/lambda)*S*Omega'*Lambda*Omega*S;
    UpsilonPlus = Upsilon;
    
    thetahat  = thetahat + Gamma1*(y-C*xhat);
    xhat      = eye(n)*xhat+Phi*thetahat+KF*(y-C*xhat)+Upsilon*Gamma1*(y-C*xhat);

end

% figure(1)
% subplot(2,2,1)
% plot(t,yArray(1,:),'-','LineWidth',10);
% hold on;
% plot(t,xbarArray(1,:),':','LineWidth',10);
% hold on;
% plot(t,xhatArray(1,:),':','LineWidth',10);
% set(gca,'color','white','LineWidth',3,'FontSize',56)
% legend('measured','AO','AKF')
% grid on;
% grid minor;
% subplot(2,2,2)
% plot(t,yArray(2,:),'-','LineWidth',10);
% hold on;
% plot(t,xbarArray(2,:),':','LineWidth',10);
% hold on;
% plot(t,xhatArray(2,:),':','LineWidth',10);
% set(gca,'color','white','LineWidth',3,'FontSize',56)
% grid on;
% grid minor;
% subplot(2,2,3)
% plot(t,yArray(3,:),'-','LineWidth',10);
% hold on;
% plot(t,xbarArray(3,:),':','LineWidth',10);
% hold on;
% plot(t,xhatArray(3,:),':','LineWidth',10);
% set(gca,'color','white','LineWidth',3,'FontSize',56)
% grid on;
% grid minor;
% xlabel('t (s)')
% subplot(2,2,4)
% plot(t,yArray(4,:),'-','LineWidth',10);
% hold on;
% plot(t,xbarArray(4,:),':','LineWidth',10);
% hold on;
% plot(t,xhatArray(4,:),':','LineWidth',10);
% set(gca,'color','white','LineWidth',3,'FontSize',56)
% grid on;
% grid minor;
% xlabel('t (s)')

figure(2)
subplot(3,2,1)
plot(t,uArray,'-k','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
legend('control')
grid on;
grid minor;
ylabel('u','FontSize',72)
subplot(3,2,2)
plot(t,alpha*ones(1,length(t)),'-','LineWidth',10);
hold on;
plot(t,thetabarArray(32,:)/dt,':','LineWidth',10);
hold on;
plot(t,thetahatArray(32,:)/dt,':','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
legend('true parameter','WyNDA - AO','WyNDA - AKF')
ylim([0 1000])
grid on;
grid minor;
ylabel('\alpha','FontSize',72)
subplot(3,2,3)
plot(t,beta*ones(1,length(t)),'-','LineWidth',10);
hold on;
plot(t,-thetabarArray(33,:)/dt,':','LineWidth',10);
hold on;
plot(t,-thetahatArray(33,:)/dt,':','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
ylim([-200 200])
grid on;
grid minor;
ylabel('\beta','FontSize',72)
subplot(3,2,4)
plot(t,gamma*ones(1,length(t)),'-','LineWidth',10);
hold on;
plot(t,-thetabarArray(47,:)/dt-thetabarArray(32,:)/dt,':','LineWidth',10);
hold on;
plot(t,-thetahatArray(47,:)/dt-thetahatArray(32,:)/dt,':','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
ylim([100 400])
grid on;
grid minor;
ylabel('\gamma','FontSize',72)
subplot(3,2,5)
plot(t,delta*ones(1,length(t)),'-','LineWidth',10);
hold on;
plot(t,thetabarArray(54,:)/dt,':','LineWidth',10);
hold on;
plot(t,thetahatArray(54,:)/dt,':','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
ylim([-400 400])
grid on;
grid minor;
ylabel('\delta','FontSize',72)
xlabel('t (s)')
subplot(3,2,6)
plot(t,omega*ones(1,length(t)),'-','LineWidth',10);
hold on;
plot(t,thetabarArray(45,:)/dt,':','LineWidth',10);
hold on;
plot(t,thetahatArray(45,:)/dt,':','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
ylim([0 200])
grid on;
grid minor;
ylabel('\omega','FontSize',72)
xlabel('t (s)')