%% Research code by Agus Hasan

clear;
clc;

%% time horizon
tf  = 6;
dt  = 0.001;
t   = dt:dt:tf;

%% number of variables and coefficients
n = 3;
r = 48;

%% system description
A = eye(n);
C = eye(n);

%% noise
R = 1;

%% state initialization
x        = [-8;7;27];
xbar     = x;
xhat     = x;
y        = x;
thetabar = zeros(r,1);
thetahat = zeros(r,1);
 
%% true parameters
sigma   = 10;
rho     = 28;
beta    = 8/3;

%% initial control inputs
%u     = [10 1]';

%% for plotting
%uArray          = [];
xArray          = [];
xbarArray       = [];
xhatArray       = [];
yArray          = [];
thetabarArray   = [];
thetahatArray   = [];

%% Initialization for estimator

lambdav = 0.995;
lambdat = 0.999;
Rx      = 1*eye(n);
Rt      = 1*eye(n);
Px      = 0.1*eye(n);
Pt      = 0.1*eye(r);
Gamma   = 1*zeros(n,r);

Pplus       = 10000*eye(n);
QF          = 1*eye(n);
RF          = 1*eye(n);
a           = 0.999;
UpsilonPlus = 0*zeros(n,r);
S           = 0.1*eye(r);
lambda      = 0.999;

%% simulation
for i=1:(tf/dt)

    %uArray         = [uArray u];
    xArray         = [xArray x];
    xbarArray      = [xbarArray xbar];    
    xhatArray      = [xhatArray xhat];    
    yArray         = [yArray y];
    thetabarArray  = [thetabarArray thetabar]; 
    thetahatArray  = [thetahatArray thetahat]; 

    % Simulate the system using Runge-Kutta
    k1 = sigma*(x(2)-x(1));
    l1 = x(1)*(rho-x(3))-x(2);
    m1 = x(1)*x(2)-beta*x(3);
    k2 = sigma*((x(2)+0.5*dt*l1)-(x(1)+0.5*dt*k1));
    l2 = (x(1)+0.5*dt*k1)*(rho-(x(3)+0.5*dt*m1))-(x(2)+0.5*dt*l1);
    m2 = (x(1)+0.5*dt*k1)*(x(2)+0.5*dt*l1)-beta*(x(3)+0.5*dt*m1);
    k3 = sigma*((x(2)+0.5*dt*l2)-(x(1)+0.5*dt*k2));
    l3 = (x(1)+0.5*dt*k2)*(rho-(x(3)+0.5*dt*m2))-(x(2)+0.5*dt*l2);
    m3 = (x(1)+0.5*dt*k2)*(x(2)+0.5*dt*l2)-beta*(x(3)+0.5*dt*m2);
    k4 = sigma*((x(2)+dt*l3)-(x(1)+dt*k3));
    l4 = (x(1)+dt*k3)*(rho-(x(3)+dt*m3))-(x(2)+dt*l3);
    m4 = (x(1)+dt*k3)*(x(2)+dt*l3)-beta*(x(3)+dt*m3);
    x(1) = x(1) + (dt/6)*(k1+2*k2+2*k3+k4);
    x(2) = x(2) + (dt/6)*(l1+2*l2+2*l3+l4);
    x(3) = x(3) + (dt/6)*(m1+2*m2+2*m3+m4);
    y = C*x+dt*R^2*randn(n,1);

    Phi = [1 y(1) y(2) y(3) y(1)^2 y(2)^2 y(3)^2 y(1)*y(2) y(1)*y(3) y(2)*y(3) sin(y(1)) sin(y(2)) sin(y(3)) cos(y(1)) cos(y(2)) cos(y(3)) zeros(32,1)';
           zeros(16,1)' 1 y(1) y(2) y(3) y(1)^2 y(2)^2 y(3)^2 y(1)*y(2) y(1)*y(3) y(2)*y(3) sin(y(1)) sin(y(2)) sin(y(3)) cos(y(1)) cos(y(2)) cos(y(3)) zeros(16,1)';
           zeros(32,1)' 1 y(1) y(2) y(3) y(1)^2 y(2)^2 y(3)^2 y(1)*y(2) y(1)*y(3) y(2)*y(3) sin(y(1)) sin(y(2)) sin(y(3)) cos(y(1)) cos(y(2)) cos(y(3))];
    
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

figure(1)
plot3(yArray(1,:),yArray(2,:),yArray(3,:),'-','LineWidth',16);
hold on;
plot3(xbarArray(1,:),xbarArray(2,:),xbarArray(3,:),':','LineWidth',16)
hold on;
plot3(xhatArray(1,:),xhatArray(2,:),xhatArray(3,:),':','LineWidth',16)
legend('measurement','WyNDA - AO','WyNDA - AKF')
set(gca,'color','white','LineWidth',3,'FontSize',56)
grid on;
grid minor;
xlabel('p','FontSize',72)
ylabel('q','FontSize',72)
zlabel('r','FontSize',72)

figure(2)
subplot(3,2,1)
plot(t,sigma*ones(1,length(t)),'-','LineWidth',10);
hold on;
plot(t,-thetabarArray(2,:)/dt,':','LineWidth',10);
hold on;
plot(t,-thetahatArray(2,:)/dt,':','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
legend('true parameter','WyNDA - AO','WyNDA - AKF')
grid on;
grid minor;
ylim([0 20])
ylabel('\sigma','FontSize',72)
subplot(3,2,2)
plot(t,-sigma*ones(1,length(t))-thetabarArray(2,:)/dt,':k','LineWidth',10);
hold on;
plot(t,-sigma*ones(1,length(t))-thetahatArray(2,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
legend('error WyNDA - AO','error WyNDA - AKF')
grid on;
grid minor;
ylabel('error \sigma','FontSize',72)
subplot(3,2,3)
plot(t,rho*ones(1,length(t)),'-','LineWidth',10);
hold on;
plot(t,thetabarArray(18,:)/dt,':','LineWidth',10);
hold on;
plot(t,thetahatArray(18,:)/dt,':','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylim([25 30])
ylabel('\rho','FontSize',72)
subplot(3,2,4)
plot(t,rho*ones(1,length(t))-thetabarArray(18,:)/dt,':k','LineWidth',10);
hold on;
plot(t,rho*ones(1,length(t))-thetahatArray(18,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('error \rho','FontSize',72)
subplot(3,2,5)
plot(t,beta*ones(1,length(t)),'-','LineWidth',10);
hold on;
plot(t,-thetabarArray(36,:)/dt,':','LineWidth',10);
hold on;
plot(t,-thetahatArray(36,:)/dt,':','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylim([0 6])
ylabel('\beta','FontSize',72)
xlabel('t (s)')
subplot(3,2,6)
plot(t,-beta*ones(1,length(t))-thetabarArray(36,:)/dt,':k','LineWidth',10);
hold on;
plot(t,-beta*ones(1,length(t))-thetahatArray(36,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('error \beta','FontSize',72)
xlabel('t (s)')

figure(3)
plot3(yArray(1,:),yArray(2,:),yArray(3,:),'-k','LineWidth',10);
%legend('measured','AO','AKF')
set(gca,'color','white','LineWidth',3,'FontSize',56)
grid on;
grid minor;
%xlabel('p','FontSize',72)
%ylabel('q','FontSize',72)
%zlabel('r','FontSize',72)

Coeff = round([(1/dt)*thetabar(1:(r/n),end)'; (1/dt)*thetabar((r/n)+1:2*(r/n),end)'; (1/dt)*thetabar(2*(r/n)+1:r,end)'])