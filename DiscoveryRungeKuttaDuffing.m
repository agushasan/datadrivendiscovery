%% Code by Agus Hasan

clear;
clc;

%% System parameters

n = 2;
r = 20;

alpha   = 0;
beta    = 1;
delta   = 0.05;
gamma   = 8;

tf = 15;
dt = 0.01;
t  = dt:dt:tf;

A   = [0 1; -alpha -delta];
B   = [0 gamma]';
C   = eye(2);
u   = 1;

R   = 2;

%% Initialization
x               = [2 3]';
xArray          = [];
xbar            = [2 3]';
xbarArray       = [];
xhat            = [2 3]';
xhatArray       = [];
theta           = zeros(1,r)';
thetaArray      = [];
thetabar        = zeros(1,r)';
thetabarArray   = [];
thetahat        = zeros(1,r)';
thetahatArray   = [];
y               = x;
yArray          = [];
uArray          = [];

lambdax = 0.9;
lambdat = 0.9985;
Rx = 1*eye(n);
Px = 1*eye(n);
Rt = 1*eye(n);
Pt = 1*eye(r);
Gamma = zeros(n,r);

Pplus       = 10000000*eye(n);
QF          = 0.1*eye(n);
RF          = 0.1*eye(n);
a           = 0.999;
UpsilonPlus = 1*zeros(n,r);
S           = 1*eye(r);
lambda      = 0.995;

%% Simulation

for i = 1:(tf/dt)

%% control

    u = cos(i*dt);
    yArray          = [yArray y];
    xArray          = [xArray x];
    xbarArray       = [xbarArray xbar];
    xhatArray       = [xhatArray xhat];    
    thetaArray      = [thetaArray theta];
    thetabarArray   = [thetabarArray thetabar];
    thetahatArray   = [thetahatArray thetahat];
    uArray          = [uArray u];

    % Simulating the system using Runge-Kutta
    k1 = x(2);
    l1 = -delta*x(2)-alpha*x(1)-beta*x(1)^3+gamma*u;
    k2 = x(2)+0.5*dt*l1;
    l2 = -delta*(x(2)+0.5*dt*l1)-alpha*(x(1)+0.5*dt*k1)-beta*(x(1)+0.5*dt*k1)^3+gamma*u;
    k3 = x(2)+0.5*dt*l2;
    l3 = -delta*(x(2)+0.5*dt*l2)-alpha*(x(1)+0.5*dt*k2)-beta*(x(1)+0.5*dt*k2)^3+gamma*u;
    k4 = x(2)+dt*l3;
    l4 = -delta*(x(2)+dt*l3)-alpha*(x(1)+dt*k3)-beta*(x(1)+dt*k3)^3+gamma*u;
    x(1) = x(1) + (dt/6)*(k1+2*k2+2*k3+k4);
    x(2) = x(2) + (dt/6)*(l1+2*l2+2*l3+l4);
    % Taking measurement
    y = C*x+dt*R^2*randn(n,1);

    Phi = [1 y(1) y(2) y(1)^2 y(2)^2 y(1)*y(2) y(1)^3 y(2)^3 u u^2 zeros(1,10);
          zeros(1,10) 1 y(1) y(2) y(1)^2 y(2)^2 y(1)*y(2) y(1)^3 y(2)^3 u u^2];
    
    % Estimation using adaptive observer
    Kx = Px*C'*inv(C*Px*C'+Rx);
    Kt = Pt*Gamma'*C'*inv(C*Gamma*Pt*Gamma'*C'+Rt);
    Gamma = (eye(n)-Kx*C)*Gamma;
    xbar = xbar+(Kx+Gamma*Kt)*(y-C*xbar);
    thetabar = thetabar-Kt*(y-C*xbar);

    xbar = eye(n)*xbar+Phi*thetabar;
    thetabar = thetabar;
    Px = (1/lambdax)*eye(n)*(eye(n)-Kx*C)*Px*eye(n);
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
subplot(2,1,1)
plot(t,yArray(1,:),'-','LineWidth',10)
set(gca,'color','white','LineWidth',3,'FontSize',36)
legend('measurement')
grid on;
grid minor;
ylabel('p [m]','FontSize',72)
xlim([0 tf]);
subplot(2,1,2)
plot(t,yArray(2,:),'-','LineWidth',10)
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
xlim([0 tf]);
ylabel('q [m/s]','FontSize',72)
xlabel('time (s)')

figure(2)
subplot(2,1,1)
plot(t,yArray(1,:),'-','LineWidth',10)
hold on;
plot(t,xbarArray(1,:),':','LineWidth',10)
hold on;
plot(t,xhatArray(1,:),':','LineWidth',10)
set(gca,'color','white','LineWidth',3,'FontSize',36)
legend('measurement','WyNDA - AO','WyNDA - AKF')
grid on;
grid minor;
ylabel('p [m]','FontSize',72)
xlim([0 tf]);
subplot(2,1,2)
plot(t,yArray(2,:),'-','LineWidth',10)
hold on;
plot(t,xbarArray(2,:),':','LineWidth',10)
hold on;
plot(t,xhatArray(2,:),':','LineWidth',10)
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
xlim([0 tf]);
ylabel('q [m/s]','FontSize',72)
xlabel('time (s)')

figure(3)
subplot(4,2,1)
plot(t,alpha*ones(length(t),1),'-','LineWidth',10)
hold on;
plot(t,-thetabarArray(12,:)/dt,':','LineWidth',10)
hold on;
plot(t,-thetahatArray(12,:)/dt,':','LineWidth',10)
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
xlim([0 tf]);
ylim([alpha-1 alpha+1]);
ylabel('\alpha','FontSize',72)
legend('true parameter','WyNDA - AO','WyNDA - AKF')
subplot(4,2,2)
plot(t,alpha*ones(length(t),1)-(-thetabarArray(12,:)/dt)',':k','LineWidth',10)
hold on;
plot(t,alpha*ones(length(t),1)-(-thetahatArray(12,:)/dt)',':g','LineWidth',10)
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
xlim([0 tf]);
ylabel('error \alpha','FontSize',72)
legend('error WyNDA - AO','error WyNDA - AKF')
subplot(4,2,3)
plot(t,beta*ones(length(t),1),'-','LineWidth',10)
hold on;
plot(t,-thetabarArray(17,:)/dt,':','LineWidth',10)
hold on;
plot(t,-thetahatArray(17,:)/dt,':','LineWidth',10)
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
xlim([0 tf]);
ylim([beta-1 beta+1]);
ylabel('\beta','FontSize',72)
subplot(4,2,4)
plot(t,beta*ones(length(t),1)-(-thetabarArray(17,:)/dt)',':k','LineWidth',10)
hold on;
plot(t,beta*ones(length(t),1)-(-thetahatArray(17,:)/dt)',':g','LineWidth',10)
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
xlim([0 tf]);
ylabel('error \beta','FontSize',72)
subplot(4,2,5)
plot(t,gamma*ones(length(t),1),'-','LineWidth',10)
hold on;
plot(t,thetabarArray(19,:)/dt,':','LineWidth',10)
hold on;
plot(t,thetahatArray(19,:)/dt,':','LineWidth',10)
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
xlim([0 tf]);
ylim([gamma-1 gamma+1]);
ylabel('\gamma','FontSize',72)
subplot(4,2,6)
plot(t,gamma*ones(length(t),1)-(thetabarArray(19,:)/dt)',':k','LineWidth',10)
hold on;
plot(t,gamma*ones(length(t),1)-(thetahatArray(19,:)/dt)',':g','LineWidth',10)
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
xlim([0 tf]);
ylabel('error \gamma','FontSize',72)
subplot(4,2,7)
plot(t,delta*ones(length(t),1),'-','LineWidth',10)
hold on;
plot(t,-thetabarArray(13,:)/dt,':','LineWidth',10)
hold on;
plot(t,-thetahatArray(13,:)/dt,':','LineWidth',10)
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
xlim([0 tf]);
ylim([delta-1 delta+1]);
ylabel('\delta','FontSize',72)
xlabel('time (s)')
subplot(4,2,8)
plot(t,delta*ones(length(t),1)-(-thetabarArray(13,:)/dt)',':k','LineWidth',10)
hold on;
plot(t,delta*ones(length(t),1)-(-thetahatArray(13,:)/dt)',':g','LineWidth',10)
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
xlim([0 tf]);
ylabel('error \delta','FontSize',72)
xlabel('time (s)')

RMSEAO  = rms(alpha*ones(length(t),1)-(-thetabarArray(12,:)/dt)')+rms(beta*ones(length(t),1)-(-thetabarArray(17,:)/dt)')+rms(gamma*ones(length(t),1)-(thetabarArray(19,:)/dt)')+rms(delta*ones(length(t),1)-(-thetabarArray(13,:)/dt)')
RMSEAKF = rms(alpha*ones(length(t),1)-(-thetahatArray(12,:)/dt)')+rms(beta*ones(length(t),1)-(-thetahatArray(17,:)/dt)')+rms(gamma*ones(length(t),1)-(thetahatArray(19,:)/dt)')+rms(delta*ones(length(t),1)-(-thetahatArray(13,:)/dt)')

Coeff1 = round([(1/dt)*thetabar(1:(r/n),end)'; (1/dt)*thetabar((r/n)+1:2*(r/n),end)'],2)
Coeff2 = round([(1/dt)*thetahat(1:(r/n),end)'; (1/dt)*thetahat((r/n)+1:2*(r/n),end)'],2)