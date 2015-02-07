clc; clear; close all;
% constants
gE = 9.81; %m/s2
%% plane specs (B737-900ER; G280; A380)
% Specific Fuel Consumption (kg/N/s)
B737.SFC = 20E-6;
G280.SFC = 20E-6;
A380.SFC = 20E-6;
% Lift to Drag Ratio
B737.L2D = 20;
G280.L2D = 20;
A380.L2D = 20;
% Dry Weight = Operation Empty Weight + Total Payload (kg)
B737.Wdry = 40E3+25E3;
G280.Wdry = 40E3+25E3;
A380.Wdry = 40E3+25E3;
% Total fuel load = Maximum Take Off Weight - Dry Weight (kg)
B737.Wfuel = 85E3 - B737.Wdry;
G280.Wfuel = 85E3 - G280.Wdry;
A380.Wfuel = 85E3 - A380.Wdry;
% Cruise speed (m/s)
B737.Vc = 217;
G280.Vc = 217;
A380.Vc = 217;

%% Assumptions/parameters

% lost aircraft model
AC = B737;
% remaining fuel ratio at incident
Fr = .5;
% communication interval/distance
interval = 15*60; %s
rint = AC.Vc*interval; %m
% probability resolution
Nrtil = 1e1; %m
% grid resolution
N1grid = 30;
N2grid = 20;
% boundary distance [-x1 +x1 -x2 +x2]
bdry = [-150 350 -150 150]*1e3; %m
% reverse probability
s = .05; q = 1/pi-2*s; 

%% incident to crash range pdf

Range = AC.Vc/gE/AC.SFC*AC.L2D*log(1+Fr*AC.Wfuel/AC.Wdry);
% engine fail, idk, idk, idk, nav error + hijack
Ncrash = [972, 548, 1226, 822, 109 + 1057];
Rcrash = [150, 50, 60, 30, Range]*1e3; %m

% stdev of the statistics (guess)
Rstdev = [.3 .2 .2 .2 .2];

Rprofile = [];
for i=1:numel(Ncrash)
    Rprofile = [Rprofile;
        randn(Ncrash(i),1)*Rstdev(i)*Rcrash(i) + Rcrash(i)];
end

% visuialization
% hist(Rprofile)
% [a,b]=ksdensity(Rprofile(Rprofile<500e3));
% plot(b,a)

%% continuous probability at x (vector)
rtil = linspace(0,rint,Nrtil);
% Rtil = @(x) norm(repmat(x,1,Nrtil)-[rtil;zeros(1,Nrtil)]);
% theta = @(x) abs(atan2(x(2),x(1)-rtil));
% ptheta = @(x) q+theta(x)*(s-q)/pi;
% integrate pdf(Rtil)*pdf(rtil)*pdf(theta) from rtil=0 to rint
% px = @(x) sum(ksdensity(Rprofile,Rtil(x),'function','pdf')*ptheta(x)/rint)*rint/Nrtil;


% Rtil = @(rtil,x) norm([x(1)-rtil, repmat(x(2),1,numel(rtil))]);
% theta = @(rtil,x) abs(atan2(x(2),x(1)-rtil));
% ptheta = @(rtil,x) q+theta(rtil,x)*(s-q)/pi;
% pdfx = @(rtil,x) ksdensity(Rprofile,Rtil(rtil,x))*ptheta(rtil,x);

Rtil = @(rtil,x) norm(x-[rtil;0]);
theta = @(rtil,x) abs(atan2(x(2),x(1)-rtil));
ptheta = @(rtil,x) q+theta(rtil,x)*(s-q)/pi;
pdfx = @(rtil,x) ksdensity(Rprofile,Rtil(rtil,x))*ptheta(rtil,x);
%% discritized probability at cell (m,n)

Agrid = (bdry(2)-bdry(1))*(bdry(4)-bdry(3));
% grid pt construction
x1grid = linspace(bdry(1),bdry(2),N1grid);
x2grid = linspace(bdry(3),bdry(4),N2grid);
[X,Y] = meshgrid(x1grid,x2grid);

P = zeros(size(X));
for i = 1:numel(X)
    for ii = 1:Nrtil
        P(i) = P(i) + pdfx(rtil(ii),[X(i);Y(i)])*rint/Nrtil;
    end
end
P = P / sum(P(:)); % renormalization

%% crash distribution graph

S = Surface(N1grid,N2grid,bdry/1e3);
figure(); hold all; grid on;
plottwoform(S,P',3); colorbar;
xlabel('Tangent Direction [km]'); ylabel('Lateral Direction [km]');
title('Probability of Aircraft Landed in Cell');
hold all;
plot3([0 rint]/1e3, [0 1], [1 1],'rx--');


%% drifting simulation




%% lateral range curve/function