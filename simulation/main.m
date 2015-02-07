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
% probability resolution (m)
drtil = 1e2;
% grid resolution
N1grid = 30;
N2grid = 20;
% boundary distance [-x1 +x1 -x2 +x2]
bdry = [-100 300 -150 150]*1e3; %m
% none-deflection probability
q = 1/pi;

%% incident to crash range pdf

Range = AC.Vc/gE/AC.SFC*AC.L2D*log(1+Fr*AC.Wfuel/AC.Wdry);
% engine fail, idk, idk, idk, nav error + hijack
Ncrash = [972, 548, 1226, 822, 109 + 1057];
Rcrash = [150, 50, 60, 30, Range]*1e3; %m

% stdev of the statistics (guess)
Rstdev = [.2 .2 .2 .2 .1];

Rprofile = [];
for i=1:numel(Ncrash)
    Rprofile = [Rprofile;
        randn(Ncrash(i),1)*Rstdev(i)*Rcrash(i) + Rcrash(i)];
end

% visuialization
% figure(); hist(Rprofile,100)

%% continuous probability at x (vector)
% rtil = 0:drtil:rint;
% Nrtil = numel(rtil);
% Rtil = @(x) norm(repmat(x,1,Nrtil)-[rtil;zeros(1,Nrtil)]);
% theta = @(x) abs(atan2(x(2),x(1)-rtil));
% ptheta = @(x) q-theta(x)*q/pi;
% % integral(pdf(Rtil)*pdf(rtil)*pdf(theta),0 to rint)
% px = @(x) sum(ksdensity(Rprofile,Rtil(x))*ptheta(x)/rint)*drtil;
% % px = @(x) sum(ksdensity(Rprofile,Rtil(x))/rint)*drtil;

Rtil = @(rtil,x) norm(x-[rtil;0]);
theta = @(rtil,x) abs(atan2(x(2),x(1)-rtil));
ptheta = @(rtil,x) q-theta(rtil,x)*q/pi;
pdfx = @(rtil,x) ksdensity(Rprofile,Rtil(rtil,x))*ptheta(rtil,x);
px = @(x) glq1d(pdfx,3,x);
%% discritized probability at cell (m,n)

Agrid = (bdry(2)-bdry(1))*(bdry(4)-bdry(3));
% grid pt construction
x1grid = linspace(bdry(1),bdry(2),N1grid)/1e3;
x2grid = linspace(bdry(3),bdry(4),N2grid)/1e3;
[X,Y] = meshgrid(x1grid,x2grid);

pxx = @(x1,x2) px([x1;x2]*1e3);
P = zeros(size(X));
for i = 1:numel(X)
%     P(i) = pxx(X(i),Y(i));
    
    P(i) = glq2ds(px);
end

figure(); hold all;
surf(X,Y,P,'edgecolor','none','facecolor','interp'); view(2); colorbar





