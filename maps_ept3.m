%=============================================================================| 
% # Copyright (C) 2023 Dr.-Ing. Arun Raina (E-Mail: arunraina@icloud.com)
%
% This matlab script is part of the code used for the paper, 
% "Analysis of hydrogen diffusion in the three stage electro-permeation test".
% DOI: https://doi.org/10.1007/s00161-023-01237-5
%
% This is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This code is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with this code. If not, see <http://www.gnu.org/licenses/>.
%=============================================================================|    
function maps_ept3
clc; clear all; close all;
format long e

tL01 = logspace(-7,-6,13);
tL02 = logspace(-6,-5,10);
tL03 = logspace(-5,-4,10);
tL04 = logspace(-4,-3,10);
tL05 = logspace(-3,-2,10);
tL0g = unique([tL01 tL02 tL03 tL04 tL05]);
bm   = size(tL0g,2);
NTg  = [1e-6; 1e-5; 1e-4; 1e-3];
btau1 = zeros(bm,4);
btau2 = zeros(bm,4);

for j = 1:4
    for i = 1:1
        [btau1v, btau2v] = bDH_bN_datF(tL0g(i),NTg(j));
        sprintf('%4d %4d',i,j)
        btau1(i,j) = btau1v;
        btau2(i,j) = btau2v;
    end
end
dlmwrite('btlag1K1e10.dat',btau1,'delimiter',' ','precision','%5.4d')
dlmwrite('btlag2K1e10.dat',btau2,'delimiter',' ','precision','%5.4d')

% =====================================================
% MAIN FUNCTION
% =====================================================
function [btau1, btau2] = bDH_bN_datF(tL0i,NTj)

% =====================================================
% define global quantities 
% =====================================================
global Q D0 R T0 l phi NL bet NT alp DH tL0 nx u02 u03 eta

% =====================================================
% List of input parameters
% =====================================================
Q   = 6680;             % Lattice enthalpy [J/mol]
D0  = 2.33e-7;          % Diffusion prefactor [m2/s]
R   = 8.3144;           % Gas constant [J/K/mol]
T0  = 293;              % Temperature [K]
m   = 0;                % PDE related plane slab
l   = 5e-3;             % Specimen thickness [m]
phi = 0.0;              % heating rate [K/sec]
NL  = 8.46e+28;         % NILS density [atoms/m3]      
bet = 1;                % No. of NILS per lattice atom
NT  = NL*NTj;           % Traps site density [atoms/m3]
K   = 1e10;             % <<<-- K = [1e10, 1e8, 1e6, 1e4]
DH  =-R*T0*log(K);      % Traps binding energy [J/mol]    
alp = 1;                % No. of atoms per trap site
tL0 = tL0i;             % Initial lat. occpancy ratio
eta = 1e-3;             % boundary condition in stage 2

% =====================================================
% approximate time [s] when steady state starts    
% =====================================================
Dapp = D0*exp(-Q/R/T0)/(1+NTj*K/(1+K*tL0));
tf   = l^2/Dapp*4;

% =====================================================
% Space & time discretization
% =====================================================
nx = 1e2;
nt = 2e4;
xa = linspace(-l/2,l/2 ,nx);
ta = linspace(0,tf,nt);

% =====================================================
% Normalised space and time
% =====================================================
x  = xa./l;
t  = (D0/l^2).*ta;
dx = x(2)-x(1);

% =====================================================
% time splitting
% =====================================================
z1 = 4000; z2 = 12000;
n1 = 1:z1; n2 = z1:z2; n3 = z2:nt;
t1 = t(n1); % stage1
t2 = t(n2); % stage2
t3 = t(n3); % stage3

% =====================================================
% boundary conditions for stages 1,2,3
% =====================================================
bc(n1) = 1; 
bc(n2) = eta;
bc(n3) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% STAGE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% =====================================================
% Solving ODE for \bar\theta_L (btL) in stage1
% =====================================================
btL(n1,:,1) = pdepe(m,@pdefun,@icfun1,@bcfun1,x,t1);

% =====================================================
% computing normalised flux at x = +l/2 at each time step
% =====================================================
% fileID = fopen('Q5','w');
for i = n1
    J(i) = -(btL(i,nx)-btL(i,nx-1))/dx;
%   fprintf(fileID,'%5.12e %5.12e\n',t(i),J(i));
end

% =====================================================
% Find normalised time constant when J=0.632
% =====================================================
btau1 = 0; jx1 = 0;
for i = n1
    if J(i)>=0.629
        btau1=t(i-1);
        jx1  =J(i-1);
        break
    end
end
if J(z1-9)-J(z1-10) > 1e-5 || J(z1) < 1e-2 || jx1 == 0
    disp('*** No steady state in stage 1 ***')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% STAGE 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% =====================================================
% initial condition for stage2
% =====================================================
u02(:,1) = btL(z1,:,1)'; 
u02(:,2) = x; 
% =====================================================
% Solving ODE for \bar\theta_L in stage2
% =====================================================
btL(n2,:,1) = pdepe(m,@pdefun,@icfun2,@bcfun2,x,t2);

% =====================================================
% computing normalised flux at x = +l/2 at each time step
% =====================================================
for i = n2
    J(i) = -(btL(i,nx)-btL(i,nx-1))/dx;
end
if J(z2-9)-J(z2-10) > 1e-5 || J(z2) > 1e-2
    disp('*** No steady state in stage 2 ***')
    sprintf('%4d ',J(z2))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% STAGE 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% =====================================================
% initial condition for stage3
% =====================================================
u03(:,1) = btL(z2,:,1)'; 
u03(:,2) = x; 

% =====================================================
% Solving ODE for \bar\theta_L in stage3
% =====================================================
btL(n3,:,1) = pdepe(m,@pdefun,@icfun3,@bcfun3,x,t3);

% =====================================================
% computing normalised flux at x = +l/2 at each time step
% =====================================================
for i = n3
    J(i) = -(btL(i,nx)-btL(i,nx-1))/dx;
end

% =====================================================
% Find normalised time constant when J=0.632
% =====================================================
btau2=0; jx2 = 0;
for i = n3
    if J(i)>=0.629
        btau2=t(i-1)-t(z2);
        jx2  =J(i-1);
        break
    end
end
btl2 = btau2+t(z2);
if J(nt-9)-J(nt-10) > 1e-5 || J(nt) < 1e-2 || jx2 == 0
    disp('*** No steady state in stage 3 ***')
end

% =====================================================
% PLOTTING
% =====================================================
L='LineWidth'; lw = 1.5';
gf1=figure('units','normalized','position',[.1 .5 .25 .4]);
subplot(2,1,1)
semilogx(t(1:z1-1),bc(1:z1-1),'r',t(z1:z2-1),bc(z1:z2-1),'b',...,
         t(z2:nt),bc(z2:nt),'g',L,lw)
ylabel('tL')
xlabel('tn')
ylim([0 1.2])
% xlim([1e-2 1e2])
% ax = gca; ax.XTick = [1e-2 1e0 1e2 1e3];
% decorate(gf1)
h1=legend('i','j','k');
h1.Position = [0.2,0.7,0.15,0.15];
h1.Box = 'off';
%--------------------------------------------------------
subplot(2,1,2)
semilogx(t(n1),J(n1),'r',t(n2),J(n2),'b',t(n3),J(n3),'g',...
         btau1,jx1,'sk', btl2,jx2,'sk',L,lw)
ylabel('Jn')
xlabel('tn')
ylim([0 1.2])
% xlim([1e-2 1e2])
% ax = gca; ax.XTick = [1e0 1e1 1e2 1e3 1e4 1e5];
% decorate(gf1)

clear global Q D0 R T0 l phi NL bet NT alp DH tL0 nx u02 u03 eta

% =====================================================
% PDE coefficients
% =====================================================
function [c,f,s] = pdefun(x,t,u,DuDx)
global Q D0 R T0 l phi NL bet NT alp DH tL0

% =====================================================
% Normalised quantities
% =====================================================
bQ   = Q/(R*T0);
bphi = phi*l^2/(T0*D0);
bT   = 1 + bphi*t;
bN   = (alp*NT)/(bet*NL);
bDH  = DH/(R*T0);
K   = exp(-bDH/bT); 
bDL  = exp(-bQ/bT);
%
ntraps = 1;
Disum  = 0; 
Stsum  = 0;
for i = 1:ntraps
    Disum = Disum + (bN*K)/(1+tL0*K*u)^2;
    Stsum = Stsum + (bN*K*bDH)/(1+tL0*K*u)^2;
end
Di = 1 + Disum;
St = -(bphi*u)/(bT^2)*Stsum;
%
c = Di; 
f = bDL*DuDx;
s = St;

% =====================================================
% Initial condition Stage 1
% =====================================================
function u0 = icfun1(x)
u0 = 0;

% =====================================================
% Initial condition Stage 2
% =====================================================
function u0 = icfun2(x)
global nx u02
for i = 1:nx
    if x == u02(i,2)
        u0 = u02(i);
    end
end

% =====================================================
% Initial condition Stage 3
% =====================================================
function u0 = icfun3(x)
global nx u03
for i = 1:nx
    if x == u03(i,2)
        u0 = u03(i);
    end
end

% =====================================================
% Boundary conditions Stage1
% =====================================================
function [pl,ql,pr,qr] = bcfun1(xl,ul,xr,ur,t)
pl = ul-1.0;
ql = 0;
pr = ur; 
qr = 0;

% =====================================================
% Boundary conditions Stage 2
% =====================================================
function [pl,ql,pr,qr] = bcfun2(xl,ul,xr,ur,t)
global eta
pl = ul-eta;
ql = 0;
pr = ur; 
qr = 0;

% =====================================================
% Boundary conditions Stage 3
% =====================================================
function [pl,ql,pr,qr] = bcfun3(xl,ul,xr,ur,t)
pl = ul-1.0;
ql = 0;
pr = ur; 
qr = 0;
