%-------------------------------------------------------------------------%
%   Script: Semi-analytical and analytical method to evaluate hydrodynamic 
%           coefficients (added mass and damping) of an surface piercing OWEC.
%
%   Article: Nguyen et al (2021) - Optimizing power generation of a 
%                                  bottom-raised oscillating surge wave 
%                                  energy converter using a theoretical model
%
%   Inputs: Environmental properties: T - wave periods
%                                     Aw - Wave height
%                                     h - water depth
%
%           OSWEC's Geometrical properties:   
%                                     w - width
%                                     c - distance from bottom
%                                     I - inertia
%                                     C - buoyant torque
%                                     
%-------------------------------------------------------------------------%
%   Written by: Jessica Nguyen, PhD 
%               University of Massachusetts Amherst
%               nvnguyen@umass.edu
%-------------------------------------------------------------------------%
clear; clc;

%% USER-DEFINED: environmental parameters;
% Environmental variables - from selected site;
envR.omega = linspace(0.1,20,400);
envR.T = 2*pi./envR.omega;

envR.Aw = 1;             % wave amplitudes;
envR.Asc = 1;
envR.beta = 0;           % wave direction [in degrees];
envR.depth = 1.0;            % water depth

% Constants;
envR.g = 9.81; 
envR.rho = 1000; 

% OSWEC geometric properties - OPTIMIZABLE PARAMETERS;
body.height = 0.5;        % Distance from bottom
body.width = 0.4;            % Device's width
body.volume = 0.001;
body.I55 = 0.07084;        % Device's moment of inertia
body.C55 = 0.3679;
body.Mass = 0.85;
body.Bv55 = 0.0;
body.Cext = 56;

% PTO system;
% nu_ptop           % Uncomments to manually set. Otherwise it is calculated based 
                    % on the optimizing values i.e. calculated from dP/dnu_pto = 0;  
pto.damping = 0;
pto.stiffness = 0;

% Convergent factors for solving dispersion (wn) and numerical scheme;
% Generally, convergence can be achieved with n = 15 and nmax = 15; 
solver.mode = 0;           % 0: no hinge. 1: hinge i.e. OSWEC device. 
solver.n = 15;             % Dispersion equation - number of frequencies to keep;
solver.nmax = 15;          % Collocation scheme for Chebyshev problem or Mathieu solutions;

%% Runs the solvers;
tic
[motion, hydroA, baseF] = optimizeHydro(solver, envR, body, pto);
toc

%% Get Capytaine data
hydroC = struct(); hydroW = struct();
wamitCase = "tr80";
if (solver.mode == 1)
    hydroW = readWAMIT(hydroW,'WAMIT\oswec_UMA\Output\'+wamitCase+'_hinge\wec.out','rao');
else
    hydroW = readWAMIT(hydroW,'WAMIT\oswec_UMA\Output\'+wamitCase+'\wec.out','rao');
end

%% Plotting results
marker = '*';
rightCo = [84.3, 18.8, 15.3]/100;
leftCo = [27.1, 45.9, 70.6]/100;

figure(1); 
set(gca,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
subplot(221);
yyaxis left; hold on;
plot(envR.omega, hydroA.A55, 'Color', leftCo)
plot(hydroW.w, squeeze(hydroW.A(5,5,:))*envR.rho, marker, 'Color', leftCo)
set(gca,'YColor','black')
ylabel('A_{55} (kg-m^2)')
xlabel('\omega (rad/s)')
ylim([-inf 1.2*max(abs(hydroA.A55))]);

yyaxis right; 
plot(envR.omega, hydroA.B55,'Color',rightCo)
plot(hydroW.w, squeeze(hydroW.B(5,5,:))*envR.rho.*hydroW.w',marker,'Color',rightCo)
ylim([0 1.2*max(abs(hydroA.B55))]);
set(gca,'YColor','black')
ylabel('B_{55} (kg-m^2/s)')
legend('','A55','','B55','Location','northeast')
hold off;

subplot(222);
yyaxis left; hold on;
plot(envR.omega, hydroA.A15,'Color', leftCo)
plot(hydroW.w, squeeze(hydroW.A(1,5,:))*envR.rho,marker, 'Color', leftCo)
set(gca,'YColor','black')
ylabel('A_{15} (kg-m^2)')
xlabel('\omega (rad/s)')
ylim([-inf 1.2*max(abs(hydroA.A15))]);

yyaxis right; 
plot(envR.omega, hydroA.B15,'Color',rightCo)
plot(hydroW.w, squeeze(hydroW.B(1,5,:))*envR.rho.*hydroW.w',marker,'Color',rightCo)
ylim([-inf 1.2*max(abs(hydroA.B15))]);
set(gca,'YColor','black')
ylabel('B_{15} (kg-m^2/s)')
legend('','A15','','B15','Location','northeast')
hold off;

subplot(223);
yyaxis left; hold on;
plot(envR.omega, abs(hydroA.X5),'Color', leftCo)
plot(hydroW.w, squeeze(hydroW.sc_ma(5,1,:))*envR.rho*envR.g,marker, 'Color', leftCo)
set(gca,'YColor','black')
ylabel('|X_5| (N-m)')
xlabel('\omega (rad/s)')
ylim([-inf 1.2*max(abs(hydroA.X5))]);

yyaxis right; hold on;
plot(envR.omega, wrapToPi(pi-angle(hydroA.X5)),'Color',rightCo)
plot(hydroW.w, squeeze(hydroW.sc_ph(5,1,:)),marker,'Color',rightCo)
set(gca,'YColor','black')
ylabel('\angle_5 (rad)')
legend('','|X_5|','','\angle_5','Location','northeast')
hold off;

subplot(224);
yyaxis left; hold on;
plot(envR.omega, abs(hydroA.X1),'Color', leftCo)
plot(hydroW.w, squeeze(hydroW.sc_ma(1,1,:))*envR.rho*envR.g,marker, 'Color', leftCo)
set(gca,'YColor','black')
ylabel('|X_1| (N))')
xlabel('\omega (rad/s)')
ylim([-inf 1.2*max(abs(hydroA.X1))]);

yyaxis right; hold on;
plot(envR.omega, wrapToPi(pi-angle(hydroA.X1)),'Color',rightCo)
plot(hydroW.w, squeeze(hydroW.sc_ph(1,1,:)),marker,'Color',rightCo)
set(gca,'YColor','black')
ylabel('\angle_1 (rad)')
legend('','|X_1|','','\angle_1','Location','northeast')
hold off;




















