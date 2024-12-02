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
% 
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
body.height = 0.5;      % OSWEC height
body.width = 0.4;       % Device's width
body.volume = 0.001;
body.I55 =  0.07084;    % Device's moment of inertia about the hinge
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
solver.mode = 1;           % 0: no hinge. 1: hinge i.e. OSWEC device. 
solver.n = 15;             % Dispersion equation - number of frequencies to keep;
solver.nmax = 15;          % Collocation scheme for Chebyshev problem or Mathieu solutions;

%% Runs the solvers;
tic
[motion, hydroA, baseF] = optimizeHydro(solver, envR, body, pto);
toc

%% Get WAMIT data
marker = '*';
leftCo = [27.1, 45.9, 70.6]/100;
rightCo = [84.3, 18.8, 15.3]/100;
lineStyle = ['-',"-.",'--',':'];
lineWidth = 0.9;
fontsize = 14;

wamitCase = ["tr10","tr20","tr40","tr80"];
wamitCase = fliplr(wamitCase);
diffPlot = 1;

for i = 1:length(wamitCase)
    hydroC = struct(); hydroW = struct();
    if (solver.mode == 1)
        hydroW = readWAMIT(hydroW,'WAMIT\oswec_UMA\Output\'+wamitCase(i)+'_hinge\wec.out','rao');
    else
        hydroW = readWAMIT(hydroW,'WAMIT\oswec_UMA\Output\'+wamitCase(i)+'\wec.out','rao');
    end

    %% Plotting results
    figure(1);
    set(figure(1),'defaultAxesColorOrder',[[0 0 0]; [0 0 0]],'Position', [100 100 1000 500]);
   
    subplot(121);
    yyaxis left; hold on;
    if (i==1 && diffPlot == 0)
        plot(envR.omega, hydroA.A55,'r', 'Marker',marker,...
            'MarkerIndices',1:15:length(hydroA.A55),'LineStyle','-','LineWidth',lineWidth)
        rightCo = leftCo;
    end
    wamitData = squeeze(hydroW.A(5,5,:))*envR.rho;
    if (diffPlot == 1)
        wamitData = (wamitData-hydroA.A55)./max(abs(wamitData))*100;
        ylabel('$\bar\varepsilon$ ($A_{55}$) (\%)','Interpreter','latex','FontSize',fontsize)
        ylim([-10 25]);
    else
        ylim([-inf 1.2*max(abs(hydroA.A55))]);
        ylabel('A_{55} (kg-m^2)','FontSize',fontsize)
    end
    plot(hydroW.w, wamitData, 'Color', leftCo, 'Marker', 'none',...
        'LineStyle',lineStyle(i),'LineWidth',lineWidth)
    xlabel('\omega (rad/s)','FontSize',fontsize)
    
    yyaxis right;
    if (i==1 && diffPlot == 0)
        plot(envR.omega, hydroA.B55,'r', 'Marker',marker,...
            'MarkerIndices',1:15:length(hydroA.B55),'LineStyle','-','LineWidth',lineWidth)
    end
    wamitData = squeeze(hydroW.B(5,5,:))*envR.rho.*hydroW.w';
    if (diffPlot == 1)
        wamitData = (wamitData-hydroA.B55)./max(abs(wamitData))*100;
        ylabel('$\bar\varepsilon$ ($B_{55}$) (\%)','Interpreter','latex','FontSize',fontsize)
        ylim([-10 25]);
    else
        ylim([-inf 1.2*max(abs(hydroA.B55))]);
        ylabel('B_{55} (kg-m^2/s)','FontSize',fontsize);
    end
    plot(hydroW.w, wamitData, 'Color', rightCo, 'Marker', 'none',...
        'LineStyle',lineStyle(i),'LineWidth',lineWidth)
    legend('$\bar\varepsilon$ ($A_{55}$)','','','','$\bar\varepsilon$ ($B_{55}$)',...
                'Interpreter','latex','Location','northeast','FontSize',fontsize)

    subplot(122);
    yyaxis left; hold on;
    if (i==1 && diffPlot == 0)
        plot(envR.omega, hydroA.A15, 'r','Marker',marker,...
            'MarkerIndices',1:15:length(hydroA.A15),'LineStyle','-','LineWidth',lineWidth)
    end
    wamitData = squeeze(hydroW.A(1,5,:))*envR.rho;
    if (diffPlot == 1)
        wamitData = (wamitData-hydroA.A15)./max(abs(wamitData))*100;
        ylabel('$\bar\varepsilon$ ($A_{15}$) (\%)','Interpreter','latex','FontSize',fontsize)
        if (solver.mode == 1)
            ylim([-10 25]);
        end
    else
        ylim([-inf 1.3*max(abs(hydroA.A15))]);
        ylabel('A_{15} (kg-m^2)','FontSize',fontsize)
    end
    plot(hydroW.w, wamitData,marker, ...
            'Color', leftCo, 'Marker', 'none', ...
            'LineStyle',lineStyle(i),'LineWidth',lineWidth)
    xlabel('\omega (rad/s)','FontSize',fontsize)

    yyaxis right;
    if (i == 1 && diffPlot == 0)
        plot(envR.omega, hydroA.B15, 'r','Marker',marker,...
            'MarkerIndices',1:15:length(hydroA.B15),'LineStyle','-','LineWidth',lineWidth)
    end
    wamitData = squeeze(hydroW.B(1,5,:))*envR.rho.*hydroW.w';
    if (diffPlot == 1)
        wamitData = (wamitData-hydroA.B15)./max(abs(wamitData))*100;
        ylabel('$\bar\varepsilon$ ($B_{15}$) (\%)','Interpreter','latex','FontSize',fontsize)
        ylim([-10 25]);
    else
        ylim([-inf 1.2*max(abs(hydroA.B15))]);
        ylabel('B_{15} (kg-m^2)','FontSize',fontsize)
    end
    plot(hydroW.w, wamitData, 'Color', rightCo, 'Marker', 'none',...
        'LineStyle',lineStyle(i),'LineWidth',lineWidth)
    xlabel('\omega (rad/s)','FontSize',fontsize)
    legend('$\bar\varepsilon$ ($A_{15}$)','','','','$\bar\varepsilon$ ($B_{15}$)',...
                'Interpreter','latex','Location','northeast','FontSize',fontsize)
    
    figure(2)
    set(figure(2),'defaultAxesColorOrder',[[0 0 0]; [0 0 0]],'Position', [100 100 1000 500]);
    subplot(121);
    yyaxis left; hold on;
    if (i==1 && diffPlot == 0)
        plot(envR.omega, abs(hydroA.X5), 'r', 'Marker',marker,...
            'MarkerIndices',1:15:length(hydroA.X5),'LineStyle','-','LineWidth',lineWidth)
    end
    wamitData = squeeze(hydroW.sc_ma(5,1,:))*envR.rho*envR.g;
    if (diffPlot == 1)
        wamitData = (wamitData-abs(hydroA.X5))./max(abs(wamitData))*100;
        ylabel('$\bar\varepsilon$ ($X_{5}$) (\%)','Interpreter','latex','FontSize',fontsize)
        xlim([0.5, inf]);
        if (solver.mode == 1)
            ylim([-10 12]);
        else
             ylim([-10 15]);           
        end
    else
        ylim([-inf 1.2*max(abs(hydroA.X5))]);
        ylabel('X_{5} (kg-m^2)','FontSize',fontsize)
    end
    plot(hydroW.w, wamitData,'Color', leftCo, 'Marker', 'none', ...
        'LineStyle',lineStyle(i),'LineWidth',lineWidth)
    xlabel('\omega (rad/s)','FontSize',fontsize)

    yyaxis right; hold on;
    anData = wrapToPi(pi-angle(hydroA.X5));
    if (i==1 && diffPlot == 0)
        plot(envR.omega, anData,'r', 'Marker',marker,...
            'MarkerIndices',1:15:length(anData),'LineStyle','-','LineWidth',lineWidth)
    end
    wamitData = squeeze(hydroW.sc_ph(5,1,:));
    if (diffPlot == 1)
        wamitData = (wamitData-anData)./max(abs(wamitData))*100;
        ylabel('$\bar\varepsilon$ ($\angle X_{5}$) (\%)','Interpreter','latex','FontSize',fontsize)
        xlim([0.5, inf]);
        if (solver.mode == 1)
            ylim([-10 12]);
        else
             ylim([-10 15]);   
        end
    else
        ylim([-inf 1.2*max(anData)]);
        ylabel('$\angle X_{5}$','Interpreter','latex','FontSize',fontsize)
    end
    plot(hydroW.w,wamitData,marker,'Color', rightCo, 'Marker', 'none', ...
        'LineStyle',lineStyle(i),'LineWidth',lineWidth)
    legend('$\bar\varepsilon$ ($X_{5}$)','','','','$\bar\varepsilon$ ($\angle X_{5}$)',...
                'Interpreter','latex','Location','northeast','FontSize',fontsize)

    subplot(122);
    yyaxis left; hold on;
    if (i==1 && diffPlot == 0)
        plot(envR.omega, abs(hydroA.X1), 'r', 'Marker',marker,...
            'MarkerIndices',1:15:length(hydroA.X1),'LineStyle','-','LineWidth',lineWidth)
    end
    wamitData = squeeze(hydroW.sc_ma(1,1,:))*envR.rho*envR.g;
    if (diffPlot == 1)
        wamitData = (wamitData-abs(hydroA.X1))./max(abs(wamitData))*100;
        ylabel('$\bar\varepsilon$ ($X_{1}$) (\%)','Interpreter','latex','FontSize',fontsize)
        if (solver.mode == 1)
            ylim([-10 12]);
        else
             ylim([-10 15]);   
        end
        xlim([0.5, inf]);
    else
        rightCo = leftCo;
        ylim([-inf 1.2*max(abs(hydroA.X1))]);
        ylabel('X_{1} (kg-m^2)','FontSize',fontsize)
    end
    plot(hydroW.w, wamitData,marker,'Color', leftCo, 'Marker', 'none', ...
        'LineStyle',lineStyle(i),'LineWidth',lineWidth)

    yyaxis right; hold on;
    anData = wrapToPi(pi-angle(hydroA.X1));
    if (i==1 && diffPlot == 0)
        plot(envR.omega, anData,'r', 'Marker',marker,...
            'MarkerIndices',1:15:length(anData),'LineStyle','-','LineWidth',lineWidth)
    end
    wamitData = squeeze(hydroW.sc_ph(1,1,:));
    if (diffPlot == 1)
        wamitData = (wamitData-anData)./max(abs(wamitData))*100;
        ylabel('$\bar\varepsilon$ ($\angle X_{1}$) (\%)','Interpreter','latex','FontSize',fontsize)
        xlim([0.5, inf]);
        if (solver.mode == 1)
            ylim([-10 12]);
        else
             ylim([-10 15]);   
        end
    else
        ylim([-inf 1.2*max(anData)]);
        ylabel('$\angle X_{1}$','Interpreter','latex','FontSize',fontsize)
    end
    plot(hydroW.w,wamitData,marker,'Color', rightCo, 'Marker', 'none',...
        'LineStyle',lineStyle(i),'LineWidth',lineWidth)
    xlabel('\omega (rad/s)','FontSize',fontsize)
    legend('$\bar\varepsilon$ ($X_{1}$)','','','','$\bar\varepsilon$ ($\angle X_{1}$)',...
                'Interpreter','latex','Location','northeast','FontSize',fontsize)
end











