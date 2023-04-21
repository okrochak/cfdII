%%% Post-processing Centerline data
close all
clear

% Plotting preferences
fntSz = 20;
legSz = 15;
lblSz = 25;
lnWd = 1.5;

% Simulation conditions to be analysed at centerline
Re = 1000; % Reynolds number
Nruns  = [15, 31, 47, 55, 63]; % Defining Post processed runs
dt_mult = 1;
tol = 1e-10; % Simulation tolerance

%% Plot components along x = 0.5, y = 0.5;
for i = 1:5
    % Aquiring post processing results
    Crun = ['r','g','b','m','k']; % Assigning a colour to each run
    N = Nruns(i);
    C = Crun(i);
    load("results/PP_N="+N+"_tol="+string(tol)+"_Re="+Re+"_dtmult="+dt_mult+'.mat') % Loading Post Processed files
    
    % Variables along line x = 0.5
    Y_midline = postProc.x; % y- coordinates 
    p_ymidline = postProc.p(:,(N/2+0.5)); % Pressure along centerline
    vort_ymidline = postProc.vort(:,(N/2+0.5)); % Vorticity along centerline
    u_ycenterline =  postProc.u(:,(N/2+0.5)); % u Velocity along centerline
    v_ycenterline =  postProc.v(:,(N/2+0.5)); % v Velocity along centerline
    
    py = figure(1);
    set(gcf,'Position',[100 100 800 700])
    hold on
    plot(Y_midline,p_ymidline,C,'DisplayName',strcat('N = ',num2str(N)),'LineWidth',lnWd)
    legend('Location','NorthEast','interpreter','latex','FontSize',legSz)
    title('Pressure at $x = 0.5$','interpreter','latex','FontSize',lblSz);
    xlabel('$y$','interpreter','latex','FontSize',lblSz);
    ylabel('$p$','interpreter','latex','FontSize',lblSz);
    set(gca,'FontSize', fntSz,'TickLabelInterpreter','latex');

    vorty = figure(2);
    set(gcf,'Position',[100 100 800 700])
    hold on
    plot(Y_midline(1:end-1),vort_ymidline(1:end-1),C,'DisplayName',strcat('N = ',num2str(N)),'LineWidth',lnWd)
    legend('Location','NorthWest','interpreter','latex','FontSize',legSz)
    title('Vorticity at $x = 0.5$','interpreter','latex','FontSize',lblSz);
    xlabel('$y$','interpreter','latex','FontSize',lblSz);
    ylabel('$\omega$','interpreter','latex','FontSize',lblSz);
    set(gca,'FontSize', fntSz,'TickLabelInterpreter','latex');
    ylim([-5 17])

    uVely = figure(3);
    set(gcf,'Position',[100 100 800 700])
    hold on
    plot(Y_midline,u_ycenterline,C,'DisplayName',strcat('N = ',num2str(N)),'LineWidth',lnWd)
    legend('Location','NorthEast','interpreter','latex','FontSize',legSz)
    title('$x-$velocity at $x = 0.5$','interpreter','latex','FontSize',lblSz);
    xlabel('$y$','interpreter','latex','FontSize',lblSz);
    ylabel('$u$','interpreter','latex','FontSize',lblSz);
    set(gca,'FontSize', fntSz,'TickLabelInterpreter','latex');

    % Variables along line y = 0.5
    X_midline = postProc.x; % x- coordinates
    p_xmidline = postProc.p(N/2+0.5,:); % Pressure along centerline
    vort_xmidline = postProc.vort((N/2+0.5),:); % Vorticity along centerline
    u_xcenterline =  postProc.u(N/2+0.5,:); % u Velocity along centerline
    v_xcenterline =  postProc.v(N/2+0.5,:); % v Velocity along centerline
    
    px = figure(4);
    set(gcf,'Position',[100 100 800 700])
    hold on
    plot(Y_midline,p_xmidline,C,'DisplayName',strcat('N =',num2str(N)),'LineWidth',lnWd)
    legend('Location','SouthEast','interpreter','latex','FontSize',legSz)
    title('Pressure at $y = 0.5$','interpreter','latex','FontSize',lblSz);
    xlabel('$x$','interpreter','latex','FontSize',lblSz);
    ylabel('$p$','interpreter','latex','FontSize',lblSz);
    set(gca,'FontSize', fntSz,'TickLabelInterpreter','latex');

    vortx = figure(5);
    set(gcf,'Position',[100 100 800 700])
    hold on
    plot(Y_midline,vort_xmidline,C,'DisplayName',strcat('N = ',num2str(N)),'LineWidth',lnWd)
    legend('Location','SouthEast','interpreter','latex','FontSize',legSz)
    title('Vorticity at $y = 0.5$','interpreter','latex','FontSize',lblSz);
    xlabel('$x$','interpreter','latex','FontSize',lblSz);
    ylabel('$\omega$','interpreter','latex','FontSize',lblSz);
    set(gca,'FontSize', fntSz,'TickLabelInterpreter','latex');

    vVelx = figure(6);
    set(gcf,'Position',[100 100 800 700])
    hold on
    plot(Y_midline,v_xcenterline,C,'DisplayName',strcat('N = ',num2str(N)),'LineWidth',lnWd)
    legend('Location','SouthEast','interpreter','latex','FontSize',legSz)
    title('$y-$velocity at $y = 0.5$','interpreter','latex','FontSize',lblSz);
    xlabel('$x$','interpreter','latex','FontSize',lblSz);
    ylabel('$v$','interpreter','latex','FontSize',lblSz);
    set(gca,'FontSize', fntSz,'TickLabelInterpreter','latex');
end

% Reference results from Botella along centerline x = 0.5
Botella_vertdata = [1.0000 -1.00000 -1.0000000 0.052987 14.7534,;
0.9766 -0.65928 -0.6644227 0.052009 12.0670,;
0.9688 -0.57492 -0.5808359 0.051514 9.49496,;
0.9609 -0.51117 -0.5169277 0.050949 6.95968,;
0.9531 -0.46604 -0.4723329 0.050329 4.85754,;
0.8516 -0.33304 -0.3372212 0.034910 1.76200,;
0.7344 -0.18719 -0.1886747 0.012122 2.09121,;
0.6172 -0.05702 -0.0570178 -0.000827 2.06539,;
0.5000 0.06080 0.0620561 0.000000 2.06722,;
0.4531 0.10648 0.1081999 0.004434 2.06215,;
0.2813 0.27805 0.2803696 0.040377 2.26772,;
0.1719 0.38289 0.3885691 0.081925 1.05467,;
0.1016 0.29730 0.3004561 0.104187 -1.63436,;
0.0703 0.22220 0.2228955 0.108566 -2.20175,;
0.0625 0.20196 0.2023300 0.109200 -2.31786,;
0.0547 0.18109 0.1812881 0.109689 -2.44960,;
0.0000 0.00000 0.0000000 0.110591 -4.16648];

figure(py);
hold on
scatter(Botella_vertdata(:,1),Botella_vertdata(:,4),'k','DisplayName','Reference','LineWidth',lnWd);
figure(vorty);
hold on
scatter(Botella_vertdata(:,1),Botella_vertdata(:,5),'k','DisplayName','Reference','LineWidth',lnWd);
figure(uVely);
hold on
scatter(Botella_vertdata(:,1),Botella_vertdata(:,3),'k','DisplayName','Reference','LineWidth',lnWd);

% Reference results from Botella along centerline y = 0.5
Botella_hordata = [0.0000 0.00000 0.0000000 0.077455 -5.46217,;
0.0312 -0.21388 -0.2279225 0.078837 -8.44350,;
0.0391 -0.27669 -0.2936869 0.078685 -8.24616,;
0.0469 -0.33714 -0.3553213 0.078148 -7.58524,;
0.0547 -0.39188 -0.4103754 0.077154 -6.50867,;
0.0937 -0.51550 -0.5264392 0.065816 0.92291,;
0.1406 -0.42665 -0.4264545 0.049029 3.43016,;
0.1953 -0.31966 -0.3202137 0.034552 2.21171,;
0.5000 0.02526 0.0257995 0.000000 2.06722,;
0.7656 0.32235 0.3253592 0.044848 2.06122,;
0.7734 0.33075 0.3339924 0.047260 2.00174,;
0.8437 0.37095 0.3769189 0.069511 0.74207,;
0.9062 0.32627 0.3330442 0.084386 -0.82398,;
0.9219 0.30353 0.3099097 0.086716 -1.23991,;
0.9297 0.29012 0.2962703 0.087653 -1.50306,;
0.9375 0.27485 0.2807056 0.088445 -1.83308,;
1.0000 0.00000 0.0000000 0.090477 -7.66369];

figure(px);
hold on
scatter(Botella_hordata(:,1),Botella_hordata(:,4),'k','DisplayName','Reference','LineWidth',lnWd);
figure(vortx);
hold on
scatter(Botella_hordata(:,1),Botella_hordata(:,5),'k','DisplayName','Reference','LineWidth',lnWd);
figure(vVelx);
hold on
scatter(Botella_hordata(:,1),Botella_hordata(:,3),'k','DisplayName','Reference','LineWidth',lnWd);

exportgraphics(figure(1),["figures/X05pressure_tol="+string(tol)+"_Re="+Re+"_dtmult="+dt_mult+".pdf"], 'Resolution', 300)
exportgraphics(figure(2),["figures/X05vorticity_tol="+string(tol)+"_Re="+Re+"_dtmult="+dt_mult+".pdf"], 'Resolution', 300)
exportgraphics(figure(3),["figures/X05xvelocity_tol="+string(tol)+"_Re="+Re+"_dtmult="+dt_mult+".pdf"], 'Resolution', 300)
exportgraphics(figure(4),["figures/Y05pressure_tol="+string(tol)+"_Re="+Re+"_dtmult="+dt_mult+".pdf"], 'Resolution', 300)
exportgraphics(figure(5),["figures/Y05vorticity_tol="+string(tol)+"_Re="+Re+"_dtmult="+dt_mult+".pdf"], 'Resolution', 300)
exportgraphics(figure(6),["figures/Y05yvelocity_tol="+string(tol)+"_Re="+Re+"_dtmult="+dt_mult+".pdf"], 'Resolution', 300)