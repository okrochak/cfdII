%%% Post-processing Centerline data
close all
clear

% Plotting preferences
fntSz = 20;
lblSz = 25;
lnWd = 2;

% Simulation conditions
tol = 1e-10;
Re = 1000;

%% Plot components along x = 0.5, y = 0.5;
for i = 1:5
    % Aquiring post processing results
    Nruns  = [15, 31, 47, 55, 63]; % Defining Post processed runs
    Crun = ['r','g','b','m','k']; % Assigning a colour to each run
    N = Nruns(i);
    C = Crun(i);
    load(strcat('results/PP_N=',num2str(N),'_tol=',num2str(tol),'_Re=',num2str(Re),'.mat')) %Loading Post Processed files
    
    % Variables along line x = 0.5
    Y_midline = postProc.x; % y- coordinates 
    p_ymidline = postProc.p(:,(N/2+0.5)); % Pressure along centerline
    vort_ymidline = postProc.vort(:,(N/2+0.5)); % Vorticity along centerline
    u_ycenterline =  postProc.u(:,(N/2+0.5)); % u Velocity along centerline
    v_ycenterline =  postProc.v(:,(N/2+0.5)); % v Velocity along centerline
    
    figure(1)
    py = subplot(1,3,1);
    hold on
    plot(Y_midline,p_ymidline,C,'DisplayName',strcat('N = ',num2str(N)))
    legend()
    title('Pressure','interpreter','latex','FontSize',lblSz);
    xlabel('$y$','interpreter','latex','FontSize',lblSz);
    ylabel('$p$','interpreter','latex','FontSize',lblSz);

    vorty = subplot(1,3,2);
    hold on
    plot(Y_midline(1:end-1),vort_ymidline(1:end-1),C,'DisplayName',strcat('N = ',num2str(N)))
    legend()
    title('Vorticity','interpreter','latex','FontSize',lblSz);
    xlabel('$y$','interpreter','latex','FontSize',lblSz);
    ylabel('$\omega$','interpreter','latex','FontSize',lblSz);

    uVely = subplot(1,3,3);
    hold on
    plot(Y_midline,u_ycenterline,C,'DisplayName',strcat('N = ',num2str(N)))
    legend()
    title('Velocity','interpreter','latex','FontSize',lblSz);
    xlabel('$y$','interpreter','latex','FontSize',lblSz);
    ylabel('$u$','interpreter','latex','FontSize',lblSz);

    sgtitle('Pressure, vorticity and horizontal velocity along centerline x = 0.5') 
    set(gcf,'Position',[100 100 1500 700])

    % Variables along line y = 0.5
    X_midline = postProc.x; % x- coordinates
    p_xmidline = postProc.p(N/2+0.5,:); % Pressure along centerline
    vort_xmidline = postProc.vort((N/2+0.5),:); % Vorticity along centerline
    u_xcenterline =  postProc.u(N/2+0.5,:); % u Velocity along centerline
    v_xcenterline =  postProc.v(N/2+0.5,:); % v Velocity along centerline
    
    figure(2)
    px = subplot(1,3,1);
    hold on
    plot(Y_midline,p_xmidline,C,'DisplayName',strcat('N =',num2str(N)))
    legend()
    title('Pressure','interpreter','latex','FontSize',lblSz);
    xlabel('$x$','interpreter','latex','FontSize',lblSz);
    ylabel('$p$','interpreter','latex','FontSize',lblSz);

    vortx = subplot(1,3,2);
    hold on
    plot(Y_midline,vort_xmidline,C,'DisplayName',strcat('N = ',num2str(N)))
    legend()
    title('Vorticity','interpreter','latex','FontSize',lblSz);
    xlabel('$x$','interpreter','latex','FontSize',lblSz);
    ylabel('$\omega$','interpreter','latex','FontSize',lblSz);

    vVelx = subplot(1,3,3);
    hold on
    plot(Y_midline,v_xcenterline,C,'DisplayName',strcat('N = ',num2str(N)))
    legend()
    title('Vertical velocity','interpreter','latex','FontSize',lblSz);
    xlabel('$x$','interpreter','latex','FontSize',lblSz);
    ylabel('$v$','interpreter','latex','FontSize',lblSz);

    sgtitle('Pressure, vorticity and vertical velocity along centerline y = 0.5') 
    set(gcf,'Position',[100 100 1500 700]);
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

axes(py)
hold on
plot(Botella_vertdata(:,1),Botella_vertdata(:,4),'--k','DisplayName','Reference results')
axes(vorty)
hold on
plot(Botella_vertdata(:,1),Botella_vertdata(:,5),'--k','DisplayName','Reference results')
axes(uVely)
hold on
plot(Botella_vertdata(:,1),Botella_vertdata(:,3),'--k','DisplayName','Reference results')

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

axes(px)
hold on
plot(Botella_hordata(:,1),Botella_hordata(:,4),'--k','DisplayName','Reference results')
axes(vortx)
hold on
plot(Botella_hordata(:,1),Botella_hordata(:,5),'--k','DisplayName','Reference results')
axes(vVelx)
hold on
plot(Botella_hordata(:,1),Botella_hordata(:,3),'--k','DisplayName','Reference results')

exportgraphics(figure(1),["figures/X05centerline.pdf"], 'Resolution', 300)
exportgraphics(figure(2),["figures/Y05centerline.pdf"], 'Resolution', 300)