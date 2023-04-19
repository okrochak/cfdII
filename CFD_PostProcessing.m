%%% Post-processing
close all
% Plotting preferences
fntSz = 20;
lblSz = 25;
lnWd = 2;

%% Calculating domain properties
    
% Domain velocities
U = u; % Flux vector over dual edges
tU = Ht11*U; % Flux vector over primal edges
[X,Y] = meshgrid(x(2:end-1)); % coordinates of primal cell centers

% Vorticity at dual points
vort = Ht02*E21*H1t1*tU; % Convert fluxes to vorticity
vortMat = reshape(vort,(N+1),(N+1))'; % Reshape vorticity vector into matrix
    
% Calculating the average velocity over primal edges
tuSize = length(tU);
tv = tU((tuSize/2)+1:end); % x and y velocity fluxes on primal mesh
tu = tU(1:(tuSize/2));
    
% tuMat, tvMat are fluxes divided by cell width, velocities through primal edges
tuMat = reshape(tu,N+1,N)' ./ repmat(th',1,N+1); % divide the fluxes by the mesh width
tvMat = reshape(tv,N,N+1)' ./ repmat(th,N+1,1); % divide the fluxes by the mesh width
        
% Calculate cell average values for primal cells
tuIntrp = interp2(tx,x(2:end-1),tuMat,X,Y);
tvIntrp = interp2(x(2:end-1),tx,tvMat,X,Y);
vortIntrp = interp2(tx,tx,vortMat,X,Y);

% Integrate velocity to get to the streamfunction (psi = 0 on y = 0)
psi = cumtrapz(x(2:end-1),tuIntrp,1);

% Integrated vorticity over domain
vortInt = trapz(tx,trapz(tx,vortMat,1),2);

% Inner dual point pressures 
p_trim = p([(2*N+1):(size(p)-2*N)]); % Pressure at inner dual points
pMat = reshape(p_trim,N,N)'; % Reshape vector into pressure matrix
pMat = pMat - 0.5*(tuIntrp.^2 + tvIntrp.^2); % Substract dynamic pressure
p_ref = pMat((N/2+0.5),(N/2+0.5)); % Reference pressure
pMat = pMat - p_ref; % Adjusting pressure

%% Integrated vorticity?
vortInt = trapz(tx,trapz(tx,vortMat,1),2);

%% Assemble post-processing results in a single structure
postProc.p = pMat;
postProc.vort = vortIntrp;
postProc.u = tuIntrp;
postProc.v = tvIntrp;
postProc.x = x(2:end-1);

save(strcat("results/PostProcessingN_",string(N)),"postProc")
%%% Plotting routine
%% Vorticity
figure(1)
set(gcf,'Position',[100 100 800 700])

contLvl = [-3 -2 -1 -0.5 0 0.5 1 2 3 4 5];
tickVals = linspace(contLvl(1), contLvl(end), length(contLvl));
% adapted from https://nl.mathworks.com/matlabcentral/answers/738107-how-to-make-a-discrete-colorbar-with-specified-intervals?s_tid=prof_contriblnk
data = vortMat; newdata = data; 
for tv = 2:length(contLvl)
    % find all data between the new tick ranges, one block at a time
    ind = data>contLvl(tv-1) & data<=contLvl(tv);
    % Change the corresponding values in the copied data to scale to the
    % block edges
    newdata(ind) = rescale(data(ind),tickVals(tv-1),tickVals(tv));
end
    
contourf(tx,tx,newdata,tickVals)
C=turbo(length(contLvl));
tickVals = linspace(contLvl(1), contLvl(end), length(contLvl)+1);
colormap(flipud(C))
colorbar('Ticks',tickVals,...
    'TickLabels',[string(contLvl(1:end)) ""],'TickLabelInterpreter','latex');
title('$\omega $ for $N = $'+string(N),'interpreter','latex','FontSize',lblSz);
xlabel('$x$','interpreter','latex','FontSize',lblSz);
ylabel('$y$','interpreter','latex','FontSize',lblSz);
set(gca,'FontSize', fntSz,'TickLabelInterpreter','latex');
hold on
p = pcolor(X,Y,X*0);
alpha(p,0); p.EdgeAlpha = 0.1; hold off;
exportgraphics(gcf,["figures/omega_N"+string(N)+".pdf"], 'Resolution', 300)



%% Plot Stream function
figure(2)
set(gcf,'Position',[100 100 800 700])

contLvl = flip([0.1175 0.115 0.11 0.1 0.09 0.07 0.05 0.03 0.01 0.001 1e-5 1e-10 ...
    0 -1e-6 -1e-5 -5e-5 -1e-4 -2.5e-4 -5e-4 -1e-3 -1.5e-3]);
tickVals = linspace(contLvl(1), contLvl(end), length(contLvl));
% adapted from https://nl.mathworks.com/matlabcentral/answers/738107-how-to-make-a-discrete-colorbar-with-specified-intervals?s_tid=prof_contriblnk
data = psi; newdata = data; 
for tv = 2:length(contLvl)
    % find all data between the new tick ranges, one block at a time
    ind = data>contLvl(tv-1) & data<=contLvl(tv);
    % Change the corresponding values in the copied data to scale to the
    % block edges
    newdata(ind) = rescale(data(ind),tickVals(tv-1),tickVals(tv));
end

contour(x(2:end-1),x(2:end-1),newdata,tickVals)
C=turbo(length(contLvl));
tickVals = linspace(contLvl(1), contLvl(end), length(contLvl)+1);
colormap(flipud(C))
cmTicks = [string(contLvl(1:end)) ""];
cmTicks(1:2:end) = "";
colorbar('Ticks',tickVals,...
    'TickLabels',cmTicks,'TickLabelInterpreter','latex');
title('Streamlines for $N = $'+string(N),'interpreter','latex','FontSize',lblSz);
xlabel('$x$','interpreter','latex','FontSize',lblSz);
ylabel('$y$','interpreter','latex','FontSize',lblSz);
set(gca,'FontSize', fntSz,'TickLabelInterpreter','latex');
hold on
p = pcolor(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),X(2:end-1,2:end-1)*0);
alpha(p,0); p.EdgeAlpha = 0.1; hold off;
exportgraphics(gcf,["figures/psi_N"+string(N)+".pdf"], 'Resolution', 300)

%% Plot pressure field;
figure(3)
set(gcf,'Position',[100 100 800 700])

contLvl = [-0.002 0.0 0.02 0.05 0.07 0.09 0.11 0.12 0.17 0.3];
tickVals = linspace(contLvl(1), contLvl(end), length(contLvl));
% adapted from https://nl.mathworks.com/matlabcentral/answers/738107-how-to-make-a-discrete-colorbar-with-specified-intervals?s_tid=prof_contriblnk
data = pMat; newdata = data; 
for tv = 2:length(contLvl)
    % find all data between the new tick ranges, one block at a time
    ind = data>contLvl(tv-1) & data<=contLvl(tv);
    % Change the corresponding values in the copied data to scale to the
    % block edges
    newdata(ind) = rescale(data(ind),tickVals(tv-1),tickVals(tv));
end

contourf(X,Y,newdata,tickVals)
C=turbo(length(contLvl));
tickVals = linspace(contLvl(1), contLvl(end), length(contLvl)+1);
colormap(flipud(C))
colorbar('Ticks',tickVals,...
    'TickLabels',[string(contLvl(1:end)) ""],'TickLabelInterpreter','latex');
title('Isobaric lines for $N = $'+string(N),'interpreter','latex','FontSize',lblSz);
xlabel('$x$','interpreter','latex','FontSize',lblSz);
ylabel('$y$','interpreter','latex','FontSize',lblSz);
set(gca,'FontSize', fntSz,'TickLabelInterpreter','latex');
hold on
p = pcolor(X,Y,X*0);
alpha(p,0); p.EdgeAlpha = 0.1; hold off;
exportgraphics(gcf,["figures/p_N"+string(N)+".pdf"], 'Resolution', 300)