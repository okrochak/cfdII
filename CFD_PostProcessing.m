%%% Post-processing
% Plotting preferences
fntSz = 20;
lblSz = 25;
lnWd = 2;

%% Calculating domain properties
U = u; 
%clear u; %U = full U vector on dual mesh
    
    % Contructing the primal mesh grid
xp = tx; 
for i = 1:N
    yp(i) = 0.5*(xp(i)+xp(i+1));
end
X = xp(2:end-1);
[X,Y] = meshgrid(X); % primal grid

    % Constructing the inner dual mesh grid
xd = x([(2):(N+1)]);
[tY,tX] = meshgrid(xd,xd); % inner dual grid

    % Velocity fluxes over primal edges vector
tU = Ht11*U; % Convert to fluxes on the primal mesh.

    % Vorticity at dual points
vort = Ht02*E21*H1t1*tU; % Convert fluxes to vorticity
vortMat = reshape(vort,(N+1),(N+1))'; % Reshape vorticity vector into matrix
    
    % Calculating the average velocity over primal edges
tuSize = length(tU);
tv = tU((tuSize/2)+1:end); % x and y velocity fluxes on primal mesh
tu = tU(1:(tuSize/2));
tuMat = reshape(tu,N+1,N)' ./ repmat(th',1,N+1); % divide the fluxes by the mesh width
tvMat = reshape(tv,N,N+1)' ./ repmat(th,N+1,1);
        
    % Interpolated primal velocity field
tuIntrp = interp2(xp,yp,tuMat,X,Y);
tvIntrp = interp2(yp,xp,tvMat,X,Y);

    % Integrate velocity to get to the streamfunction
[X,Y] = meshgrid(xp,xp);
psi = cumtrapz(xp(2:end-1),tuIntrp,1) - cumtrapz(xp(2:end-1),tvIntrp,2);

    % Integrated vorticity over domain
vortInt = trapz(xp,trapz(xp,vortMat,1),2);

    % Interpolated Dual point velocity
uIntrp = interp2(xp,yp,tuMat,tX,tY); % Velocity in x-direction
vIntrp = interp2(yp,xp,tvMat,tX,tY); % Velocity in y-direction

    % Inner dual point pressures 
p_trim = p([(2*N+1):(size(p)-2*N)]); % Pressure at inner dual points
p_ref = p_trim((N*N/2+0.5)); % Reference pressure
pdata = p_trim - p_ref; % Adjusting pressure
pMat = reshape(pdata,N,N); % Reshape vector into pressure matrix
pMat = pMat - 0.5*(uIntrp.^2 + vIntrp.^2); % Substract dynamic pressure


%%% Plotting routine
%% Vorticity plot
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
    
contourf(xp,xp,newdata,tickVals)
C=turbo(length(contLvl));
tickVals = linspace(contLvl(1), contLvl(end), length(contLvl)+1);
colormap(flipud(C))
colorbar('Ticks',tickVals,...
    'TickLabels',[string(contLvl(1:end)) ""],'TickLabelInterpreter','latex');
title('$\omega $ for $N = $'+string(N),'interpreter','latex','FontSize',lblSz);
xlabel('$x$','interpreter','latex','FontSize',lblSz);
ylabel('$y$','interpreter','latex','FontSize',lblSz);
set(gca,'FontSize', fntSz,'TickLabelInterpreter','latex');
exportgraphics(gcf,["figures/omega_N"+string(N)+".pdf"], 'Resolution', 300)



%% Stream function
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

contour(xp(2:end-1),xp(2:end-1),newdata,tickVals)
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
exportgraphics(gcf,["figures/psi_N"+string(N)+".pdf"], 'Resolution', 300)


%% Plot components along x = 0.5, y = 0.5;
% use linear interpolation for intermediate fluxes. maybe use circulation
% isntead?
figure(3)
xline.tu = tuMat(ceil(N/2),:);
xline.tv = (tvMat(ceil(N/2),:) + tvMat(ceil(N/2)+1,:)) / 2;

yline.tu = (tuMat(:,ceil(N/2)) + tuMat(:,ceil(N/2)+1)) / 2;
yline.tv = tvMat(:,ceil(N/2));

contourf(xp(2:end-1),xp(2:end-1),psi,-0.5:0.05:0.5);

%% Plot pressure field;
figure(4)
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

contourf(tX,tY,pMat,tickVals)
C=turbo(length(contLvl));
tickVals = linspace(contLvl(1), contLvl(end), length(contLvl)+1);
colormap(flipud(C))
colorbar('Ticks',tickVals,...
    'TickLabels',[string(contLvl(1:end)) ""],'TickLabelInterpreter','latex');
title('Isobaric lines for $N = $'+string(N),'interpreter','latex','FontSize',lblSz);
xlabel('$x$','interpreter','latex','FontSize',lblSz);
ylabel('$y$','interpreter','latex','FontSize',lblSz);
set(gca,'FontSize', fntSz,'TickLabelInterpreter','latex');
exportgraphics(gcf,["figures/p_N"+string(N)+".pdf"], 'Resolution', 300)
