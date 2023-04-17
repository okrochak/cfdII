%%
clear all
close all
clc

% The system that you need to solve will be singular. Matlab gives you a
% warning at each time step. To switch of this warning, remove the comment
% in the next line

warning on

% 00D#MMXXI#

% This file contains the skeleton of the program which will solve the lid
% driven cavity problem on a unit square. The parts that have to be
% supplemented are described in the assignment.
%
% The pieces that need to be filled in are indicated
%

%
% When running the code, determine a suitable time step. A too small time
% step will make the calculation very long, while for a too large time step
% the solution will blow up due to numerical instability.
%

Re = 1000;              % Reynolds number
N = 31;      %  N = 15, 31, 47, 55, 63.          % Number of volumes in the x- and y-direction
Delta = 1/N;            % uniform spacing to be used in the mapping to compute tx


% Determine a suitable time step and stopping criterion, tol

tol =   1e-10;             % tol determines when steady state is reached and the program terminates

% wall velocities
U_wall_top = -1;
U_wall_bot = 0;
U_wall_left = 0;
U_wall_right = 0;
V_wall_top = 0;
V_wall_bot = 0;
V_wall_left = 0;
V_wall_right = 0;

%
%   Generation of a non-uniform mesh
%

%
%   tx are the coordinates of the nodal points on the outer-oriented 
%   primal mesh
%

tx = zeros(1,N+1);
for i=1:N+1
    xi = (i-1)*Delta;
    tx(i) = 0.5*(1. - cos(pi*xi));
end

% Mesh width on the outer-oriented primal mesh
th = zeros(N,1);
th = tx(2:N+1) - tx(1:N);

%
%  x are the coordinates of the nodal points on the dual inner-orineted 
%  mesh (including endpoints 0 and 1)
%  h contains the edge lengths on the inner-oriented dual mesh
%
x = 0.5*(tx(1:N) + tx(2:N+1));
x = [0 x 1];

h = zeros(N+1,1);
h = x(2:N+2) - x(1:N+1);

% Stable time step. Note that this is a conservative estimate, so it is
% possible to increase this time step somewhat.

dt = min(min(h),0.5*Re*min(h)^2);
%
%   Initial condition u=v=0
%
%   Both u and v will be stored in one big vector called 'u'
%
%   The vector u only contains the true unknowns, not the velocities that
%   are prescribed by the boundary conditions
%
%   The vector u contains the *inner-oriented* circulations as unknowns

u = zeros(2*N*(N+1),1);


%% Set up the Incindence matrix 'tE21' which connects the fluxes to the volumes.
K = N^2 + 4*N;   % number of rows of Etilde21 incidence matrix with additional cells/fluxes
tE21Ncols = length(u) + 4*N; % number of columns of Etilde21 
tE21indI = zeros(1,4*tE21Ncols); 
tE21indJ = zeros(1,4*tE21Ncols); 
tE21val = zeros(1,4*tE21Ncols); 
counter = 1;

% In those loops I assemble I index, J index and value arrays, which
% will be later assembled in a sparse matrix. This makes the code longer
% but decreases memory consumption. in hindsight i,j approch should have
% been used
for i = 1:(2*N)     % write the first block of rows (green in excel)
    if mod(i,2) == 1
        tE21indI(counter) = i;
        tE21indJ(counter) = i;
        tE21val(counter) = -1;
        counter = counter + 1; % stupid matlab lol
        tE21indI(counter) = i;
        tE21indJ(counter) = 2*N + (floor(i/2)*(N+1)) + 1;
        tE21val(counter) = 1;
        counter = counter + 1; 
    else
        tE21indI(counter) = i; %this is equiv. to         tE21(i,i) = 1;
        tE21indJ(counter) = i;
        tE21val(counter) = 1;
        counter = counter + 1;
        tE21indI(counter) = i;
        tE21indJ(counter) = 3*N + (floor(i/2)-1)*(N+1) + 1;
        tE21val(counter) = -1;        
        counter = counter + 1; 
    end
end
for i = (2*N+1):(2*N+N^2)   % write the middle block (blue one)
    tE21indI(counter) = i;
    tE21indJ(counter) = i + floor((i - 2*N - 1)/N);
    tE21val(counter) = -1;
    counter = counter + 1; 

    tE21indI(counter) = i;
    tE21indJ(counter) = i + floor((i - 2*N - 1)/N) + 1;
    tE21val(counter) =  1;
    counter = counter + 1;  

    tE21indI(counter) = i;
    tE21indJ(counter) =  N*(N+1) + i;
    tE21val(counter) = -1;
    counter = counter + 1; 

    tE21indI(counter) = i;
    tE21indJ(counter) =  N*(N+1) + i + N;
    tE21val(counter) =  1;
    counter = counter + 1; 
end
for i = (2*N+N^2+1):K   % write the last block (orange one)

    if floor((i - (2*N+N^2)-1)/N) == 0 
        tE21indI(counter) = i;
        tE21indJ(counter) = 2*N+N*(N+1) + i - (2*N+N^2);
        tE21val(counter) = 1;
        counter = counter + 1; 
        
        tE21indI(counter) = i;
        tE21indJ(counter) = tE21Ncols - 2*N + i - (2*N+N^2);
        tE21val(counter) = -1;
        counter = counter + 1;
    else
        tE21indI(counter) = i;
        tE21indJ(counter) = length(u) + i - (2*N+N^2);
        tE21val(counter) = -1;
        counter = counter + 1; 

        tE21indI(counter) = i;
        tE21indJ(counter) =  tE21Ncols - 2*N + i - (2*N+N^2);
        tE21val(counter) = 1;
        counter = counter + 1;  
    end
end

% Remove unneeded allocated zeros
tE21indI = nonzeros(tE21indI); 
tE21indJ = nonzeros(tE21indJ); 
tE21val = nonzeros(tE21val); 

tE21 = sparse(tE21indI,tE21indJ,tE21val);
clear tE21val tE21indJ tE21indI counter K 

%% Trim tE21 and create u_norm
% Now we will have to trim tE21. Only the middle block will stay, the other
% will be 'reshuffled'
ind_rem = [1:(2*N) (tE21Ncols-(N*2)+1):tE21Ncols];
M_u_norm = tE21(:,ind_rem);
tE21(:,ind_rem) = []; % the matrix is now trimmed


%  Inserting boundary conditions for normal velocity components and store
%  this part in the vector u_norm, see assignment.

% The problem now became tE21*u+ u_norm = 0; 
% where u_norm = M_u_norm*u_boundary;

% Assemble the normal boundary vector u_Nbc
u_Nbc = zeros(4*N,1);
u_Nbc(1:2:(2*N)) = U_wall_left.*th;
u_Nbc(2:2:(2*N)) = U_wall_right.*th;
u_Nbc((2*N+1):(3*N)) = V_wall_bot.*th;
u_Nbc((3*N+1):(4*N)) = V_wall_top.*th;


% Compute u_norm

u_norm = M_u_norm*u_Nbc;

clear ind_rem
%%  Set up the (extended) sparse, inner-oriented incidence matrix E21
%  The same approach as for tE21 is followed here: first build index
%  arrays, then assemble. 
E21rows = (N+1)^2;
E21indI = zeros(1,4*E21rows); 
E21indJ = zeros(1,4*E21rows); 
E21val = zeros(1,4*E21rows); 
counter = 1;
for i = 1:(N+1) % row index (bottom to top)
    for j = 1:(N+1) %loop cell by cell  %(left to right)
        cellN = (i-1)*(N+1) + j; 
        if i == 1
            E21indI(counter) = cellN;   % first diagonal entries in top left
            E21indJ(counter) = cellN;
            E21val(counter) = 1;
            counter = counter + 1;
        elseif i == (N+1)
            E21indI(counter) = cellN;   % bottom diagonal entries in bottom left
            E21indJ(counter) = j + (N+1);
            E21val(counter) = -1;
            counter = counter + 1;        
        end
        if i ~= (N+1) % the -1 diagonal
            E21indI(counter) = cellN;
            E21indJ(counter) = cellN + 2*(N+1);
            E21val(counter) = -1;
            counter = counter + 1;   
        end
        if i ~= 1 % the 1 diagonal
            E21indI(counter) = cellN;
            E21indJ(counter) = cellN + (N+1);
            E21val(counter) = 1;
            counter = counter + 1;   
        end
        if j == 1
            E21indI(counter) = cellN;   %alone 1s in 1,-1 ; 1,-1 pattern
            E21indJ(counter) = 2*(N+1) + N*(N+1) + (i-1)*N + 1;
            E21val(counter) = 1;
            counter = counter + 1;

            E21indI(counter) = cellN;   %-1s in the rightmost block
            E21indJ(counter) = 2*(N+1) + 2*N*(N+1) + (i-1)*2 + 1;
            E21val(counter) = -1;
            counter = counter + 1;
        elseif j == (N+1)
            E21indI(counter) = cellN;   %alone -1s in 1,-1 ; 1,-1 pattern
            E21indJ(counter) = 2*(N+1) + N*(N+1) + i*N;
            E21val(counter) = -1;
            counter = counter + 1;       

            E21indI(counter) = cellN;    %1s in the rightmost block
            E21indJ(counter) = 2*(N+1) + 2*N*(N+1) + i*2;
            E21val(counter) = 1;
            counter = counter + 1;  
        else        
            E21indI(counter) = cellN;   %pairs of 1 and -1 in 1,-1 ; 1,-1 pattern
            E21indJ(counter) = 2*(N+1) + N*(N+1) + (i-1)*N + (j-1);
            2*(N+1) + N*(N+1) + (i-1)*N + (j-1);
            E21val(counter) = -1;
            counter = counter + 1;  

            E21indI(counter) = cellN;   
            E21indJ(counter) = 2*(N+1) + N*(N+1) + (i-1)*N + (j-1) + 1;
            E21val(counter) = 1;
            counter = counter + 1;              
        end
        
    end
end

E21 = sparse(E21indI,E21indJ,E21val);
clear E21val E21indJ E21indI counter cellN E21rows

%% Trim matrix E21
%  Split off the prescribed tangential velocity and store this in 
%  the vector u_pres
ind_rem = [1:(2*(N+1)) (2*(N+1)+2*N*(N+1)+1):(4*(N+1)+2*N*(N+1))];  %outer parts columns index
M_u_tan = E21(:,ind_rem);
E21(:,ind_rem) = []; % the matrix is now trimmed

% assemble tangential velocity boundary vector u_Tbc (u_prescr)
u_Tbc = zeros(4*(N+1),1);
u_Tbc(1:(N+1)) = U_wall_bot.*h;
u_Tbc((N+1+1):(2*(N+1))) = U_wall_top.*h;
u_Tbc((2*(N+1)+1):2:(4*(N+1))) = V_wall_left.*h;
u_Tbc((2*(N+1))+2:2:(4*(N+1))) = V_wall_right.*h;

u_pres = M_u_tan*u_Tbc;
%%  Set up the sparse, outer-oriented incidence matrix tE10. 

tE10 = E21';    % yet to check for N=3
%%  Set up the sparse, inner-oriented  incidence matrix E10
E10 = -tE21';   % checked for N=3 

%%  Set up the Hodge matrices Ht11 and H1t1

%H1t1 and H1t1 will have size of len(u) x len(u)

% build the diagonal terms

ind_th = zeros(1,length(u));    % which th should this edge take
ind_h = zeros(1,length(u));      % which h should this edge take

% first write the u-part, first half 
for i = 1:N     
    for j = 1:(N+1)
        k = (i-1) *(N+1) + j;
        % (i,j) corresponds to primal mesh
        ind_th(k) = i;
        ind_h(k) = j;
    end
end
% now the v-part, second half
for i = 1:(N+1)     
    for j = 1:N
        k = N*(N+1) + (i-1) *(N) + j;
        % (i,j) corresponds to the primal mesh. Now the sing is flipped
        % though. 
        ind_th(k) = j;
        ind_h(k) = i;
    end
end
% Now with indices known, setting up the Hodge diagonals is trivial
H1t1diag = (h(ind_h)./th(ind_th));
Ht11diag = (th(ind_th)./h(ind_h));
H1t1 = spdiags(H1t1diag',0,length(u),length(u));
Ht11 = spdiags(Ht11diag',0,length(u),length(u));

clear H1t1diag Ht11diag ind_h ind_rem ind_th 
%  Set up the Hodge matrix Ht02

%%
ind_h = zeros(1,length(h));    % which width h should this edge take
ind_h_2 = zeros(1,length(h));      % which height h should this edge take

%  Determining indexes row by row
for j = 1:(N+1)     
    for i = 1:(N+1)
        k = (i-1) *(N+1) + j;
        % (i,j) corresponds to primal mesh
        ind_h(k) = i;
        ind_h_2(k) = j;
    end
end

% Now with indices known, setting up the Hodge diagonals is trivial
Ht02diag = (1./(h(ind_h).*h(ind_h_2)));
Ht02 = spdiags(Ht02diag',0,length(h)^2,length(h)^2);

%
% The prescribed velocties will play a role in the momentum equation
%
u_pres_vort=Ht02*u_pres; %U_pres to outer oriented 0 form representing contribution of boundary conditions to point wise vorticity
u_pres = H1t1*E21'*Ht02*u_pres; %U_pres to inner oriented 1 forms


% Now all matrices are set up and the time stepping can start. 'iter' will
% record the number of time steps. This allows you to give output after a
% preselected number of time steps.
%
% 'diff' will be the maximal du/dt or dv/dt. If 'diff' is sufficiently
% small, steady state has been reached. Determine a suitable value for
% 'tol'
%

convective = zeros(2*N*(N+1),1);
ux_xi = zeros((N+1)*(N+1),1);
uy_xi = zeros((N+1)*(N+1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    From here the code enters the time-stepping loop and no new parts in
%    the code need to be inserted. If done correctly, everything should
%    work now.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diff = 1;
iter = 1;


% Set up the matrix for the Poisson equation    

A = -tE21*Ht11*tE21'; 
if sum(sum(A,2)) > N*(10^(-15))
    warning("Poisson matrix row-sum too large!");
end
% Perform an LU-decomposition for the pressure matrix A

[L,U] = lu(A);

% Abbreviation for some matrix products which are constant in the time loop

VLaplace = H1t1*E21'*Ht02*E21;
DIV = tE21*Ht11;

while diff > tol
%while iter < 10
    %Vector xi is obtained. It corresponds with the point-wise vorticity
    %at each cell
    
    %Vectors ux_xi and uy_xi correspond with the multiplications of
    %xi with the horizontal and vertical velocity components at each cell.
    %Only the cells required for vector convective are calculated. The
    %ordering of each vector with respect to the ordering of cells in the
    %grid is different (left to right for ux_xi and bottom to top for
    %uy_xi)
    
    xi = Ht02*E21*u + u_pres_vort;
    
    for i=1:N+1
        for j=1:N+1
            k = j + (i-1)*(N+1); 
            if j==1
                ux_xi(k) = U_wall_bot*xi(i+(j-1)*(N+1));    
                uy_xi(k) = V_wall_left*xi(j+(i-1)*(N+1));
            elseif j==N+1
                ux_xi(k) = U_wall_top*xi(i+(j-1)*(N+1));
                uy_xi(k) = V_wall_right*xi(j+(i-1)*(N+1));
            else
                ux_xi(k) = (u(i+(j-1)*(N+1))+u(i+(j-2)*(N+1)))*xi(i+(j-1)*(N+1))/(2.*h(i));                      
                uy_xi(k) = (u(N*(N+1)+j+(i-1)*N) + u(N*(N+1)+j-1+(i-1)*N))*xi(j+(i-1)*(N+1))/(2.*h(i));  
            end
        end
    end

    for  i=1:N
        for j=1:N+1
            convective(j+(i-1)*(N+1)) = -(uy_xi(j+(i-1)*(N+1))+uy_xi(j+i*(N+1)))*h(j)/2.;
            convective(N*(N+1)+i+(j-1)*N) = (ux_xi(j+(i-1)*(N+1))+ux_xi(j+i*(N+1)))*h(j)/2.;
        end
    end
    

    % Set up the right hand side for the Poisson equation for the pressure
    
    rhs_Poisson  =   DIV*(u/dt  - convective - VLaplace*u/Re - u_pres/Re) + u_norm/dt; 
    
    % Solve for the new pressure
    
    temp = L\rhs_Poisson;
    p = U\temp;
    
    % Store the velocity from the previous time step in the vector u_old
    
    uold = u;
    
    % Update the velocity field
    
    u = u - dt* (convective - tE21'*p + VLaplace*u/Re + u_pres/Re); 
    
    %
    %  Every other 1000 iterations check whether you approach steady state
    %  and check whether yopu satisfy conservation of mass. The largest
    %  rate at which mass is destroyed or created is denoted by 'maxdiv'.
    %  This number should be very small, in the order of machine precision.
    
    if mod(iter,1000) == 0
    
        maxdiv = max(DIV*u + u_norm) 
        
        diff = max(abs(u-uold))/dt
        
    end
    iter = iter + 1;
end

%% Save results
clear U

filename = "results/results_N="+N+"_tol="+string(tol)+"_Re="+Re+".mat"; %filename to save workspace to for post-processing
save(filename)

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

    % Calculating velocities at dual edges
h_vec = [repmat((1./h),1,N), repelem((1./h), N)]; % Vector containing edge widths
Uedge = U.*h_vec'; % Converting fluxes to velocities along dual edges
    
    % Calculating velocities at dual points in x- direction
Uedge_x = Uedge(1:size(Uedge)/2); % Edge velocities along x- direction
Upoint_x = (Uedge_x(1:size(Uedge_x)-1) + Uedge_x(2:end)) / 2; % Averaging over edge velocities
Upoint_x(N+1 : N+1 : size(Upoint_x,1)) = []; % Eliminating unphysical averages
Upoint_xMat = reshape(Upoint_x,N,N); % Dual point x- velocity matrix
    
    % Calculating velocities at dual points in y- direction
Uedge_y = Uedge(size(Uedge)/2+1:end); % Edge velocities along y- direction
UedgeMat_y = reshape(Uedge_y,N,N+1); % Edge velocity matrix
Upoint_yMat = zeros (N,N); % Dual point y- velocity matrix
for i = 1:N
    for j = 1:N
        Upoint_yMat(i,j) = ( UedgeMat_y(i,j+1) + UedgeMat_y(i,j) )/ 2;
    end
end
Upoint_y = reshape(Upoint_yMat,N*N,1);

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

pMat = reshape(p_trim,N,N); % Reshape vector into pressure matrix
pMat = pMat - 0.5*(Upoint_xMat.^2 + Upoint_yMat.^2); % Substract dynamic pressure
p_ref = pMat((N/2+0.5),(N/2+0.5)); % Reference pressure
pMat = pMat - p_ref; % Adjusting pressure

% p_ref = p_trim((N*N/2+0.5)); % Reference pressure
% pdata = p_trim - p_ref; % Adjusting pressure
% pMat = reshape(pdata,N,N); % Reshape vector into pressure matrix
% pMat = pMat - 0.5*(Upoint_xMat.^2 + Upoint_yMat.^2); % Substract dynamic pressure

%%% Plotting routine
%% Plot Vorticity
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

contourf(tX,tY,newdata,tickVals)
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
p = pcolor(tX,tY,tX*0);
alpha(p,0); p.EdgeAlpha = 0.1; hold off;
exportgraphics(gcf,["figures/p_N"+string(N)+".pdf"], 'Resolution', 300)

%% Plot components along x = 0.5, y = 0.5;
% use linear interpolation for intermediate fluxes. maybe use circulation
% isntead?
% xline.tu = tuMat(ceil(N/2),:);
% xline.tv = (tvMat(ceil(N/2),:) + tvMat(ceil(N/2)+1,:)) / 2;+
% yline.tu = (tuMat(:,ceil(N/2)) + tuMat(:,ceil(N/2)+1)) / 2;
% yline.tv = tvMat(:,ceil(N/2))
%contourf(xp(2:end-1),xp(2:end-1),psi,-0.5:0.05:0.5);

    % Variables along line x = 0.5
Y_midline = xd;
p_ymidline = pMat((N/2+0.5),:);
vort_ymidline = vortMat(:,(N/2+0.5));
u_ycenterline =  Upoint_xMat(N/2+0.5,:);
v_ycenterline =  Upoint_yMat(N/2+0.5,:);
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

figure()
subplot(1,3,1)
hold on
plot(Y_midline,p_ymidline,'k')
plot(Botella_vertdata(:,1),Botella_vertdata(:,4),'--k')
legend('Calculation results','Reference results')
title('Pressure','interpreter','latex','FontSize',lblSz);
xlabel('$y$','interpreter','latex','FontSize',lblSz);
ylabel('$p$','interpreter','latex','FontSize',lblSz);
hold off

subplot(1,3,2)
hold on
plot(Y_midline,vort_ymidline(1:end-1),'k')
plot(Botella_vertdata(:,1),Botella_vertdata(:,5),'--k')
legend('Calculation results','Reference results')
title('Vorticity','interpreter','latex','FontSize',lblSz);
xlabel('$y$','interpreter','latex','FontSize',lblSz);
ylabel('$\omega$','interpreter','latex','FontSize',lblSz);
hold off

subplot(1,3,3)
hold on
plot(Y_midline,u_ycenterline,'k')
plot(Botella_vertdata(:,1),Botella_vertdata(:,3),'--k')
legend('Calculation results','Reference results')
title('Velocity','interpreter','latex','FontSize',lblSz);
xlabel('$y$','interpreter','latex','FontSize',lblSz);
ylabel('$u$','interpreter','latex','FontSize',lblSz);
hold off
sgtitle('Pressure, vorticity and horizontal velocity along centerline x = 0.5') 
set(gcf,'Position',[100 100 1500 700])

    % Variables along line y = 0.5
X_midline = xd;
p_xmidline = pMat(:,(N/2+0.5));
vort_xmidline = vortMat((N/2+0.5),:);
u_xcenterline =  Upoint_xMat(:,N/2+0.5);
v_xcenterline =  Upoint_yMat(:,N/2+0.5);
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
figure()

subplot(1,3,1)
hold on
plot(Y_midline,p_xmidline,'r')
plot(Botella_hordata(:,1),Botella_hordata(:,4),'--r')
legend('Calculation results','Reference results')
title('Pressure','interpreter','latex','FontSize',lblSz);
xlabel('$x$','interpreter','latex','FontSize',lblSz);
ylabel('$p$','interpreter','latex','FontSize',lblSz);
hold off

subplot(1,3,2)
hold on
plot(Y_midline,vort_xmidline(1:end-1),'r')
plot(Botella_hordata(:,1),Botella_hordata(:,5),'--r')
legend('Calculation results','Reference results')
title('Vorticity','interpreter','latex','FontSize',lblSz);
xlabel('$x$','interpreter','latex','FontSize',lblSz);
ylabel('$\omega$','interpreter','latex','FontSize',lblSz);
hold off

subplot(1,3,3)
hold on
plot(Y_midline,v_xcenterline,'r')
plot(Botella_hordata(:,1),Botella_hordata(:,3),'--r')
legend('Calculation results','Reference results')
title('Vertical velocity','interpreter','latex','FontSize',lblSz);
xlabel('$x$','interpreter','latex','FontSize',lblSz);
ylabel('$v$','interpreter','latex','FontSize',lblSz);
hold off
sgtitle('Pressure, vorticity and vertical velocity along centerline y = 0.5') 
set(gcf,'Position',[100 100 1500 700])
