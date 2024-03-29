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
Re = 10000;   % Reynolds number
N = 63  ;      % N = 15, 31, 47, 55, 63.          % Number of volumes in the x- and y-direction
dt_mult = 5; % Tested to work at 5 for all N, can maybe go higher? Tested at 10 and it was unstable already
Delta = 1/N; % uniform spacing to be used in the mapping to compute tx

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

dt = min(min(h),0.5*Re*min(h)^2)*dt_mult;
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

% In those loops we assemble I index, J index and value arrays, which
% will be later assembled in a sparse matrix. This makes the code longer
% but decreases memory consumption. in hindsight i,j approch should have
% been used

for i = 1:(2*N)     % write the first block of rows (green in excel)
    if mod(i,2) == 1
        tE21indI(counter) = i;
        tE21indJ(counter) = i;
        tE21val(counter) = -1;
        counter = counter + 1;
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
for i = (2*N+N^2+1):K   % Write the last block (orange one)

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

% Inserting boundary conditions for normal velocity components and store
% this part in the vector u_norm, see assignment.
% The problem now became tE21*u+ u_norm = 0; 
% where u_norm = M_u_norm*u_boundary;

% Assemble the normal boundary vector u_Nbc
u_Nbc = zeros(4*N,1);
u_Nbc(1:2:(2*N)) = U_wall_left.*th; % Left wall BC
u_Nbc(2:2:(2*N)) = U_wall_right.*th; % Bottom wall BC
u_Nbc((2*N+1):(3*N)) = V_wall_bot.*th; % Bottom wall BC
u_Nbc((3*N+1):(4*N)) = V_wall_top.*th; % Upper BC

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
u_Tbc(1:(N+1)) = U_wall_bot.*h; % Lower wall BC
u_Tbc((N+1+1):(2*(N+1))) = U_wall_top.*h; % Upper BC
u_Tbc((2*(N+1)+1):2:(4*(N+1))) = V_wall_left.*h; % Left wall BC
u_Tbc((2*(N+1))+2:2:(4*(N+1))) = V_wall_right.*h; % Right wall BC

u_pres = M_u_tan*u_Tbc;

%%  Set up the sparse, outer-oriented incidence matrix tE10. 
tE10 = E21';    % tE10 is related to E21 by the transpose of E21

%%  Set up the sparse, inner-oriented  incidence matrix E10
E10 = -tE21';   % E10 is related to tE21 by the transpose of -tE21

%%  Set up the Hodge matrices Ht11 and H1t1
% H1t1 and H1t1 will have size of len(u) x len(u)
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
        % (i,j) corresponds to the primal mesh. Now the indexing is flipped
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
%%  Set up the Hodge matrix Ht02
% build the diagonal terms
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
Ht02diag = (1./(h(ind_h).*h(ind_h_2))); % 1/ Area of dual cell
Ht02 = spdiags(Ht02diag',0,length(h)^2,length(h)^2);

%%
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
A_sum = sum(sum(A,2));
A_sym = issymmetric(A);
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

filename = "results/results_N="+N+"_tol="+string(tol)+"_Re="+Re+"_dtmult="+dt_mult+".mat"; %filename to save workspace to for post-processing
save(filename)