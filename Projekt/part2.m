% Load matrices p, e, t from PDETOOL mesh save

% load('pet1grov.mat')
load('pet2med.mat')

% PREPROCESSOR

% Symmetry line at x = 0.

% SUBDOMAINS
%
% 1: Silver-epoxy
% 2 & 3: Copper
% 4: Silicon

% -- Define constants --
% The order in the arrays are defined by which constant
% corresponds to which subdivision in the mesh.

% subd. 2 & 3 are the same material

% Conductivity constants

load('constants.mat');

Q = 5 * 10^7; % W / m^3

% BOUNDARY VALUES

% Bottom boundaries are fixated in x- y direction.
% u_x = 0, u_y = 0.
% Needed for later part of project
% Symmetry line at x = 0 is fixated to other side, i.e
% u_x = 0.

% Plane strain conditions

% Newton convection along top boundary
% q_n = alpha_c* (T - T_inf), alpha_c = 40 W / (m^2 K).
alpha_c = 40;

% Rest of boundaries insulated.
% Initial temperature T_0 = 30 C.


% Calculate EDOF matrix from pdetool arrays
% Assumes arrays p, e, t are given from pdetool

coord=p';

coord = coord * 10^(-3);

enod=t(1:3,:)'; % nodes of elements
nelm=size(enod,1); % number of elements
nnod=size(coord,1); % number of nodes
dof=(1:nnod)'; % dof number is node number
dof_S=[(1:nnod)',(nnod+1:2*nnod)']; % give each dof a number

for ie=1:nelm
    edof_S(ie,:)=[ie dof_S(enod(ie,1),:), dof_S(enod(ie,2),:),dof_S(enod(ie,3),:)];
    edof(ie,:)=[ie,enod(ie,:)];
end

nen = 3; % Nbr of nodes per element

% edof = zeros(nelm, nen + 1);

% Extract pdetool mesh into CALFEM format

edof(:,1)=1:nelm ;
edof(:,2:4)=t(1:3,:)';

ndof=max(max(t(1:3,:)));
[Ex,Ey]=coordxtr(edof,coord,(1:ndof)',3);
% eldraw2(Ex,Ey,[1,4,1])


% Obtain a list of all CONVECTIVE BOUNDARIES --------------

% Check which segments should have convections

er = e([1 2 5],:); % Reduced e
%conv_segments = [10 11 12]; % Choosen boundary segments
conv_segments = [14 1 17]; % According to edge numbering
edges_conv = [];
for i = 1:size(er,2)
    if ismember(er(3,i),conv_segments)
        edges_conv = [edges_conv er(1:2,i)];
    end
end

% ---------------------------------------------------------

% Build K matrix from Kc and Ke.

K = zeros(ndof);
fb = zeros(ndof, 1);
C = zeros(ndof);

% Calc Kc for each edge and add to K

for i = 1:length(edges_conv)
    n1 = edges_conv(1, i); % Node 1
    n2 = edges_conv(2, i); % Node 2
    
    p1 = coord(n1, :); % coordinate for point 1
    p2 = coord(n2, :); % for point 2
    
    L = dist(p1,p2'); % Distance
    
    Kce = thickness * alpha_c * (L/6) * [2, 1; 1, 2];
    fbe = T_inf * alpha_c * thickness * [L/2 ; L/2];
    
    % K = assem(edof(i,:), K, Kc);
    
    % Kt = assem(edof(el,:),Kt,Kte);
    
    % Insert Kc at correct position
    K([n1,n2],[n1,n2]) = K([n1,n2],[n1,n2]) + Kce;
    
    fb([n1, n2]) = fb([n1, n2]) + fbe;
    
end

Q_vec = [0, 0 ,0 , Q]; % Only in the silicon part, subd 4

fl = zeros(ndof, 1);

for elnr = 1:nelm
    % Get correct k constant-matrix
    sd = t(4,elnr); % subdomain
    D = eye(2)*k_const(sd); % cond. matrix
    eq = Q_vec(sd); % heat source
    
    [Ke, fe] = flw2te(Ex(elnr, :), Ey(elnr, :), thickness, D, eq);
    
    % assemble K
    [K, fl]= assem(edof(elnr,:), K, Ke, fl, fe);
    
    % Build C matrix (time stepping)
    temp_const = ro_const(sd)*c_const(sd)*thickness;
    
    Ce = plantml(Ex(elnr,:), Ey(elnr,:), temp_const);
    
    indx = edof(elnr, 2:end);
    C(indx, indx) = C(indx, indx) + Ce;
    
end

f = fb + fl;

% ---- Time stepping ----

dt = 1; % Size of time step
n_steps = 0;
t_length = 60*20; % 20 minutes

a_old = ones(ndof, 1) * T_0; % Initial temperatures

a_next = (C + dt*K) \ (C*a_old + dt * f);

n_steps = n_steps + 1;

% Create array for max temperature
temp_max = zeros(t_length, 1);
temp_max(1) = T_0;
temp_max(2) = max(a_next);

while n_steps < t_length
    
    a_old = a_next;
    
    a_next = (C + dt*K) \ (C*a_old + dt * f);
    
    n_steps = n_steps + 1;
    
    % add max temp to array
    temp_max(n_steps) = max(a_next);
    
end

%% ----- PLOT ------

eT = extract(edof, a_next);
patch(Ex', Ey',eT')
hold on
patch(-Ex', Ey',eT')

colorbar
colormap('hot')
