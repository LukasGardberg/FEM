%% Solution to PART 1, HEAT FLOW %%

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

% Conductivity constants

k_epoxy = 5;
k_silicon = 149;
k_copper = 385;

k_const = [k_epoxy, k_copper, k_copper, k_silicon];

% Thickness: 10mm
thickness = 10 * 10^(-3); % in meters

Q = 5 * 10^7; % W / m^3
%Q = 5 * 10^7 * 0.75;

% Surrounding temperature
T_inf = 18;

% BOUNDARY VALUES

% Bottom boundaries are fixated in x- y direction.
% u_x = 0, u_y = 0.
% Symmetry line at x = 0 is fixated to other side, i.e
% u_x = 0.

% Plane strain conditions

% Newton convection along top boundary
% q_n = alpha_c* (T - T_inf), alpha_c = 40 W / (m^2 K).
alpha_c = 40;

% Rest of boundaries insulated.

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

% Extract pdetool mesh into CALFEM format

edof(:,1)=1:nelm ;
edof(:,2:4)=t(1:3,:)';

ndof=max(max(t(1:3,:)));
[Ex,Ey]=coordxtr(edof,coord,(1:ndof)',3);
% eldraw2(Ex,Ey,[1,4,1])


% Obtain a list of all CONVECTIVE BOUNDARIES --------------

% Check which segments that should have convections

er = e([1 2 5],:); % Reduced e
%conv_segments = [10 11 12]; % Choosen boundary segments
conv_segments = [14 1 17]; % Corresponds to our specific mesh
edges_conv = [];
for i = 1:size(er,2)
    if ismember(er(3,i),conv_segments)
        edges_conv = [edges_conv er(1:2,i)];
    end
end

% ---------------------------------------------------------

% Initialize matrices
K = zeros(ndof);
fb = zeros(ndof, 1);

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

% Global load vector
fl = zeros(ndof, 1);

for elnr = 1:nelm
    % Get correct k constant-matrix
    sd = t(4,elnr); % subdomain
    D = eye(2)*k_const(sd); % cond. matrix
    eq = Q_vec(sd); % heat source
    
    [Ke, fe] = flw2te(Ex(elnr, :), Ey(elnr, :), thickness, D, eq);
    
    [K, fl]= assem(edof(elnr,:), K, Ke, fl, fe);
    
end

f = fb + fl;

a = solveq(K, f);

%% ----- PLOT ------

eT = extract(edof, a);


figure()
patch(Ex', Ey',eT')
hold on
patch(-Ex', Ey',eT')
title('Temperature distribution [C]')
colormap(hot);
colorbar;
xlabel('x-position [m]')
ylabel('y-position [m]')
axis equal

%Tmin = 66.4;
%Tmin = 54.3;
%Tmax = 68.3;
%caxis([Tmin Tmax])


