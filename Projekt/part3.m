% Solves the last subproblem for the project in FHLF10
%
% a) Stationary heat distributin
%
% b) transient heat distribution
% 
% - c) the von Mises stress field and the displacement
% field due to thermal expansion

% ---- PREPROCESSOR ----

% Load matrices p, e, t from PDETOOL mesh save

%load('pet1grov.mat')
load('pet2med.mat')

% Symmetry line at x = 0.

% SUBDOMAINS
%
% 1: Silver-epoxy
% 2 & 3: Copper
% 4: Silicon

% Load constants

load('constants.mat');

Q = 5 * 10^7; % W / m^3

% BOUNDARY VALUES

% Bottom boundaries are fixated in x- y direction.
% u_x = 0, u_y = 0.

% Plane strain conditions

% Newton convection along top boundary
% q_n = alpha_c* (T - T_inf), alpha_c = 40 W / (m^2 K).
alpha_c = 40;

% Rest of boundaries insulated.
% Initial temperature T_0 = 30 C.

% Calculate EDOF matrix from pdetool arrays p, e, t.

coord=p';

% Change to mm
coord = coord * 10^(-3);

enod=t(1:3,:)'; % nodes of elements
nelm=size(enod,1); % number of elements
nnod=size(coord,1); % number of nodes
dof=(1:nnod)'; % dof number is node number
dof_S=[(1:nnod)',(nnod+1:2*nnod)']; % give each dof a number

for ie=1:nelm
    % edof_S - vector of degrees of freedoms for nodes
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
conv_segments = [14 1 17];
edges_conv = [];
for i = 1:size(er,2)
    if ismember(er(3, i), conv_segments)
        edges_conv = [edges_conv er(1:2,i)];
    end
end

% ---------------------------------------------------------

% Calculate head distribution stiffness matrix

K = zeros(ndof);
fb = zeros(ndof, 1);

% Add convection matrix Kc to K

for i = 1:length(edges_conv)
    n1 = edges_conv(1, i); % edge node 1
    n2 = edges_conv(2, i); % edge node 2
    
    p1 = coord(n1, :); % edge point 1
    p2 = coord(n2, :); % edge point 2
    
    L = dist(p1,p2'); % distance between points
    
    Kce = thickness * alpha_c * (L/6) * [2, 1; 1, 2]; % Element conv. matrix
    fbe = T_inf * alpha_c * thickness * [L/2 ; L/2]; % Element boundary force
    
    % K = assem(edof(i,:), K, Kc);
    
    % Kt = assem(edof(el,:),Kt,Kte);
    
    % Insert Kc at correct position
    K([n1,n2],[n1,n2]) = K([n1,n2],[n1,n2]) + Kce;
    
    fb([n1, n2]) = fb([n1, n2]) + fbe;
    
end

Q_vec = [0, 0 ,0 , Q]; % Only in the silicon part, subd. 4

fl = zeros(ndof, 1); % global load vector

% Build K matrix
for elnr = 1:nelm
    % Get correct k constant-matrix
    sd = t(4,elnr); % subdomain
    D = eye(2)*k_const(sd); % cond. matrix
    eq = Q_vec(sd); % heat source
    
    [Ke, fe] = flw2te(Ex(elnr, :), Ey(elnr, :), thickness, D, eq);
    
    % Assemble K
    [K, fl]= assem(edof(elnr,:), K, Ke, fl, fe);
    
end

f = fb + fl;

a = solveq(K, f);


% ---- SOLVE TORSION PROBLEM ---- % -------------

ndof_S = ndof * 2;

% Find which edges that are locked (ux = 0, uy = 0)

%conv_segments = [10 11 12]; % Choosen boundary segments
lock_segment = 18;
edges_lock = [];
for i = 1:size(er,2)
    if ismember(er(3,i),lock_segment)
        edges_lock = [edges_lock er(1:2,i)];
    end
end


% Find which edges on symmetry line (ux = 0)

sym_segments = [21, 22, 23, 24];
edges_sym = [];

for i = 1:size(er, 2)
   if ismember(er(3,i), sym_segments)
       edges_sym = [edges_sym er(1:2, i)];
   end
end

% degrees of freedom at bottom boundary, ux=0, uy=0
dof_lock = unique(dof_S(edges_lock, :));

% nodes which should have u_x = 0;
edges_sym = unique(edges_sym);

% degrees of freedom at left boundary, ux = 0. ACTUALLY SAME
dof_sym = dof_S(edges_sym, 1);

% Boundary conditins for solveq
bc = [dof_sym zeros(length(dof_sym), 1)
      dof_lock zeros(length(dof_lock), 1)];


% Calculate average temperature for each element
% AND calc elongatin due to temperature for each element

delta_T = zeros(nelm, 1);
temp_T0 = ones(nen,1)*T_0;

% eps(elnr, i): elongation for elnr in i-direction
epsilon_dT = zeros(nelm, 3);

for elnr = 1:nelm
   sd = t(4,elnr); % subdomain
   
   indx = edof(elnr, 2:end);
   
   T_temp = a(indx) - temp_T0;

   delta_T(elnr) = sum(T_temp) / length(T_temp);
   
   % Add to epsilon matrix
   c = (1 + ny_const(sd)) * alpha_const(sd) * delta_T(elnr);
   epsilon_dT(elnr, :) = c * [1, 1, 0]; % elongation, epsilon_0
   
end


% Assemble K_S matrix

K_S = zeros(ndof_S);

% input to element stiffness matrix function
% 2 indicates strain
ep = [2, thickness];

% element loads
f_S = zeros(ndof_S, 1);

for elnr = 1:nelm
    
   sd = t(4,elnr); % subdomain
   
   % Calc constitutive matrix
   D_temp = calc_dstrain(e_const(sd), ny_const(sd));
   
   Ke_S = plante(Ex(elnr, :), Ey(elnr, :),  ep, D_temp);
   
   % plantf, (D*epsilon^(Dt))^T, not actually for this purpose but works.
   fe_dT = plantf(Ex(elnr, :), Ey(elnr, :), ep, (D_temp * epsilon_dT(elnr,:)')');
   
   % Assemble K_S
   [K_S, f_S]= assem(edof_S(elnr,:), K_S, Ke_S, f_S, fe_dT);
    
end

u = solveq(K_S, f_S, bc); % node displacements

% ---- PLOT ---- %

ed = extract(edof_S, u);

% Calculate displaced coordinates
mag = 100; % Magnification (due to small deformations)
exd = Ex + mag*ed(:,1:2:end);
eyd = Ey + mag*ed(:,2:2:end);


% figure()
% patch(Ex',Ey',[0 0 0],'EdgeColor','none','FaceAlpha',0.3)
% patch(-Ex',Ey',[0 0 0],'EdgeColor','none','FaceAlpha',0.3)
% hold on
% patch(exd',eyd',[0 0 0],'FaceAlpha',0.3)
% patch(-exd',eyd',[0 0 0],'FaceAlpha',0.3)
% axis equal
% title('Displacement field [Magnitude enhancement 100]')

% ---- FIND STRESSES ---- %

element_stresses = zeros(nelm, 3);

for elnr = 1:nelm
    sd = t(4,elnr); % subdomain
    
    D_temp = calc_dstrain(e_const(sd), ny_const(sd)); % Const. matrix
    
    % Store element stresses, D*epsilon
    element_stresses(elnr, :) = plants(Ex(elnr, :), Ey(elnr, :), ep, D_temp, ed(elnr, :));
    
    % TODO subtract D*eps_0 from el_stresses
    
    element_stresses(elnr, :) = ... 
        element_stresses(elnr,:) - (D_temp*epsilon_dT(elnr,:)')';
    
end

% Calculate von mises stress for each element

Seff_el = zeros(nelm, 1); % array with element von m. stresses
for elnr = 1:nelm
    sd = t(4,elnr); % subdomain
    
    Seff_el(elnr) = calc_mises(element_stresses(elnr, :), ... 
        alpha_const(sd), e_const(sd), delta_T(sd), ny_const(sd)); % store

end

% average out stresses in elements for nodes
Seff_nod = zeros(size(coord, 1), 1);
for i=1:size(coord,1)
    
    [c0, c1] = find(edof(:, 2:4) == i);
    Seff_nod(i, 1) = sum(Seff_el(c0)) / size(c0, 1);
    
end
%% plot
ed_Seff = extract(edof, Seff_nod);
figure
patch(Ex', Ey', ed_Seff')
hold on
patch(-Ex', Ey', ed_Seff')
colorbar
xlabel('x-position [m]')
ylabel('y-position [m]')
title('von Mises stress field [Pa]')

Smin = 5.84e+06;
Smax = 1.48e+08;
caxis([Smin Smax])


