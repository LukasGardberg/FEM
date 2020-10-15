
%% PART 2

% geom2

k = 1;
ep = 1; % Thickness
D = eye(2)*k; % Homogenous material

F = zeros(ndof,1); % No heat sources?
nen = 3; % number of element nodes

% Symmetric boundary conditions?
% zero heat flux on shared boundaries?

efl = zeros(nen, nelm); % Element load forces
efb = zeros(nen, nelm); % Element boundary forces

fb = zeros(ndof,1);
fl = zeros(ndof,1);

fb = insert(edof,fb,efb);
fl = insert(edof,fl,efl);

f = fb + fl;

% Solver

K = zeros(ndof);

for i = 1:nelm
    Ke = flw2te(ex(i,:), ey(i,:), ep, D);
    K = assem(edof(i,:), K, Ke);
end

a = solveq(K, f, bc);

% Plot

figure(2)
ed = extract(edof, a);
patch(ex', ey', ed');
colormap jet
colorbar
