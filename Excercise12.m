%% Part 1, 3 elements

% Snittarea
A = 10;
% Konduktivitet
k = 5;
% V?rmef?rdelning (last?)
Q = 100;

x_L = 2;
x_R = 8;

% Antal element
nelm = 20;
len = x_R - x_L;

% L?ngd p? element
L = len / nelm;

% Koordinatvektor f?r noder
% Coord = [2; 4; 6; 8];

Coord = x_L:L:x_R;

Coord = Coord';

% Frihetsgradet f?r elementen
% Edof = [1  1 2
%         2  2 3
%         3  3 4];
    
Edof = zeros(nelm, 3);

for i = 1:nelm
   
    Edof(i,:) = [i, i, i + 1];
    
end

% Antal noder per element
nen = 2;
ndof = nelm + 1;

% Frihetsgrader f?r noder
% Dof = [1; 2; 3; 4];

% Koordinater f?r elementen
% Ex = coordxtr(Edof, Coord, Dof, nen);

% Styvhetsmatris
K = zeros(ndof);

% Element-lastvektor (V?rmek?llor)

efl = ones(nen, nelm) * Q;

% Fl?desvektor
qxf = -15; % Varf?r negativ?
q = zeros(ndof,1);
q(end) = qxf * A;

% (Element) boundary-vektor (fl?de)
efb = zeros(nen, nelm);

% Assign correct boundary flow values
for i = 1:length(efb)
efb(:,i) = [q(i); q(i+1)];
end

% Assemble into global force vectors
fb = zeros(ndof,1);
fl = zeros(ndof,1);

% ???? Vad ?r detta, fick fel innan utan
% Kraft per element
fe = Q*L/2;

% Solver %%%%%%%%%%%%%%%%%%%

% Assemble global stiffness matrix

Ke = spring1e((k*A)/L);

for elnr = 1:nelm
    
    [K, fl] = assem(Edof(elnr,:), K, Ke, fl, fe);
    
end

fb = insert(Edof,fb,efb);

f = fb + fl;

% Randvillkor? Endast nod 1 som har randv. T = 0.



bc = [1 0];

a = solveq(K,f,bc)

% plot exact vs fem

T_x = @(x) -x.^2+13.*x-22;
x = linspace(x_L,x_R,ndof);
x_fine = linspace(x_L,x_R,1000);

hold on
plot(x,a,'or');
plot(x_fine,T_x(x_fine));
