%% FEM labb uppg1
A = 10;
k = 5;
Le = 2; %length of the element
Q = 100; % givet
K = zeros(4);
fl = zeros(4, 1); %f load vector
edof = [1 1 2; 2 2 3; 3 3 4];

X = [2; 4; 6; 8];
ep = A*k/Le; % spring stiffness - men vi har ingen s? ersatts med detta varde
Ke = spring1e(ep); % K for every node per element
fe = Q*Le/2; % load for every element? 

 for i = 1:3
     %K = assem(edof(i,:), K, Ke);
     [K,fl] = assem(edof(i, :),K,Ke,fl,fe); % calculates stiffness matrix, K, and the load vector, fl (integralen)
 end

q4 = 15; %givet
fb = [0; 0; 0; -A*q4]; % vet ej fb(1), s?tter den bara till 0
f = fl + fb;
 
bc = [1 0]; % berattar att vi inte vet fb(1) ?  
[a, qA] = solveq(K, f, bc); % calculates T-vector, and q*A

q = qA/A + fb/A; %hela q-vektorn

%KT = K*a; %vet inte hela f, EMMAS
%q0 = (KT(1) - fl(1))/A;

%% 20 element
nelm = 20; % nbr of elements 
A = 10;
k = 5;
L = 6; %total length
Le = L/nelm; % length of each element
Q = 100; %givet
K = zeros(nelm + 1); % nelm + 1 = nbr of nodes
fl = zeros(nelm + 1, 1);
edof = zeros(nelm, 3); 

 for i = 1:nelm % which elements are connected to which nodes
     edof(i, :) = [i i (i+1)];
 end

ep = A*k/Le; 
Ke = spring1e(ep); 
fe = Q*Le/2; % OKLART

 for i = 1:nelm % For every loop something adds to both K and fl
     %K = assem(edof(i,:), K, Ke);
     [K,fl] = assem(edof(i, :),K,Ke,fl,fe);
 end

qend = 15; % boundary value
fb = zeros(nelm + 1, 1); % boundary f-vector
fb(end) = -A*qend; % adds the boundary value
f = fl + fb; % f-vector
 
bc = [1 0]; % Vill veta flodet i forsta noden. Det har visar det  
[a qA]= solveq(K, f, bc);

q = qA(1)/A;


%KT = K*a;
%q0 = (KT(1) - fl(1))/A;

T_x = @(x) -x.^2+13.*x-22;
x = linspace(2,8,nelm + 1);
x_fine = linspace(2,8,1000);

hold on
plot(x,a,'or');
plot(x_fine,T_x(x_fine));