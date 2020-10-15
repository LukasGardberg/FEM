%% FEM labb uppg2

clf

% defined coordinates

nf=7;
nr=6;  % must be even
R1=2;
R2=10;
coord=[];
for R=R1:(R2-R1)/nf:R2
for fi=0:pi/2/nr:pi/2
   if fi<pi/4
     Rf=(R1-R)/(R1-R2)*R/cos(fi)+(R2-R)/(R2-R1)*R;
     coord=[coord;Rf*cos(fi) Rf*sin(fi)];
   else
     Rf=(R1-R)/(R1-R2)*R/sin(fi)+(R2-R)/(R2-R1)*R;
     coord=[coord;Rf*cos(fi) Rf*sin(fi)];
   end
end
   
end
plot(coord(:,1),coord(:,2),'*')
axis equal
hold on

% global dof for nodal points

ndof=0;
dof=[];
for R=R1:(R2-R1)/nf:R2
  for fi=0:pi/2/nr:pi/2
    ndof=ndof+1;
    dof=[dof; ndof];
  end
end

% topology matriz

nelm=0;

edof=[];
for fn=0:nf-1
 
  for ir=0:nr-1
     n1=fn*(nr+1)+1+ir;
     n2=fn*(nr+1)+2+ir;
     n3=(fn+1)*(nr+1)+1+ir;
     n4=(fn+1)*(nr+1)+2+ir;    

     nelm=nelm+1;
     edof=[edof; nelm n3 n2 n1];
     nelm=nelm+1;
     edof=[edof; nelm n3 n4 n2];
  end
end

[ex,ey]=coordxtr(edof,coord,dof,3);
eldraw2(ex,ey,[1 2 0],edof(:,1))
%pause

% Define boundary conditions
%
% find nodal points on boundaries

bc3=find(abs(coord(:,1)-R2)<1e-3);
bc4=find(abs(coord(:,2)-R2)<1e-3);
bcarc=find(abs(coord(:,1).^2+coord(:,2).^2-R1^2)<1e-3);


bcarc(:,2)=1000; % circle
bc3(:,2)=100;    % top boundary
bc4(:,2)=100;    % right boundary

bc=[bcarc;bc3;bc4];


bc(11,:)=[]; % Remove 


%% Min kod ...
k = 1; %konduktivitet
t = 1; %tjocklek, 
Tin = 1000;
Tout = 100;
R1 = 2;
R2 = 10;

nodes = length(coord); % also ndof = number of degrees of freedom
D = k*eye(2); % hur styvhetsmatrisen ser ut
K = zeros(nodes); 
fl = zeros(nodes, 1); 

ep = t; % tjocklek
eq = 0;

 for i = 1:nelm     
     [Ke, fe] = flw2te(ex(i, :), ey(i, :), ep, D, eq);
     %[Ke] = flw2te(ex(i, :), ey(i, :), ep, D);
     [K, fl]= assem(edof(i,:), K, Ke, fl, fe);
 end

AreaSquare = R2*R2;
AreaCircle = R1^2*pi;
fb = zeros(nodes, 1);
fb(1) = AreaCircle*Tin;
fb(end) = -AreaSquare*Tout;

f = fb + fl;
 
[a, q] = solveq(K, f, bc);

ed = extract(edof, a);

figure(1)
patch(ex', ey', ed')
colormap jet
colorbar