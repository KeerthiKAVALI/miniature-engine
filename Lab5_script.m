
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close
clear all
clc
%discritization
Nx=21;
Ny=Nx;
L=1;
H=1;
x=linspace(0,L,Nx)';
y=linspace(0,H,Ny)';

%Intitialization of Tolerance
Nmodes = 20;
Niter = 10; 

%Scale factor
s=5*1000;

% compute the K - stiffness
op1 = 1;
op2 = 1;
Kx = FEM_mat_1D(x,op1,op2);
Ky = FEM_mat_1D(y,op1,op2);

% compute M - mass
op1 = 0;
op2 = 0;
Mx = FEM_mat_1D(x,op1,op2);
My = FEM_mat_1D(y,op1,op2);
%Compute the Mixed Matrix Cx (C1x represents Cx "first row second column" and C2x represents Cx "second row first column" )
op1=1;
op2=0;
C1x=FEM_mat_1D(x,op1,op2);
op1=0;
op2=1;
C2x=FEM_mat_1D(x,op1,op2);
%compute the Mixed Matrix Cy (C1y represents Cy "first row second column" and C2y represents Cy "second row first column" )
op1 = 0;
op2 = 1;
C1y = FEM_mat_1D(y,op1,op2);
op1 = 1;
op2 = 0;
C2y = FEM_mat_1D(y,op1,op2);
%defining Nx and Ny and AA 
Nx = numel(x);
Ny = numel(y);
%D is the number of Dimensions (Rows Number)
D=2;
%N is number of Terms (Columns Number)
N=8;
AA=cell(D,N);
%Determining Zeros Matrices (size of Kx or Ky)
z=zeros(size(Ky));
%Material Properties
E = 10^4;
v = 0.3;
k = E/(1-v^2);
p =(1-v)/2;

% Filling the Matrix AA
AA{1,1} = k*[Kx,z;z,z];
AA{2,1} = [My,z;z,z];

AA{1,2} = k*v*[z C1x;z z];
AA{2,2} = [z,C1y;z,z];

AA{1,3} =  k*v*[z,z;C2x,z];
AA{2,3} = [z,z;C2y,z];

AA{1,4} = k*[z,z;z,Mx];
AA{2,4} = [z,z;z,Ky];

AA{1,5} = k*p*[Mx,z;z,z];
AA{2,5} = [Ky,z;z,z];

AA{1,6} = k*p*[z,C2x;z,z];
AA{2,6} = [z,C2y;z,z];

AA{1,7} = k*p*[z,z;C1x,z];
AA{2,7} = [z,z;C1y,z];

AA{1,8} = k*p*[z,z;z,Kx];
AA{2,8} = [z,z;z,My];  


%Filling the Matrix BB
%Case 1: bx=by=1;
% f11 = FEM_mat_1D(x,0,0,[])*ones(Nx,1);
% f21 = FEM_mat_1D(y,0,0,[])*ones(Ny,1);
% f12 = FEM_mat_1D(x,0,0,[])*ones(Nx,1);
% f22 = FEM_mat_1D(y,0,0,[])*ones(Ny,1);

%Case 2: bx=x.^2*y, by=(y-1).^2;
f11 = FEM_mat_1D(x,0,0,[])*(x.*x);
f21 = FEM_mat_1D(y,0,0,[])*y;
f12 = FEM_mat_1D(x,0,0,[])*ones(Nx,1);
f22 = FEM_mat_1D(y,0,0,[])*((y-ones(Ny,1)).*(y-ones(Ny,1)));

BB{1,1} = [f11,zeros(Nx,1);zeros(Nx,1),f12];
BB{2,1} = [f21,zeros(Ny,1);zeros(Ny,1),f22];

%Filling the Matrix GG
F1 = zeros(Nx*2,0);
F2 = zeros(Ny*2,0);
GG{1} = F1;
GG{2} = F2;

%Bords (Boundary Condition)
Bords{1}=1;
Bords{2}=Ny+1;

%Calling the function
FF = easy_PGD(AA,BB,GG,Nmodes,Niter,Bords);

%Displacement along X
ux = FF{1}(1:Nx,:);
vx = FF{1}(Nx+1:Nx*2,:);

%Displacement along Y
uy = FF{2}(1:Ny,:);
vy = FF{2}(Ny+1:Ny*2,:);

%Horizontal and vertical components
u = ux*uy';
v = vx*vy';

%Plot for horizontal displacement
[X,Y] = meshgrid(x,y);
figure(1)
quiver(X,u);
xlabel('Y','fontsize',14)
ylabel('X','fontsize',14)
title('X-displacement displacement plot','fontsize',14)
set(gca,'fontsize',14)
X=X+u'*s;
figure(2)
surface(x,y,u);
colormap;
colorbar;
xlabel('Y','fontsize',14)
ylabel('X','fontsize',14)
title('X-displacement plot','fontsize',14)
set(gca,'fontsize',14)
  
%Plot for vertical displacement
[X,Y] = meshgrid(x,y);
figure(3)
quiver(Y,v);
xlabel('Y','fontsize',14)
ylabel('X','fontsize',14)
title('Y-displacement plot','fontsize',14)
set(gca,'fontsize',14)
Y=Y+v'*s;
figure(4)
surface(X,Y,v);
colormap;
colorbar;
xlabel('Y','fontsize',14)
ylabel('X','fontsize',14)
title('Y-displacement plot','fontsize',14)
set(gca,'fontsize',14)
 
%Magnitude Displacement Plot
Magnitude_disp=sqrt(u.*u+v.*v);
[X,Y] = meshgrid(x,y);
figure(5)
surf(x,y,Magnitude_disp);
xlabel('X','fontsize',14)
ylabel('Y','fontsize',14)
title('Displacement field plot','fontsize',14)
set(gca,'fontsize',14)
figure(6)
surface(X,Y,Magnitude_disp);
xlabel('X','fontsize',14)
ylabel('Y','fontsize',14)
title('Displacement field','fontsize',14)
set(gca,'fontsize',14)
colormap;
colorbar;

%Deformed shape Plot
[x,y] = meshgrid(x,y);
x=x+u'*s;
y=y+v'*s;
figure(7)
quiver(x,y,u,v)
xlabel('X','fontsize',14)
ylabel('Y','fontsize',14)
title('Vector field Deformation','fontsize',14)
set(gca)
figure(8)
surface(x,y,Magnitude_disp);
xlabel('X','fontsize',14)
ylabel('Y','fontsize',14)
title('Shape Deformation','fontsize',14)
set(gca)
colormap;
colorbar;











