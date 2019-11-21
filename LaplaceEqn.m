%cifying parameters
function out=LaplaceEqn(U)
nx=8;                          %Number of steps in space(x)
ny=8;                          %Number of steps in space(y)       

%ny = size(U,1);

dx=1/(nx-1);                     %Width of space step(x)
dy=1/(ny-1);                     %Width of space step(y)
x=0:dx:1;                        %Range of x(0,2) and specifying the grid points
y=0:dy:1;                        %Range of y(0,2) and specifying the grid points
UW=0;                            %x=0 Dirichlet B.C 
UE=0;                            %x=L Dirichlet B.C 
US=0;                            %y=0 Dirichlet B.C 
UN=0;                            %y=L Dirichlet B.C 
%UnW=0;                           %x=0 Neumann B.C (du/dn=UnW)
%UnE=0;                           %x=L Neumann B.C (du/dn=UnE)
%UnS=0;                           %y=0 Neumann B.C (du/dn=UnS)
%UnN=0;                           %y=L Neumann B.C (du/dn=UnN)
u=zeros(nx,ny); 
%U=zeros(nx,ny);
%Pre-allocating u
%%
%B.C vector
%bc=zeros(nx-2,ny-2);
bc(1,:)=UW ;
bc(nx,:)=UE ; %Dirichlet B.Cs
bc(:,1)=US;
bc(:,ny)=UN ; %Dirichlet B.Cs
%bc(1,:)=-UnW/dx; bc(nx-2,:)=UnE/dx;  %Neumann B.C
%bc(:,1)=-UnS/dy; bc(:,ny-2)=UnN/dy;  %Neumann B.Cs
%B.Cs at the corners:
%bc(1,1)=UW/dx^2; bc(nx,1)=UE/dx^2;
%bc(1,ny)=UW/dx^2; bc(nx,ny)=UE/dx^2;
%Calculating the coefficient matrix for the implicit scheme
Ex=sparse(2:nx-2,1:nx-3,1,nx-2,nx-2);
Ax=Ex+Ex'-2*speye(nx-2);        %Dirichlet B.Cs
%Ax(1,1)=-1; Ax(nx-2,nx-2)=-1;  %Neumann B.Cs
Ey=sparse(2:ny-2,1:ny-3,1,ny-2,ny-2);
Ay=Ey+Ey'-2*speye(ny-2);        %Dirichlet B.Cs
%Ay(1,1)=-1; Ay(ny-2,ny-2)=-1;  %Neumann B.Cs
A=kron(Ay/dy^2,speye(nx-2))+kron(speye(ny-2),Ax/dx^2);
%%

%Source term


f=@(x,y)(2.*pi^2.*sin(pi.*(x)).*sin(pi.*(y))+sign(-sin(8.*pi.*(x)).*sin(8.*pi.*(y))));
S = f(x',y);
S(:,nx)=0;
S(nx,:)=0;
S= S(any(S,2),any(S));

%S=reshape(S,nx-2,ny-2)
%S = U;
%S=zeros(nx-2,ny-2);

S=reshape(S-bc,[],1);
S=A\S
S=reshape(S,nx-2,ny-2);
u(2:nx-1,2:ny-1)=S;
%Boundary conditions
%Dirichlet:
u(1,:)=UW;
u(:,1)=US;
u(nx,:)=UE;
u(:,ny)=UN;
                        %Range of x(0,2) and specifying the grid points

M=@(x,y)(sin(pi.*x).*sin(pi.*y)+sin(8.*pi.*x).*sin(8.*pi.*y));
pw = M(x',y);
% 
out= sum(sum(0.5*(u - pw).^2*dx*dy));

%u(1,:)=u(2,:)-UnW*dx;
%u(nx,:)=u(nx-1,:)+UnE*dx;
%u(:,1)=u(:,2)-UnS*dy;
%u(:,ny)=u(:,ny-1)+UnN*dy;
%%
%Plotting the solution
% 
% surf(x,y,u','EdgeColor','none');
% shading interp
% title('2-D Laplace''s equation')
% xlabel('Spatial co-ordinate (x) \rightarrow')
% ylabel('{\leftarrow} Spatial co-ordinate (y)')
% zlabel('Solution profile (P) \rightarrow')
end
