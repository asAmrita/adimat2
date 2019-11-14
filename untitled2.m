function out=laplaceEqn(U)

% A = 10*rand(size(U));
% b = 10*rand(size(U,1));
% 
% y = A \ (U*b);
% 
% out = sum(sum(y.^2));



pw=25;

%out=0;
U0=0;
lembda=1;
nx=60;                           %Number of steps in space(x)
ny=60;                           %Number of steps in space(y)       
niter=1;                     %Number of iterations 
dx=2/(nx-1);                     %Width of space step(x)
dy=2/(ny-1);                     %Width of space step(y)
x=0:dx:2;                        %Range of x(0,2) and specifying the grid points
y=0:dy:2;                        %Range of y(0,2) and specifying the grid points
%%
%Initial Conditions
p=zeros(ny,nx);                  %Preallocating p
pn=zeros(ny,nx);                 %Preallocating pn
%%
%Boundary conditions
p(:,1)=0;
p(:,nx)=0;
p(1,:)=0;                   %Neumann conditions
p(ny,:)=0;               ...same as above
%%
%Explicit iterative scheme with C.D in space (5-point difference)
j=2:nx-1;
i=2:ny-1;
for it=1:niter
     
    pn=p;
    p(i,j)=((dy^2*(pn(i+1,j)+pn(i-1,j)))+(dx^2*(pn(i,j+1)+pn(i,j-1))))/(2*(dx^2+dy^2))-(U(i,j))
    %Boundary conditions (Neumann conditions)
    p(:,1)=0;
    p(:,nx)=0;
    p(1,:)=0;
    p(ny,:)=0;
   % p=U;
   
end

out=(sum(abs(p-pw).^2*dx*dy)+lembda/2*sum(abs(U).^2)*dx*dy);
 
% 
% %%
% %Plotting the solution
% plot(x,out);       
% %shading interp
% title({'2-D Laplace''s equation';['{\itNumber of iterations} = ',num2str(it)]})
% xlabel('Spatial co-ordinate (x) \rightarrow')
% %ylabel('{\leftarrow} Spatial co-ordinate (y)')
% %zlabel('Solution profile (P) \rightarrow')
end