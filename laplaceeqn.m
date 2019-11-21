function out=laplaceeqn(U)

nx = size(U,2);
ny = size(U,1);

%nx=25;                           %Number of steps in space(x)
%ny=25;                           %Number of steps in space(y)       
niter=1000;                     %Number of iterations 
dx=1/(nx-1);                     %Width of space step(x)
dy=1/(ny-1);                     %Width of space step(y)
x=0:dx:1;                        %Range of x(0,2) and specifying the grid points
y=0:dy:1;                        %Range of y(0,2) and specifying the grid points
%p=@(x,y)(min(max(1/4.*pi^2.*sin(2.*pi.*x).*sin(2.*pi.*y),ymin),ymax));
%U=@(x,y)(-sign(-1/128.*pi^2*sin(8.*pi.*x)*sin(8.*pi.*y)));
M=@(x,y)(sin(pi.*x).*sin(pi.*y)+sin(8.*pi.*x).*sin(8.*pi.*y));
f=@(x,y)(2.*pi^2.*sin(pi.*x).*sin(pi.*y)+sign(-sin(8.*pi.*x).*sin(8.*pi.*y)));
%Initial Conditions
 P=zeros(ny,nx);                  %Preallocating p
 pn=zeros(ny,nx);
 %M=zeros(ny,nx);
% pw=zeros(ny,nx);
 %pw=zeros(ny,nx);
 %Preallocating pn
% % %
% % %Boundary conditions
  P(:,1)=0;
  P(:,nx)=0;
  P(1,:)=0;                       %Neumann conditions
  P(ny,:)=0; 
%  
%  
%                 ...same as abovey
% % % %
% % % Explicit iterative scheme with C.D in space (5-point difference)
  j=2:nx-1;
  i=2:ny-1;
  for it=1:niter
        
    pn=P;
   
    % P(i,j)=((dy^2*(pn(i+1,j)+pn(i-1,j)))+(dx^2*(pn(i,j+1)+pn(i,j-1))))/(2*(dx^2+dy^2))-(U(i,j)*dx^2*dy^2)/(2*(dx^2+dy^2))-(f(x(i)',y(j))*dx^2*dy^2)/(2*(dx^2+dy^2));
    P(i,j)=((dy^2*(pn(i+1,j)+pn(i-1,j)))+(dx^2*(pn(i,j+1)+pn(i,j-1)))-(U(i,j)+f(x(i)',y(j))*dx^2*dy^2))/(2*(dx^2+dy^2));
  
%  
     
% %  %Boundary conditions (Neumann conditions)
      P(:,1)=0;
      P(:,nx)=0;
      P(1,:)=0;
      P(ny,:)=0;
% 
%    
  end
%   j=1:nx;
%  i=1:ny;
%  for it =1:niter
%      pw(i,j)=M(i,j);
%  end
pw = M(x',y);

% % % 
  out= sum(sum(0.5*(P - pw).^2*dx*dy));
% % %  Plotting the solution
% %  h=surf(x,y,U','EdgeColor','none');       
% %  shading interp
% % % axis([-0.5 2.5 -0.5 2.5 -100 100])
% %  title({'2-D Poisson equation';['{\itNumber of iterations} = ',num2str(it)]})
% %  xlabel('Spatial co-ordinate (x) \rightarrow')
% %  ylabel('{\leftarrow} Spatial co-ordinate (y)')
% %  zlabel('Solution profile (P) \rightarrow')
% 
% %surf(x,y,U)
 end