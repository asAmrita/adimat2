
function out=laplaceeqn1(U)

% A = 10*rand(size(U));
% b = 10*rand(size(U,1));
% 
% y = A \ (U*b);
% 
% out = sum(sum(y.^2));
%p=1;
%out=p+U;

%out=U;



%U=ones(25,25)

 nx=5;                           %Number of steps in space(x)
 ny=5;                           %Number of steps in space(y)       
 niter=2;                     %Number of iterations 
 dx=2/(nx-1);                     %Width of space step(x)
 dy=2/(ny-1);                     %Width of space step(y)
 x=0:dx:2;                        %Range of x(0,2) and specifying the grid points
 y=0:dy:2;                        %Range of y(0,2) and specifying the grid points
% %%
% %Initial Conditions
  %P=zeros(ny,nx);                  %Preallocating p
% pn=zeros(ny,nx);                 %Preallocating pn
% %%
% %Boundary conditions
 P(:,1)=0;
 P(:,nx)=0;
 P(1,:)=0;
 P(2,:)=0;                   %Neumann conditions
 P(ny,:)=0;
 P(ny-1,:)=0;               ...same as above
%%
% %Explicit iterative scheme with C.D in space (5-point difference)
 j=2:nx-1;
 i=2:ny-1;
 for it=1:niter
      
     pn=P;
    P(i,j)=((dy^2*(pn(i+1,j)+pn(i-1,j)))+(dx^2*(pn(i,j+1)+pn(i,j-1))))/(2*(dx^2+dy^2))-(P(i,j));
% %     %Boundary conditions (Neumann conditions)
     P(:,1)=0;
     P(:,nx)=0;
     P(1,:)=0;
     P(2,:)=0;
    P(ny,:)=0;
    P(ny-1,:)=0;
    
%    %U(:,1)=0;
%     %U(:,nx)=0;
%     %U(1,:)=U(2,:);
%     %U(ny,:)=U(ny-1,:);
%    % p=U;
%    
%    
%

 
 end
 out=U+P;
end
% 
% %%
% %Plotting the solution
% plot(x,out);       
% %shading interp
% title({'2-D Laplace''s equation';['{\itNumber of iterations} = ',num2str(it)]})
% xlabel('Spatial co-ordinate (x) \rightarrow')
% %ylabel('{\leftarrow} Spatial co-ordinate (y)')
% %zlabel('Solution profile (P) \rightarrow') 
