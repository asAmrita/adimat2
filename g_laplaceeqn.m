% Generated by ADiMat 0.6.6-5530 (00419e1f)
% Copyright 2009-2013 Johannes Willkomm, Fachgebiet Scientific Computing,
% TU Darmstadt, 64289 Darmstadt, Germany
% Copyright 2001-2008 Andre Vehreschild, Institute for Scientific Computing,
% RWTH Aachen University, 52056 Aachen, Germany.
% Visit us on the web at http://www.adimat.de
% Report bugs to johannes@johannes-willkomm.de
%
%
%                             DISCLAIMER
%
% ADiMat was prepared as part of an employment at the Institute
% for Scientific Computing, RWTH Aachen University, Germany and is
% provided AS IS. NEITHER THE AUTHOR(S), THE GOVERNMENT OF THE FEDERAL
% REPUBLIC OF GERMANY NOR ANY AGENCY THEREOF, NOR THE RWTH AACHEN UNIVERSITY,
% INCLUDING ANY OF THEIR EMPLOYEES OR OFFICERS, MAKES ANY WARRANTY,
% EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY
% FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION OR
% PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
% PRIVATELY OWNED RIGHTS.
%
% Global flags were:
% FORWARDMODE -- Apply the forward mode to the files.
% NOOPEROPTIM -- Do not use optimized operators. I.e.:
%		 g_a*b*g_c -/-> mtimes3(g_a, b, g_c)
% NOLOCALCSE  -- Do not use local common subexpression elimination when
%		 canonicalizing the code.
% NOGLOBALCSE -- Prevents the application of global common subexpression
%		 elimination after canonicalizing the code.
% NOPRESCALARFOLDING -- Switch off folding of scalar constants before
%		 augmentation.
% NOPOSTSCALARFOLDING -- Switch off folding of scalar constants after
%		 augmentation.
% NOCONSTFOLDMULT0 -- Switch off folding of product with one factor
%		 being zero: b*0=0.
% FUNCMODE    -- Inputfile is a function (This flag can not be set explicitly).
% NOTMPCLEAR  -- Suppress generation of clear g_* instructions.
% UNBOUND_XML  -- Write list of unbound identifiers in XML format.
% DEPENDENCIES_XML  -- Write list of functions in XML format.
% UNBOUND_ERROR	-- Stop with error if unbound identifiers found (default).
% FUNCTION_LIST_XML	-- Write list of functions to XML file.
% VERBOSITYLEVEL=5
% AD_IVARS= U
% AD_DVARS= out

function [g_out, out]= g_laplaceeqn(g_U, U)
   
   % A = 10*rand(size(U));
   % b = 10*rand(size(U,1));
   % 
   % y = A \ (U*b);
   % 
   % out = sum(sum(y.^2));
   p= 1; 
   g_out= g_U+ g_zeros(size(p));
   out= p+ U; 
   %ymin= -1;
   %ymax= 1;
   
   % pw=25;
   % U=zeros(5,5);
   % nx = size(U,2);
   % ny = size(U,1);
   % 
   % 
   % lambda=1;
   % %nx=25;                           %Number of steps in space(x)
   % %ny=25;                           %Number of steps in space(y)       
   % niter=1;                     %Number of iterations 
   % dx=2/(nx-1);                     %Width of space step(x)
   % dy=2/(ny-1);                     %Width of space step(y)
   % x=0:dx:2;                        %Range of x(0,2) and specifying the grid points
   % y=0:dy:2;                        %Range of y(0,2) and specifying the grid points
   % %p=@(x,y)(min(max(1/4.*pi^2.*sin(2.*pi.*x).*sin(2.*pi.*y),ymin),ymax));
   % %U=@(x,y)(-sign(-1/128.*pi^2*sin(8.*pi.*x)*sin(8.*pi.*y)));
   % %f=@(x,y)(2.*pi^2.*sin(pi.*x).*sin(pi.*y)+sign(-sin(8.*pi.*x).*sin(8.*pi.*y)))
   % %Initial Conditions
   %  P=zeros(ny,nx);                  %Preallocating p
   %  pn=zeros(ny,nx);                 %Preallocating pn
   % % %
   % % %Boundary conditions
   %  P(:,1)=0;
   %  P(:,nx)=0;
   %  P(1,:)=0;    %P(2,:);                   %Neumann conditions
   %  P(ny,:)=0;   %P(ny-1,:);               ...same as abovey
   % % % %
   % % % Explicit iterative scheme with C.D in space (5-point difference)
   %  j=2:nx-1;
   %  i=2:ny-1;
   %  for it=1:niter
   % %      
   %     pn=P;
   %    P(i,j)=((dy^2*(pn(i+1,j)+pn(i-1,j)))+(dx^2*(pn(i,j+1)+pn(i,j-1))))/(2*(dx^2+dy^2))-(P(i,j));
   % %     
   % %     %Boundary conditions (Neumann conditions)
   %      P(:,1)=0;
   %      P(:,nx)=0;
   %      P(1,:)=0;
   %      P(ny,:)=0;
   % 
   %    
   %  end
   % % 
   % % out= sum(sum(0.5*(P - pw).^2*dx*dy+lambda/2*sum(U).^2*dx*dy));
   % % %  Plotting the solution
   % % %  h=surf(x,y,U','EdgeColor','none');       
   % % %  shading interp
   % % % % axis([-0.5 2.5 -0.5 2.5 -100 100])
   % % %  title({'2-D Poisson equation';['{\itNumber of iterations} = ',num2str(it)]})
   % % %  xlabel('Spatial co-ordinate (x) \rightarrow')
   % % %  ylabel('{\leftarrow} Spatial co-ordinate (y)')
   % % %  zlabel('Solution profile (P) \rightarrow')
   % % 
   % % %surf(x,y,U)
end
