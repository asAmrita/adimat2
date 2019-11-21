 clear all 
 clc
 UMin= -1;
 UMax= 1;
 U=zeros(8,8);
 nx = size(U,2);
 ny = size(U,1);
% 
% dx=1/(nx-1);                     %Width of space step(x)
% dy=1/(ny-1);                     %Width of space step(y)
% x=0:dx:1;                        %Range of x(0,2) and specifying the grid points
% y=0:dy:1;

 %dk=admDiffVFor(@LaplaceEqn, 1,U)
 %B = reshape(dk, 8, 8)
 alpha_guess=100;
 F= @LaplaceEqn;

 

while 1

   B=(admDiffVFor(@LaplaceEqn, 1,U));
   dk = reshape(B, 8, 8);
   [alpha] = backtr(alpha_guess,U,dk,F);
   
   if (norm(alpha*dk)/norm(U) < 1e-5) 
       break
   end
   
   % U = U + alpha*(dk)
  
   U = U - alpha*(-dk);
   U=min(U,UMax);
   U=max(U,UMin);
   
   
end
% surf(x,y,U)

%   end 