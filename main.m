
clear all 
clc
UMin= -1;
UMax= 1;
U=zeros(5,5);
% dk=admDiffVFor(@laplaceeqn, 1,U);

 alpha_guess=100;
 F=@laplaceeqn;

 

while 1

   B=-(admDiffVFor(@laplaceeqn, 1,U));
   dk = reshape(B, 5, 5);
   [alpha] = backtr(alpha_guess,U,dk,F);
   
   if (norm(alpha*dk)/norm(U) < 1e-3) 
       break
   end
   
   %U = U + alpha*(dk)
  
   U = U + alpha*(dk);
   U=min(U,UMax);
   U=max(U,UMin);
   
   
end

%   end 
surf(x,y,U)
laplaceeqn(U)

