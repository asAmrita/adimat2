
clear all 
clc
UMin= -1;
UMax= 1;
%nx = 8;
%ny = 8;

U=zeros(8,8);
% dk=admDiffVFor(@laplaceeqn, 1,U);

alpha_guess=100;
F=@laplaceeqn;
nx = size(U,2);
ny = size(U,1);
dx=1/(nx-1);                     
dy=1/(ny-1);                    
x=0:dx:1;                        
y=0:dy:1;  
 

while 1

   B=(admDiffVFor(@laplaceeqn, 1,U));
   dk =reshape(B, 8, 8);
    [alpha] = backtr(alpha_guess,U,dk,F);
   
   
   if ((norm(alpha*dk)/norm(U) <1e-5))
       break
   end
   
   %U = U - alpha*(dk)
  
   U = U - alpha*(-dk);
   U=min(U,UMax);
   U=max(U,UMin)
   
   
end

%   end 
 surf(x,y,U)
laplaceeqn(U)

