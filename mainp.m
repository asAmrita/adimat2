clear all
clc
x=-1.5;
y=20.5;
U=[x,y];

alpha_guess=1;

F=@para;

while 1


   dk=-(admDiffVFor(@para, 1,U));

   [alpha] = backtr(alpha_guess,U,dk,F)
   
   if (norm(alpha*dk)/norm(U) < 1e-10)
       break
   end
   
   U = U + alpha*(dk);

 
 
 
end
