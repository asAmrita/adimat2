U=zeros(5,5);
  alpha_guess=1;
 F=@laplaceEqn;

 while 1

   B=-(admDiffVFor(@laplaceEqn, 1,U));
   dk = reshape(B, 5, 5);
   [alpha] = backtr(alpha_guess,U,dk,F);
   
   if (norm(alpha*dk)/norm(U) < 1e-5) 
       break
   end
   
   U = U + alpha*(dk)

end