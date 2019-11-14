
clear all 
%alpha_guess=1;

%Xk=1;
 %U0=0;
 %dk=1;
 %F=@laplaceeqn;
  %gamma = 1e-4;
   %  delta = 0.5;
    % rhok  = 1e-8;
 %[alpha] = backtr(alpha_guess,Xk,dk,F,gamma,delta,rhok);
  
 % alpha=.1;
U=1:1:25
out=laplaceeqn1(U);

   % admTransform(@laplaceeqn)
opts=admOptions('i',U);

Jac = admDiffFor(@laplaceeqn1, 1,opts);
%[a_out]=creatFullGradients(out)
%[a_U,r]=a_laplaceeqn(U,a_out)
%addiff(@laplaceeqn1,U)
%U=U0+alpha*(-diff(U))
%out=laplaceeqn(U)
%out=laplaceeqn(alpha_guess)