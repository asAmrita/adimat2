
clear all 
%alpha_guess=1;
%Xk=1;
 U0=0;
 %dk=1;
 %F=@laplaceeqn;
  %gamma = 1e-4;
   %  delta = 0.5;
    % rhok  = 1e-8;
 %[alpha] = backtr(alpha_guess,Xk,dk,F,gamma,delta,rhok);
  %out=laplaceeqn(U);
  alpha=.1;
U=ones(60,60);
Jac = admDiffFor(@laplaceeqn, 1, U)


%U=U0+alpha*(-diff(U))
%out=laplaceeqn(U)
%out=laplaceeqn(alpha_guess)