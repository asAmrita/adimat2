function [alpha] = backtr(alpha_guess,Xk,dk,F)
%function [alpha] = backtr(alpha_guess,Xk,dk,F,gamma,delta,rhok)
% INPUT:
%       NOTE: (*) indicates necessary input, the other variables are optional 
%       (*) alpha _guess - current steplength (1*1) [>0];
%       (*) Xk           - current iterate    (N*1);
%       (*) dk           - search direction   (N*1);
%           gamma        - constant provided by the user (1*1) [>0];
%           delta        - constant provided by the user (1*1) into the range [0,  1];
%           rhok         - constant provided by the user (1*1) into the range [0,  1];
%       (*) F            - function handle of the objective function (RN->R );


if (nargin < 5)
    gamma = 1e-4;
    delta = 0.5;
    rhok  = 1e-8;
elseif (nargin < 6)
    delta = 0.5;
    rhok = 1e-8;
elseif (nargin < 7)
    rhok = 1e-8;
end



        % positive direction (+)alpha
    alpha = alpha_guess;
    while (F(Xk+alpha.*dk)>F(Xk)-gamma*alpha^2*(norm(dk))^2) 
        if (alpha*norm(dk) < rhok)   
            alpha  = 0;              % <-- failure to search for a value of alpha nonzero
        else
            alpha = alpha*delta;     % <-- reduction of the steplength
        end
    end 
    alpha1 = alpha;
    F1     = F(Xk+alpha1.*dk)-(F(Xk)-gamma*alpha1^2*(norm(dk))^2);
        % negative direction (-)alpha
    alpha = alpha_guess;
    while (F(Xk-alpha.*dk)>F(Xk)-gamma*alpha^2*(norm(dk))^2)  
        if (alpha*norm(dk) < rhok)
            alpha   = 0;              % <-- failure to search for a value of alpha nonzero
        else
            alpha = alpha*delta;      % <-- reduction of the steplength
        end
    end
    alpha2 = -alpha;
    F2     = F(Xk+alpha2.*dk)-(F(Xk)-gamma*alpha2^2*(norm(dk))^2);
    % choice of the value of alpha for which it is provided with sufficient reduction 
    if (F1<F2)           
        alpha = alpha1;
    else
        alpha = alpha2;
    end  
    %out=laplaceeqn(alpha);
    %U=U0+alpha*(out);
   
    
end
