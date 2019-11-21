% Generated by ADiMat 0.6.6-5530 (00419e1f)
% © 2001-2008 Andre Vehreschild <vehreschild@sc.rwth-aachen.de>
% © 2009-2018 Johannes Willkomm <johannes@johannes-willkomm.de>
% TU Darmstadt, 64289 Darmstadt, Germany
% Visit us on the web at http://www.adimat.de/
% Report bugs to johannes@johannes-willkomm.de
%
%                             DISCLAIMER
% 
% ADiMat was prepared as part of an employment at the Institute for Scientific Computing,
% RWTH Aachen University, Germany and at the Institute for Scientific Computing,
% TU Darmstadt, Germany and is provided AS IS. 
% NEITHER THE AUTHOR(S), THE GOVERNMENT OF THE FEDERAL REPUBLIC OF GERMANY
% NOR ANY AGENCY THEREOF, NOR THE RWTH AACHEN UNIVERSITY, NOT THE TU DARMSTADT,
% INCLUDING ANY OF THEIR EMPLOYEES OR OFFICERS, MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
% OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS,
% OR USEFULNESS OF ANY INFORMATION OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE
% WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
%
% Parameters:
%  - dependents=out
%  - independents=U
%  - inputEncoding=ISO-8859-1
%
% Functions in this file: d_LaplaceEqn
%

function [d_out out] = d_LaplaceEqn(d_U, U)
   nx = 8;
   ny = 8;
   dx = 1 / (nx - 1);
   dy = 1 / (ny - 1);
   x = 0:dx:1;
   y = 0:dy:1;
   UW = 0;
   UE = 0;
   US = 0;
   UN = 0;
   u = zeros(nx, ny);
   d_u = d_zeros(u);
   bc(1, :) = UW;
   bc(nx, :) = UE;
   bc(:, 1) = US;
   bc(:, ny) = UN;
   Ex = sparse(2 : nx - 2, 1 : nx - 3, 1, nx - 2, nx - 2);
   Ax = Ex + Ex' - 2*speye(nx - 2);
   Ey = sparse(2 : ny - 2, 1 : ny - 3, 1, ny - 2, ny - 2);
   Ay = Ey + Ey' - 2*speye(ny - 2);
   A = kron(Ay / dy^2, speye(nx - 2)) + kron(speye(ny - 2), Ax / dx^2);
   f = @(x, y) 2.*pi^2.*sin(pi .* x).*sin(pi .* y) + sign(-sin(8 .* pi .* x) .* sin(8 .* pi .* y));
   tmpda1 = x';
   d_S = adimat_opdiff_sum(d_U, d_zeros(f(tmpda1, y)));
   S = f(tmpda1, y) + U;
   tmpva1 = 0;
   S(:, nx) = tmpva1;
   d_S = adimat_opdiff_subsasgn(d_S, struct('type', {'()'}, 'subs', {{':' nx}}), d_zeros(tmpva1));
   tmpva1 = 0;
   S(nx, :) = tmpva1;
   d_S = adimat_opdiff_subsasgn(d_S, struct('type', {'()'}, 'subs', {{nx ':'}}), d_zeros(tmpva1));
   tmpda3 = S;
   tmpda2 = S;
   d_tmpca1 = adimat_opdiff_subsref(d_S, struct('type', '()', 'subs', {{any(tmpda2, 2) any(tmpda3)}}));
   tmpca1 = S(any(tmpda2, 2), any(tmpda3));
   d_S = d_tmpca1;
   S = tmpca1;
   tmpda2 = [];
   d_tmpca1 = adimat_opdiff_sum(d_S, d_zeros(-bc));
   tmpca1 = S - bc;
   [d_S S] = adimat_diff_reshape(d_tmpca1, tmpca1, tmpda2, 1);
   d_tmpca1 = d_S;
   tmpca1 = S;
   d_S = adimat_opdiff_sol_left(A, d_tmpca1, tmpca1);
   S = A \ tmpca1;
   tmpda3 = ny - 2;
   tmpda2 = nx - 2;
   d_tmpca1 = d_S;
   tmpca1 = S;
   [d_S S] = adimat_diff_reshape(d_tmpca1, tmpca1, tmpda2, tmpda3);
   tmpda2 = 2 : ny - 1;
   tmpda1 = 2 : nx - 1;
   d_u = adimat_opdiff_subsasgn(d_u, struct('type', {'()'}, 'subs', {{tmpda1 tmpda2}}), d_S);
   u(tmpda1, tmpda2) = S;
   tmpva1 = UW;
   u(1, :) = tmpva1;
   d_u = adimat_opdiff_subsasgn(d_u, struct('type', {'()'}, 'subs', {{1 ':'}}), d_zeros(tmpva1));
   tmpva1 = US;
   u(:, 1) = tmpva1;
   d_u = adimat_opdiff_subsasgn(d_u, struct('type', {'()'}, 'subs', {{':' 1}}), d_zeros(tmpva1));
   tmpva1 = UE;
   u(nx, :) = tmpva1;
   d_u = adimat_opdiff_subsasgn(d_u, struct('type', {'()'}, 'subs', {{nx ':'}}), d_zeros(tmpva1));
   tmpva1 = UN;
   u(:, ny) = tmpva1;
   d_u = adimat_opdiff_subsasgn(d_u, struct('type', {'()'}, 'subs', {{':' ny}}), d_zeros(tmpva1));
   M = @(x, y) sin(pi .* x).*sin(pi .* y) + sin(8 .* pi .* x).*sin(8 .* pi .* y);
   pw = M(x', y);
   d_tmpca6 = adimat_opdiff_sum(d_u, d_zeros(-pw));
   tmpca6 = u - pw;
   d_tmpca5 = adimat_opdiff_epow_right(d_tmpca6, tmpca6, 2);
   tmpca5 = tmpca6 .^ 2;
   d_tmpca4 = adimat_opdiff_mult_left(0.5, d_tmpca5, tmpca5);
   tmpca4 = 0.5 * tmpca5;
   d_tmpca3 = adimat_opdiff_mult_right(d_tmpca4, tmpca4, dx);
   tmpca3 = tmpca4 * dx;
   d_tmpca2 = adimat_opdiff_mult_right(d_tmpca3, tmpca3, dy);
   tmpca2 = tmpca3 * dy;
   d_tmpca1 = adimat_diff_sum1(d_tmpca2, tmpca2);
   tmpca1 = sum(tmpca2);
   d_out = adimat_diff_sum1(d_tmpca1, tmpca1);
   out = sum(tmpca1);
end