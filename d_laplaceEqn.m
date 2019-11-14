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
% Functions in this file: d_laplaceEqn
%

function [d_out out] = d_laplaceEqn(d_U, U)
% A = 10*rand(size(U));
% b = 10*rand(size(U,1));
% 
% y = A \ (U*b);
% 
% out = sum(sum(y.^2));
%p=1;
%out=p+U;
%UMin= -1;
%UMax= 1;
   pw = 25;
   nx = size(U, 2);
   ny = size(U, 1);
   lambda = 1;
   niter = 1;
   dx = 2 / (nx - 1);
   dy = 2 / (ny - 1);
   x = 0:dx:2;
   y = 0:dy:2;
   M = @(x, y) sin(pi .* x).*sin(pi .* y) + sin(8 .* pi .* x).*sin(8 .* pi .* y);
   P = zeros(ny, nx);
   pn = zeros(ny, nx);
   M = zeros(ny, nx);
   pw = zeros(ny, nx);
   P(:, 1) = 0;
   P(:, nx) = 0;
   P(1, :) = 0;
   P(ny, :) = 0;
   j = 2 : nx - 1;
   i = 2 : ny - 1;
   for it=1 : niter
      pn = P;
      P(i, j) = (dy^2*(pn(i + 1, j) + pn(i - 1, j)) + dx^2*(pn(i, j + 1) + pn(i, j - 1)))/(2 * (dx^2 + dy^2)) - P(i, j);
      P(:, 1) = 0;
      P(:, nx) = 0;
      P(1, :) = 0;
      P(ny, :) = 0;
   end
   j = 1 : nx;
   i = 1 : ny;
   for it=1 : niter
      pw(i, j) = M(i, j);
   end
   d_tmpca2 = adimat_opdiff_sum(adimat_opdiff_mult_right(adimat_opdiff_mult_right(adimat_opdiff_mult_left(lambda / 2, adimat_opdiff_epow_right(adimat_diff_sum1(d_U, U), sum(U), 2), sum(U) .^ 2), lambda/2 * sum(U).^2, dx), lambda/2 * sum(U).^2 * dx, dy), d_zeros(0.5*(P - pw).^2*dx*dy));
   tmpca2 = 0.5*(P - pw).^2*dx*dy + lambda/2*sum(U).^2*dx*dy;
   d_tmpca1 = adimat_diff_sum1(d_tmpca2, tmpca2);
   tmpca1 = sum(tmpca2);
   d_out = adimat_diff_sum1(d_tmpca1, tmpca1);
   out = sum(tmpca1);
end
