function dersv = basisfunder (mu, p, x, knots, nders)

% BASISFUNDER:  B-Spline Basis function derivatives.
%
% Calling Sequence:
% 
%   ders = basisfunder (mu, p, x, knots, nders)
%
%    INPUT:
%   
%      mu    - knot span index (see findspan)
%      p     - degree of curve
%      x     - parametric points
%      knots - knot vector
%      nders - number of derivatives to compute
%
%    OUTPUT:
%   
%      ders - ders(n, i, :) (i-1)-th derivative at n-th point
%   
%    Adapted from Algorithm A2.3 from 'The NURBS BOOK' pg72.
%
% See also: 
%
%    numbasisfun, basisfun, findspan
%
%    Copyright (C) 2009,2011 Rafael Vazquez
%    Copyright (C) 2022 Yannis Voet, Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

x = x(:);
mu = mu(:);

r = min([nders p]);
mu = mu+1;
N = length(x);
D = cell(1, r+1);
B = basisfun(mu-1, x, p-r, knots);
D{r+1} = B;

for k = p-r+1:p

  i1 = mu + ((1-k):0);
  i2 = mu + (1:k);
  t1 = reshape(knots(i1(:)), N, k);
  t2 = reshape(knots(i2(:)), N, k);

  omega = (x-t1)./(t2-t1);
  omega(t2-t1<1e-15) = 0; % "Anything divided by zero is zero" convention
  B = [(1-omega).*B zeros(N,1)]+[zeros(N,1) omega.*B];

  % Save for the next iteration
  D{p-k+1} = B;

  % Compute the derivative
  omega = 1./(t2-t1);
  omega(t2-t1<1e-15) = 0; % "Anything divided by zero is zero" convention
  for j = 1:k-p+r
    D{end-j+1} = [-omega.*D{end-j+1} zeros(N,1)]+[zeros(N,1) omega.*D{end-j+1}];
  end
end

for k = 1:nders-p
  D{r+1+k} = zeros(N, p+1);
end

D = reshape(cell2mat(D), N, p+1, nders+1);
v = 0:nders;
v(v>r) = r;
f = reshape(factorial(p)./factorial(p-v), 1, 1, nders+1);
dersv = f.*D;
dersv = permute(dersv, [1 3 2]);


%!test
%! k    = [0 0 0 0 1 1 1 1];
%! p    = 3;
%! u    = rand (1);
%! i    = findspan (numel(k)-p-2, p, u, k);
%! ders = basisfunder (i, p, u, k, 1);
%! sumders = sum (squeeze(ders), 2);
%! assert (sumders(1), 1, 1e-15);
%! assert (sumders(2:end), 0, 1e-15);

%!test
%! k    = [0 0 0 0 1/3 2/3 1 1 1 1];
%! p    = 3;
%! u    = rand (1);
%! i    = findspan (numel(k)-p-2, p, u, k);
%! ders = basisfunder (i, p, u, k, 7); 
%! sumders = sum (squeeze(ders), 2);
%! assert (sumders(1), 1, 1e-15);
%! assert (sumders(2:end), zeros(rows(squeeze(ders))-1, 1), 1e-13);

%!test
%! k    = [0 0 0 0 1/3 2/3 1 1 1 1];
%! p    = 3;
%! u    = rand (100, 1);
%! i    = findspan (numel(k)-p-2, p, u, k);
%! ders = basisfunder (i, p, u, k, 7);
%! for ii=1:10
%!   sumders = sum (squeeze(ders(ii,:,:)), 2);
%!   assert (sumders(1), 1, 1e-15);
%!   assert (sumders(2:end), zeros(rows(squeeze(ders(ii,:,:)))-1, 1), 1e-13);
%! end
%! assert (ders(:, (p+2):end, :), zeros(numel(u), 8-p-1, p+1), 1e-13)
%! assert (all(all(ders(:, 1, :) <= 1)), true)
