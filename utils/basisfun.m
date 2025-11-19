function B = basisfun (mu, x, p, knots)

% BASISFUN:  Basis function for B-Spline
%
% Calling Sequence:
% 
%   N = basisfun(mu,x,p,knots)
%   
%    INPUT:
%   
%      mu    - knot span  ( from FindSpan() )
%      x     - parametric points
%      p     - spline degree
%      knots - knot sequence
%   
%    OUTPUT:
%   
%      N - Basis functions vector(numel(uv)*(p+1))
%   
%    Adapted from Algorithm A2.2 from 'The NURBS BOOK' pg70.
%
% See also: 
%
%    numbasisfun, basisfunder, findspan
%
% Copyright (C) 2000 Mark Spink
% Copyright (C) 2007 Daniel Claxton
% Copyright (C) 2009 Carlo de Falco
% Copyright (C) 2022 Yannis Voet
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

mu = mu + 1;
% Making sure x and mu are both column vectors
x = x(:);
mu = mu(:);

Np = length(x);
B = ones(Np,1);

for k = 1:p
  i1 = mu + ((1-k):0);
  i2 = mu + (1:k);
  t1 = reshape(knots(i1(:)), Np, k);
  t2 = reshape(knots(i2(:)), Np, k);
  omega = (x-t1)./(t2-t1);
  omega(t2-t1<1e-15) = 0; % "Anything divided by zero is zero" convention
  B = [(1-omega).*B zeros(Np,1)] + [zeros(Np,1) omega.*B];
end

end

%!test
%!  n = 3;
%!  U = [0 0 0 1/2 1 1 1];
%!  p = 2; 
%!  u = linspace (0, 1, 10);
%!  s = findspan (n, p, u, U);
%!  Bref = [1.00000   0.00000   0.00000
%!          0.60494   0.37037   0.02469
%!          0.30864   0.59259   0.09877
%!          0.11111   0.66667   0.22222
%!          0.01235   0.59259   0.39506
%!          0.39506   0.59259   0.01235
%!          0.22222   0.66667   0.11111
%!          0.09877   0.59259   0.30864
%!          0.02469   0.37037   0.60494
%!          0.00000   0.00000   1.00000];
%!  B = basisfun (s, u, p, U);
%!  assert (B, Bref, 1e-5);
