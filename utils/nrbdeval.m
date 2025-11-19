function varargout = nrbdeval (nurbs, dnurbs, varargin)

% NRBDEVAL: Evaluation of the derivative and second derivatives of NURBS curve, surface or volume.
%
%     [pnt, jac] = nrbdeval (crv, dcrv, tt)
%     [pnt, jac] = nrbdeval (srf, dsrf, {tu tv})
%     [pnt, jac] = nrbdeval (vol, dvol, {tu tv tw})
%     [pnt, jac, hess] = nrbdeval (crv, dcrv, dcrv2, tt)
%     [pnt, jac, hess] = nrbdeval (srf, dsrf, dsrf2, {tu tv})
%     [pnt, jac, hess] = nrbdeval (vol, dvol, dvol2, {tu tv tw})
%     [pnt, jac, hess, third_der] = nrbdeval (crv, dcrv, dcrv2, dcrv3, tt)
%     [pnt, jac, hess, third_der] = nrbdeval (srf, dsrf, dsrf2, dsrf3, {tu tv})
%     [pnt, jac, hess, third_der, fourth_der] = nrbdeval (crv, dcrv, dcrv2, dcrv3, dcrv4, tt)
%     [pnt, jac, hess, third_der, fourth_der] = nrbdeval (srf, dsrf, dsrf2, dsrf3, dsrf4, {tu tv})
%
% INPUTS:
%
%   crv,   srf,   vol   - original NURBS curve, surface or volume.
%   dcrv,  dsrf,  dvol  - NURBS derivative representation of crv, srf 
%                          or vol (see nrbderiv2)
%   dcrv2, dsrf2, dvol2 - NURBS second derivative representation of crv,
%                          srf or vol (see nrbderiv2)
%   dcrv3, dsrf3, dvol3 - NURBS third derivative representation of crv,
%                          srf or vol (see nrbderiv)
%   dcrv4, dsrf4, dvol4 - NURBS fourth derivative representation of crv,
%                          srf or vol (see nrbderiv)
%
%   tt     - parametric evaluation points
%            If the nurbs is a surface or a volume then tt is a cell
%            {tu, tv} or {tu, tv, tw} are the parametric coordinates
%
% OUTPUT:
%
%   pnt  - evaluated points.
%   jac  - evaluated first derivatives (Jacobian).
%   hess - evaluated second derivatives (Hessian).
%   third_der - evaluated third derivatives
%   fourth_der - evaluated fourth derivatives
%
% Copyright (C) 2000 Mark Spink 
% Copyright (C) 2010 Carlo de Falco
% Copyright (C) 2010, 2011 Rafael Vazquez
% Copyright (C) 2018 Luca Coradello
% Copyright (C) 2023 Pablo Antolin
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

if (nargin == 3)
  tt = varargin{1};
  dnurbs2 = [];
  dnurbs3 = [];
  dnurbs4 = [];
elseif (nargin == 4)
  dnurbs2 = varargin{1};
  dnurbs3 = [];
  dnurbs4 = [];
  tt = varargin{2};
elseif (nargin == 5)
  dnurbs2 = varargin{1};
  dnurbs3 = varargin{2};       
  dnurbs4 = [];
  tt = varargin{3}; 
elseif (nargin == 6)
  dnurbs2 = varargin{1};
  dnurbs3 = varargin{2};
  dnurbs4 = varargin{3};       
  tt = varargin{4}; 
else 
  error ('nrbrdeval: wrong number of input parameters')
end

if (~isstruct(nurbs))
  error('NURBS representation is not structure!');
end

if (~strcmp(nurbs.form,'B-NURBS'))
  error('Not a recognised NURBS representation');
end

max_der = 1;

if (nargout >= 3)
  if (isempty (dnurbs2))
    warning ('nrbdeval: dnurbs4 missing. The second derivative is not computed');
  else
    max_der = 2;
  end
end

if (nargout >= 4)
  if (isempty (dnurbs3))
    warning ('nrbdeval: dnurbs4 missing. The third derivative is not computed');
  else
    max_der = 3;
  end
end

if (nargout >= 5)
  if (isempty  (dnurbs4))
    warning ('nrbdeval: dnurbs4 missing. The fourth derivative is not computed');
  else
    max_der = 4;
  end
end


% The evaluation of the NURBS can be computed as
% X = P / W, where P and W are the cp and cw just above, respectively.
% Thus, its four first derivatives are computed as:
%
% Xi    = (Pi - X * Wi) / W
%
% Xij   = (Pij - Xi * Wj - Xj * Wi - X * Wij) / W
%
% Xijk  = (Pijk - Xij * Wk   - Xik * Wj -  Xjk * Wi
%               - Xi  * Wjk  - Xj  * Wik - Xk  * Wij
%               - X   * Wijk) / W
%
% Xijkl = (Pijkl - Xijk * Wl   - Xijl * Wk -   Xikl * Wj   - Xjkl * Wi
%                - Xij  * Wkl  - Xik  * Wjl -  Xil  * Wjk  - Xjk  * Wil - Xjl * Wik - Xkl * Wij
%                - Xi   * Wjkl - Xj   * Wikl - Xk   * Wijl - Xl   * Wijk
%                - X    * Wijkl) / W
%
% where the subindices i,j,k,l mean partial derivatives respect to the i (j,k,l) coordinates, respectively.

if (is_rational(nurbs))

  if (iscell(nurbs.knots))
    ndim = size(nurbs.knots,2);
  else
    ndim = 1;
  end

  if (ndim == 1)
    [X, dX, d2X, d3X, d4X] = nrbdeval_crv_rational (nurbs, dnurbs, dnurbs2, dnurbs3, dnurbs4, tt, max_der);
  elseif (ndim == 2)
    [X, dX, d2X, d3X, d4X] = nrbdeval_srf_rational (nurbs, dnurbs, dnurbs2, dnurbs3, dnurbs4, tt, max_der);
  else % if (ndim == 3)
    [X, dX, d2X, d3X, d4X] = nrbdeval_vol_rational (nurbs, dnurbs, dnurbs2, dnurbs3, dnurbs4, tt, max_der);
  end
else
  [X, dX, d2X, d3X, d4X] = nrbdeval_non_rational (nurbs, dnurbs, dnurbs2, dnurbs3, dnurbs4, tt, max_der);
end

varargout{1} = X;
varargout{2} = dX;
if (max_der >= 2)
  varargout{3} = d2X;
end
if (max_der >= 3)
  varargout{4} = d3X;
end
if (max_der == 4)
  varargout{5} = d4X;
end


end

function [X, dX, d2X, d3X, d4X] = nrbdeval_non_rational (nurbs, dnurbs, dnurbs2, dnurbs3, dnurbs4, tt, max_der)

  if (iscell(nurbs.knots))
    ndim = size(nurbs.knots,2);
  else
    ndim = 1;
  end

  [X, ~] = nrbeval(nurbs, tt);

  if (~iscell(dnurbs))
    dnurbs = {dnurbs};
  end

  for i = 1 : ndim
    [dX{i}, ~] = nrbeval(dnurbs{i}, tt);
  end

  d2X = []; d3X = []; d4X = [];

  if (max_der < 2)
      return;
  end

  if (~iscell(dnurbs2))
    dnurbs2 = {dnurbs2};
  end

  for i = 1 : ndim
    for j = 1 : ndim
      [d2X{i,j}, ~] = nrbeval(dnurbs2{i,j}, tt);
   end
  end

  if (max_der < 3)
      return;
  end

  if (~iscell(dnurbs3))
    dnurbs3 = {dnurbs3};
  end

  for i = 1 : ndim
    for j = 1 : ndim
      for k = 1 : ndim
        [d3X{i,j,k}, ~] = nrbeval(dnurbs3{i,j,k}, tt);
      end
    end
  end

  if (max_der < 4)
      return;
  end

  if (~iscell(dnurbs4))
    dnurbs4 = {dnurbs4};
  end

  for i = 1 : ndim
    for j = 1 : ndim
      for k = 1 : ndim
        for l = 1 : ndim
          [d4X{i,j,k,l}, ~] = nrbeval(dnurbs4{i,j,k,l}, tt);
        end
      end
    end
  end
end



function [X, dX, d2X, d3X, d4X] = nrbdeval_crv_rational (nurbs, dnurbs, dnurbs2, dnurbs3, dnurbs4, tt, max_der)

[P, W] = nrbeval(nurbs, tt);
W = W(ones(3,1), :, :, :);
one_div_W = 1.0 ./ W;

X = one_div_W .* P;

clear P;


[dP, dW] = nrbeval(dnurbs, tt);
dW = dW(ones(3,1), :, :, :);

dX{1} = one_div_W .* (dP - X .* dW);

clear dP;

d2X = []; d3X = []; d4X = [];

if (max_der < 2)
    return;
end

[d2P, d2W] = nrbeval(dnurbs2, tt);
d2W = d2W(ones(3,1), :, :, :);

d2X{1} = one_div_W .* (d2P - dX{1} .* dW - dX{1} .* dW - X .* d2W);

if (max_der < 3)
  return;
end

clear d2P;


[d3P, d3W] = nrbeval(dnurbs3, tt);
d3W = d3W(ones(3,1), :, :, :);

d3X{1} = one_div_W .* (d3P - 3 * d2X{1} .* dW - 3 * dX{1} .* d2W - X .* d3W);


if (max_der < 4)
  return;
end


[d4P, d4W] = nrbeval(dnurbs4, tt);
d4W = d4W(ones(3,1), :, :, :);

d4X{1} = one_div_W .* (d4P - 4 * d3X{1} .* dW - 6 * d2X{1} .* d2W - 4 * dX{1} .* d3W - X .* d4W);


end

function [X, dX, d2X, d3X, d4X] = nrbdeval_srf_rational (nurbs, dnurbs, dnurbs2, dnurbs3, dnurbs4, tt, max_der)

[P, W] = nrbeval(nurbs, tt);
W = W(ones(3,1), :, :, :);
one_div_W = 1.0 ./ W;

X = one_div_W .* P;

clear P;

dP = cell(1, 2); dW = cell(1, 2);
[dP{1}, dW{1}] = nrbeval(dnurbs{1}, tt);
[dP{2}, dW{2}] = nrbeval(dnurbs{2}, tt);
dW{1} = dW{1}(ones(3,1), :, :, :);
dW{2} = dW{2}(ones(3,1), :, :, :);

dX{1} = one_div_W .* (dP{1} - X .* dW{1});
dX{2} = one_div_W .* (dP{2} - X .* dW{2});

d2X = []; d3X = []; d4X = [];

if (max_der < 2)
    return;
end

clear dP;

d2P = cell(2, 2); d2W = cell(2, 2);
[d2P{1,1}, d2W{1,1}] = nrbeval(dnurbs2{1,1}, tt);
[d2P{1,2}, d2W{1,2}] = nrbeval(dnurbs2{1,2}, tt);
[d2P{2,2}, d2W{2,2}] = nrbeval(dnurbs2{2,2}, tt);
d2W{1,1} = d2W{1,1}(ones(3,1), :, :, :);
d2W{1,2} = d2W{1,2}(ones(3,1), :, :, :);
d2W{2,2} = d2W{2,2}(ones(3,1), :, :, :);

d2X{1, 1} = one_div_W .* (d2P{1, 1} - dX{1} .* dW{1} - dX{1} .* dW{1} - X .* d2W{1, 1});
d2X{1, 2} = one_div_W .* (d2P{1, 2} - dX{1} .* dW{2} - dX{2} .* dW{1} - X .* d2W{1, 2});
d2X{2, 2} = one_div_W .* (d2P{2, 2} - dX{2} .* dW{2} - dX{2} .* dW{2} - X .* d2W{2, 2});

d2X{2, 1} = d2X{1, 2};

if (max_der < 3)
    return;
end

clear d2P;

d3P = cell(2, 2, 2); d3W = cell(2, 2, 2);
[d3P{1,1,1}, d3W{1,1,1}] = nrbeval(dnurbs3{1,1,1}, tt);
[d3P{1,1,2}, d3W{1,1,2}] = nrbeval(dnurbs3{1,1,2}, tt);
[d3P{1,2,2}, d3W{1,2,2}] = nrbeval(dnurbs3{1,2,2}, tt);
[d3P{2,2,2}, d3W{2,2,2}] = nrbeval(dnurbs3{2,2,2}, tt);
d3W{1,1,1} = d3W{1,1,1}(ones(3,1), :, :, :);
d3W{1,1,2} = d3W{1,1,2}(ones(3,1), :, :, :);
d3W{1,2,2} = d3W{1,2,2}(ones(3,1), :, :, :);
d3W{2,2,2} = d3W{2,2,2}(ones(3,1), :, :, :);

d3X{1, 1, 1} = one_div_W .* (d3P{1, 1, 1} - 3 * d2X{1, 1} .* dW{1} - 3 * dX{1} .* d2W{1, 1} - X .* d3W{1, 1, 1});
d3X{1, 1, 2} = one_div_W .* (d3P{1, 1, 2} - d2X{1, 1} .* dW{2} - 2 * d2X{1, 2} .* dW{1} - 2 * dX{1} .* d2W{1, 2} - dX{2} .* d2W{1, 1} - X .* d3W{1, 1, 2});
d3X{1, 2, 2} = one_div_W .* (d3P{1, 2, 2} - 2 * d2X{1, 2} .* dW{2} - d2X{2, 2} .* dW{1} - dX{1} .* d2W{2, 2} - 2 * dX{2} .* d2W{1, 2} - X .* d3W{1, 2, 2});
d3X{2, 2, 2} = one_div_W .* (d3P{2, 2, 2} - 3 * d2X{2, 2} .* dW{2} - 3 * dX{2} .* d2W{2, 2} - X .* d3W{2, 2, 2});

d3X{2, 1, 1} = d3X{1, 1, 2};
d3X{1, 2, 1} = d3X{1, 1, 2};
d3X{2, 2, 1} = d3X{1, 2, 2};
d3X{2, 1, 2} = d3X{1, 2, 2};

if (max_der < 4)
    return;
end

clear d3P;


d4P = cell(2, 2, 2, 2); d4W = cell(2, 2, 2, 2);
[d4P{1,1,1,1}, d4W{1,1,1,1}] = nrbeval(dnurbs4{1,1,1,1}, tt);
[d4P{1,1,1,2}, d4W{1,1,1,2}] = nrbeval(dnurbs4{1,1,1,2}, tt);
[d4P{1,1,2,2}, d4W{1,1,2,2}] = nrbeval(dnurbs4{1,1,2,2}, tt);
[d4P{1,2,2,2}, d4W{1,2,2,2}] = nrbeval(dnurbs4{1,2,2,2}, tt);
[d4P{2,2,2,2}, d4W{2,2,2,2}] = nrbeval(dnurbs4{2,2,2,2}, tt);
d4W{1,1,1,1} = d4W{1,1,1,1}(ones(3,1), :, :, :);
d4W{1,1,1,2} = d4W{1,1,1,2}(ones(3,1), :, :, :);
d4W{1,1,2,2} = d4W{1,1,2,2}(ones(3,1), :, :, :);
d4W{1,2,2,2} = d4W{1,2,2,2}(ones(3,1), :, :, :);
d4W{2,2,2,2} = d4W{2,2,2,2}(ones(3,1), :, :, :);

d4X{1, 1, 1, 1} = one_div_W .* (d4P{1, 1, 1, 1} - 4 * d3X{1, 1, 1} .* dW{1} - 6 * d2X{1, 1} .* d2W{1, 1} - 4 * dX{1} .* d3W{1, 1, 1} - X .* d4W{1, 1, 1, 1});
d4X{1, 1, 1, 2} = one_div_W .* (d4P{1, 1, 1, 2} - d3X{1, 1, 1} .* dW{2} - 3 * d3X{1, 1, 2} .* dW{1} - 3 * d2X{1, 1} .* d2W{1, 2} - 3 * d2X{1, 2} .* d2W{1, 1} - 3 * dX{1} .* d3W{1, 1, 2} - dX{2} .* d3W{1, 1, 1} - X .* d4W{1, 1, 1, 2});
d4X{1, 1, 2, 2} = one_div_W .* (d4P{1, 1, 2, 2} - 2 * d3X{1, 1, 2} .* dW{2} - 2 * d3X{1, 2, 2} .* dW{1} - d2X{1, 1} .* d2W{2, 2} - 4 * d2X{1, 2} .* d2W{1, 2} - d2X{2, 2} .* d2W{1, 1} - 2 * dX{1} .* d3W{1, 2, 2} - 2 * dX{2} .* d3W{1, 1, 2} - X .* d4W{1, 1, 2, 2});
d4X{1, 2, 2, 2} = one_div_W .* (d4P{1, 2, 2, 2} - 3 * d3X{1, 2, 2} .* dW{2} - d3X{2, 2, 2} .* dW{1} - 3 * d2X{1, 2} .* d2W{2, 2} - 3 * d2X{2, 2} .* d2W{1, 2} - dX{1} .* d3W{2, 2, 2} - 3 * dX{2} .* d3W{1, 2, 2} - X .* d4W{1, 2, 2, 2});
d4X{2, 2, 2, 2} = one_div_W .* (d4P{2, 2, 2, 2} - 4 * d3X{2, 2, 2} .* dW{2} - 6 * d2X{2, 2} .* d2W{2, 2} - 4 * dX{2} .* d3W{2, 2, 2} - X .* d4W{2, 2, 2, 2});

d4X{2, 1, 1, 1} = d4X{1, 1, 1, 2};
d4X{1, 2, 1, 1} = d4X{1, 1, 1, 2};
d4X{2, 2, 1, 1} = d4X{1, 1, 2, 2};
d4X{1, 1, 2, 1} = d4X{1, 1, 1, 2};
d4X{2, 1, 2, 1} = d4X{1, 1, 2, 2};
d4X{1, 2, 2, 1} = d4X{1, 1, 2, 2};
d4X{2, 2, 2, 1} = d4X{1, 2, 2, 2};
d4X{2, 1, 1, 2} = d4X{1, 1, 2, 2};
d4X{1, 2, 1, 2} = d4X{1, 1, 2, 2};
d4X{2, 2, 1, 2} = d4X{1, 2, 2, 2};
d4X{2, 1, 2, 2} = d4X{1, 2, 2, 2};

end

function [X, dX, d2X, d3X, d4X] = nrbdeval_vol_rational (nurbs, dnurbs, dnurbs2, dnurbs3, dnurbs4, tt, max_der)

[P, W] = nrbeval(nurbs, tt);
W = W(ones(3,1), :, :, :);
one_div_W = 1.0 ./ W;

X = one_div_W .* P;


dP = cell(1, 3); dW = cell(1, 3);
[dP{1}, dW{1}] = nrbeval(dnurbs{1}, tt);
[dP{2}, dW{2}] = nrbeval(dnurbs{2}, tt);
[dP{3}, dW{3}] = nrbeval(dnurbs{3}, tt);
dW{1} = dW{1}(ones(3,1), :, :, :);
dW{2} = dW{2}(ones(3,1), :, :, :);
dW{3} = dW{3}(ones(3,1), :, :, :);

dX{1} = one_div_W .* (dP{1} - X .* dW{1});
dX{2} = one_div_W .* (dP{2} - X .* dW{2});
dX{3} = one_div_W .* (dP{3} - X .* dW{3});


d2X = []; d3X = []; d4X = [];

if (max_der < 2)
    return;
end

clear dP;

d2P = cell(3, 3); d2W = cell(3, 3);
[d2P{1,1}, d2W{1,1}] = nrbeval(dnurbs2{1,1}, tt);
[d2P{1,2}, d2W{1,2}] = nrbeval(dnurbs2{1,2}, tt);
[d2P{1,3}, d2W{1,3}] = nrbeval(dnurbs2{1,3}, tt);
[d2P{2,2}, d2W{2,2}] = nrbeval(dnurbs2{2,2}, tt);
[d2P{2,3}, d2W{2,3}] = nrbeval(dnurbs2{2,3}, tt);
[d2P{3,3}, d2W{3,3}] = nrbeval(dnurbs2{3,3}, tt);
d2W{1,1} = d2W{1,1}(ones(3,1), :, :, :);
d2W{1,2} = d2W{1,2}(ones(3,1), :, :, :);
d2W{1,3} = d2W{1,3}(ones(3,1), :, :, :);
d2W{2,2} = d2W{2,2}(ones(3,1), :, :, :);
d2W{2,3} = d2W{2,3}(ones(3,1), :, :, :);
d2W{3,3} = d2W{3,3}(ones(3,1), :, :, :);

d2X{1, 1} = one_div_W .* (d2P{1, 1} - dX{1} .* dW{1} - dX{1} .* dW{1} - X .* d2W{1, 1});
d2X{1, 2} = one_div_W .* (d2P{1, 2} - dX{1} .* dW{2} - dX{2} .* dW{1} - X .* d2W{1, 2});
d2X{1, 3} = one_div_W .* (d2P{1, 3} - dX{1} .* dW{3} - dX{3} .* dW{1} - X .* d2W{1, 3});
d2X{2, 2} = one_div_W .* (d2P{2, 2} - dX{2} .* dW{2} - dX{2} .* dW{2} - X .* d2W{2, 2});
d2X{2, 3} = one_div_W .* (d2P{2, 3} - dX{2} .* dW{3} - dX{3} .* dW{2} - X .* d2W{2, 3});
d2X{3, 3} = one_div_W .* (d2P{3, 3} - dX{3} .* dW{3} - dX{3} .* dW{3} - X .* d2W{3, 3});

d2X{2, 1} = d2X{1, 2};
d2X{3, 1} = d2X{1, 3};
d2X{3, 2} = d2X{2, 3};

if (max_der < 3)
    return;
end

clear d2P;


d3P = cell(3, 3, 3); d3W = cell(3, 3, 3);
[d3P{1,1,1}, d3W{1,1,1}] = nrbeval(dnurbs3{1,1,1}, tt);
[d3P{1,1,2}, d3W{1,1,2}] = nrbeval(dnurbs3{1,1,2}, tt);
[d3P{1,1,3}, d3W{1,1,3}] = nrbeval(dnurbs3{1,1,3}, tt);
[d3P{1,2,2}, d3W{1,2,2}] = nrbeval(dnurbs3{1,2,2}, tt);
[d3P{1,2,3}, d3W{1,2,3}] = nrbeval(dnurbs3{1,2,3}, tt);
[d3P{1,3,3}, d3W{1,3,3}] = nrbeval(dnurbs3{1,3,3}, tt);
[d3P{2,2,2}, d3W{2,2,2}] = nrbeval(dnurbs3{2,2,2}, tt);
[d3P{2,2,3}, d3W{2,2,3}] = nrbeval(dnurbs3{2,2,3}, tt);
[d3P{2,3,3}, d3W{2,3,3}] = nrbeval(dnurbs3{2,3,3}, tt);
[d3P{3,3,3}, d3W{3,3,3}] = nrbeval(dnurbs3{3,3,3}, tt);
d3W{1,1,1} = d3W{1,1,1}(ones(3,1), :, :, :);
d3W{1,1,2} = d3W{1,1,2}(ones(3,1), :, :, :);
d3W{1,1,3} = d3W{1,1,3}(ones(3,1), :, :, :);
d3W{1,2,2} = d3W{1,2,2}(ones(3,1), :, :, :);
d3W{1,2,3} = d3W{1,2,3}(ones(3,1), :, :, :);
d3W{1,3,3} = d3W{1,3,3}(ones(3,1), :, :, :);
d3W{2,2,2} = d3W{2,2,2}(ones(3,1), :, :, :);
d3W{2,2,3} = d3W{2,2,3}(ones(3,1), :, :, :);
d3W{2,3,3} = d3W{2,3,3}(ones(3,1), :, :, :);
d3W{3,3,3} = d3W{3,3,3}(ones(3,1), :, :, :);

d3X{1, 1, 1} = one_div_W .* (d3P{1, 1, 1} - 3 * d2X{1, 1} .* dW{1} - 3 * dX{1} .* d2W{1, 1} - X .* d3W{1, 1, 1});
d3X{1, 1, 2} = one_div_W .* (d3P{1, 1, 2} - d2X{1, 1} .* dW{2} - 2 * d2X{1, 2} .* dW{1} - 2 * dX{1} .* d2W{1, 2} - dX{2} .* d2W{1, 1} - X .* d3W{1, 1, 2});
d3X{1, 1, 3} = one_div_W .* (d3P{1, 1, 3} - d2X{1, 1} .* dW{3} - 2 * d2X{1, 3} .* dW{1} - 2 * dX{1} .* d2W{1, 3} - dX{3} .* d2W{1, 1} - X .* d3W{1, 1, 3});
d3X{1, 2, 2} = one_div_W .* (d3P{1, 2, 2} - 2 * d2X{1, 2} .* dW{2} - d2X{2, 2} .* dW{1} - dX{1} .* d2W{2, 2} - 2 * dX{2} .* d2W{1, 2} - X .* d3W{1, 2, 2});
d3X{1, 2, 3} = one_div_W .* (d3P{1, 2, 3} - d2X{1, 2} .* dW{3} - d2X{1, 3} .* dW{2} - d2X{2, 3} .* dW{1} - dX{1} .* d2W{2, 3} - dX{2} .* d2W{1, 3} - dX{3} .* d2W{1, 2} - X .* d3W{1, 2, 3});
d3X{1, 3, 3} = one_div_W .* (d3P{1, 3, 3} - 2 * d2X{1, 3} .* dW{3} - d2X{3, 3} .* dW{1} - dX{1} .* d2W{3, 3} - 2 * dX{3} .* d2W{1, 3} - X .* d3W{1, 3, 3});
d3X{2, 2, 2} = one_div_W .* (d3P{2, 2, 2} - 3 * d2X{2, 2} .* dW{2} - 3 * dX{2} .* d2W{2, 2} - X .* d3W{2, 2, 2});
d3X{2, 2, 3} = one_div_W .* (d3P{2, 2, 3} - d2X{2, 2} .* dW{3} - 2 * d2X{2, 3} .* dW{2} - 2 * dX{2} .* d2W{2, 3} - dX{3} .* d2W{2, 2} - X .* d3W{2, 2, 3});
d3X{2, 3, 3} = one_div_W .* (d3P{2, 3, 3} - 2 * d2X{2, 3} .* dW{3} - d2X{3, 3} .* dW{2} - dX{2} .* d2W{3, 3} - 2 * dX{3} .* d2W{2, 3} - X .* d3W{2, 3, 3});
d3X{3, 3, 3} = one_div_W .* (d3P{3, 3, 3} - 3 * d2X{3, 3} .* dW{3} - 3 * dX{3} .* d2W{3, 3} - X .* d3W{3, 3, 3});

d3X{2, 1, 1} = d3X{1, 1, 2};
d3X{3, 1, 1} = d3X{1, 1, 3};
d3X{1, 2, 1} = d3X{1, 1, 2};
d3X{2, 2, 1} = d3X{1, 2, 2};
d3X{3, 2, 1} = d3X{1, 2, 3};
d3X{1, 3, 1} = d3X{1, 1, 3};
d3X{2, 3, 1} = d3X{1, 2, 3};
d3X{3, 3, 1} = d3X{1, 3, 3};
d3X{2, 1, 2} = d3X{1, 2, 2};
d3X{3, 1, 2} = d3X{1, 2, 3};
d3X{3, 2, 2} = d3X{2, 2, 3};
d3X{1, 3, 2} = d3X{1, 2, 3};
d3X{2, 3, 2} = d3X{2, 2, 3};
d3X{3, 3, 2} = d3X{2, 3, 3};
d3X{2, 1, 3} = d3X{1, 2, 3};
d3X{3, 1, 3} = d3X{1, 3, 3};
d3X{3, 2, 3} = d3X{2, 3, 3};

if (max_der < 4)
    return;
end


clear d3P;


d4P = cell(3, 3, 3, 3); d4W = cell(3, 3, 3, 3);
[d4P{1,1,1,1}, d4W{1,1,1,1}] = nrbeval(dnurbs4{1,1,1,1}, tt);
[d4P{1,1,1,2}, d4W{1,1,1,2}] = nrbeval(dnurbs4{1,1,1,2}, tt);
[d4P{1,1,1,3}, d4W{1,1,1,3}] = nrbeval(dnurbs4{1,1,1,3}, tt);
[d4P{1,1,2,2}, d4W{1,1,2,2}] = nrbeval(dnurbs4{1,1,2,2}, tt);
[d4P{1,1,2,3}, d4W{1,1,2,3}] = nrbeval(dnurbs4{1,1,2,3}, tt);
[d4P{1,1,3,3}, d4W{1,1,3,3}] = nrbeval(dnurbs4{1,1,3,3}, tt);
[d4P{1,2,2,2}, d4W{1,2,2,2}] = nrbeval(dnurbs4{1,2,2,2}, tt);
[d4P{1,2,2,3}, d4W{1,2,2,3}] = nrbeval(dnurbs4{1,2,2,3}, tt);
[d4P{1,2,3,3}, d4W{1,2,3,3}] = nrbeval(dnurbs4{1,2,3,3}, tt);
[d4P{1,3,3,3}, d4W{1,3,3,3}] = nrbeval(dnurbs4{1,3,3,3}, tt);
[d4P{2,2,2,2}, d4W{2,2,2,2}] = nrbeval(dnurbs4{2,2,2,2}, tt);
[d4P{2,2,2,3}, d4W{2,2,2,3}] = nrbeval(dnurbs4{2,2,2,3}, tt);
[d4P{2,2,3,3}, d4W{2,2,3,3}] = nrbeval(dnurbs4{2,2,3,3}, tt);
[d4P{2,3,3,3}, d4W{2,3,3,3}] = nrbeval(dnurbs4{2,3,3,3}, tt);
[d4P{3,3,3,3}, d4W{3,3,3,3}] = nrbeval(dnurbs4{3,3,3,3}, tt);
d4W{1,1,1,1} = d4W{1,1,1,1}(ones(3,1), :, :, :);
d4W{1,1,1,2} = d4W{1,1,1,2}(ones(3,1), :, :, :);
d4W{1,1,1,3} = d4W{1,1,1,3}(ones(3,1), :, :, :);
d4W{1,1,2,2} = d4W{1,1,2,2}(ones(3,1), :, :, :);
d4W{1,1,2,3} = d4W{1,1,2,3}(ones(3,1), :, :, :);
d4W{1,1,3,3} = d4W{1,1,3,3}(ones(3,1), :, :, :);
d4W{1,2,2,2} = d4W{1,2,2,2}(ones(3,1), :, :, :);
d4W{1,2,2,3} = d4W{1,2,2,3}(ones(3,1), :, :, :);
d4W{1,2,3,3} = d4W{1,2,3,3}(ones(3,1), :, :, :);
d4W{1,3,3,3} = d4W{1,3,3,3}(ones(3,1), :, :, :);
d4W{2,2,2,2} = d4W{2,2,2,2}(ones(3,1), :, :, :);
d4W{2,2,2,3} = d4W{2,2,2,3}(ones(3,1), :, :, :);
d4W{2,2,3,3} = d4W{2,2,3,3}(ones(3,1), :, :, :);
d4W{2,3,3,3} = d4W{2,3,3,3}(ones(3,1), :, :, :);
d4W{3,3,3,3} = d4W{3,3,3,3}(ones(3,1), :, :, :);

d4X{1, 1, 1, 1} = one_div_W .* (d4P{1, 1, 1, 1} - 4 * d3X{1, 1, 1} .* dW{1} - 6 * d2X{1, 1} .* d2W{1, 1} - 4 * dX{1} .* d3W{1, 1, 1} - X .* d4W{1, 1, 1, 1});
d4X{1, 1, 1, 2} = one_div_W .* (d4P{1, 1, 1, 2} - d3X{1, 1, 1} .* dW{2} - 3 * d3X{1, 1, 2} .* dW{1} - 3 * d2X{1, 1} .* d2W{1, 2} - 3 * d2X{1, 2} .* d2W{1, 1} - 3 * dX{1} .* d3W{1, 1, 2} - dX{2} .* d3W{1, 1, 1} - X .* d4W{1, 1, 1, 2});
d4X{1, 1, 1, 3} = one_div_W .* (d4P{1, 1, 1, 3} - d3X{1, 1, 1} .* dW{3} - 3 * d3X{1, 1, 3} .* dW{1} - 3 * d2X{1, 1} .* d2W{1, 3} - 3 * d2X{1, 3} .* d2W{1, 1} - 3 * dX{1} .* d3W{1, 1, 3} - dX{3} .* d3W{1, 1, 1} - X .* d4W{1, 1, 1, 3});
d4X{1, 1, 2, 2} = one_div_W .* (d4P{1, 1, 2, 2} - 2 * d3X{1, 1, 2} .* dW{2} - 2 * d3X{1, 2, 2} .* dW{1} - d2X{1, 1} .* d2W{2, 2} - 4 * d2X{1, 2} .* d2W{1, 2} - d2X{2, 2} .* d2W{1, 1} - 2 * dX{1} .* d3W{1, 2, 2} - 2 * dX{2} .* d3W{1, 1, 2} - X .* d4W{1, 1, 2, 2});
d4X{1, 1, 2, 3} = one_div_W .* (d4P{1, 1, 2, 3} - d3X{1, 1, 2} .* dW{3} - d3X{1, 1, 3} .* dW{2} - 2 * d3X{1, 2, 3} .* dW{1} - d2X{1, 1} .* d2W{2, 3} - 2 * d2X{1, 2} .* d2W{1, 3} - 2 * d2X{1, 3} .* d2W{1, 2} - d2X{2, 3} .* d2W{1, 1} - 2 * dX{1} .* d3W{1, 2, 3} - dX{2} .* d3W{1, 1, 3} - dX{3} .* d3W{1, 1, 2} - X .* d4W{1, 1, 2, 3});
d4X{1, 1, 3, 3} = one_div_W .* (d4P{1, 1, 3, 3} - 2 * d3X{1, 1, 3} .* dW{3} - 2 * d3X{1, 3, 3} .* dW{1} - d2X{1, 1} .* d2W{3, 3} - 4 * d2X{1, 3} .* d2W{1, 3} - d2X{3, 3} .* d2W{1, 1} - 2 * dX{1} .* d3W{1, 3, 3} - 2 * dX{3} .* d3W{1, 1, 3} - X .* d4W{1, 1, 3, 3});
d4X{1, 2, 2, 2} = one_div_W .* (d4P{1, 2, 2, 2} - 3 * d3X{1, 2, 2} .* dW{2} - d3X{2, 2, 2} .* dW{1} - 3 * d2X{1, 2} .* d2W{2, 2} - 3 * d2X{2, 2} .* d2W{1, 2} - dX{1} .* d3W{2, 2, 2} - 3 * dX{2} .* d3W{1, 2, 2} - X .* d4W{1, 2, 2, 2});
d4X{1, 2, 2, 3} = one_div_W .* (d4P{1, 2, 2, 3} - d3X{1, 2, 2} .* dW{3} - 2 * d3X{1, 2, 3} .* dW{2} - d3X{2, 2, 3} .* dW{1} - 2 * d2X{1, 2} .* d2W{2, 3} - d2X{1, 3} .* d2W{2, 2} - d2X{2, 2} .* d2W{1, 3} - 2 * d2X{2, 3} .* d2W{1, 2} - dX{1} .* d3W{2, 2, 3} - 2 * dX{2} .* d3W{1, 2, 3} - dX{3} .* d3W{1, 2, 2} - X .* d4W{1, 2, 2, 3});
d4X{1, 2, 3, 3} = one_div_W .* (d4P{1, 2, 3, 3} - 2 * d3X{1, 2, 3} .* dW{3} - d3X{1, 3, 3} .* dW{2} - d3X{2, 3, 3} .* dW{1} - d2X{1, 2} .* d2W{3, 3} - 2 * d2X{1, 3} .* d2W{2, 3} - 2 * d2X{2, 3} .* d2W{1, 3} - d2X{3, 3} .* d2W{1, 2} - dX{1} .* d3W{2, 3, 3} - dX{2} .* d3W{1, 3, 3} - 2 * dX{3} .* d3W{1, 2, 3} - X .* d4W{1, 2, 3, 3});
d4X{1, 3, 3, 3} = one_div_W .* (d4P{1, 3, 3, 3} - 3 * d3X{1, 3, 3} .* dW{3} - d3X{3, 3, 3} .* dW{1} - 3 * d2X{1, 3} .* d2W{3, 3} - 3 * d2X{3, 3} .* d2W{1, 3} - dX{1} .* d3W{3, 3, 3} - 3 * dX{3} .* d3W{1, 3, 3} - X .* d4W{1, 3, 3, 3});
d4X{2, 2, 2, 2} = one_div_W .* (d4P{2, 2, 2, 2} - 4 * d3X{2, 2, 2} .* dW{2} - 6 * d2X{2, 2} .* d2W{2, 2} - 4 * dX{2} .* d3W{2, 2, 2} - X .* d4W{2, 2, 2, 2});
d4X{2, 2, 2, 3} = one_div_W .* (d4P{2, 2, 2, 3} - d3X{2, 2, 2} .* dW{3} - 3 * d3X{2, 2, 3} .* dW{2} - 3 * d2X{2, 2} .* d2W{2, 3} - 3 * d2X{2, 3} .* d2W{2, 2} - 3 * dX{2} .* d3W{2, 2, 3} - dX{3} .* d3W{2, 2, 2} - X .* d4W{2, 2, 2, 3});
d4X{2, 2, 3, 3} = one_div_W .* (d4P{2, 2, 3, 3} - 2 * d3X{2, 2, 3} .* dW{3} - 2 * d3X{2, 3, 3} .* dW{2} - d2X{2, 2} .* d2W{3, 3} - 4 * d2X{2, 3} .* d2W{2, 3} - d2X{3, 3} .* d2W{2, 2} - 2 * dX{2} .* d3W{2, 3, 3} - 2 * dX{3} .* d3W{2, 2, 3} - X .* d4W{2, 2, 3, 3});
d4X{2, 3, 3, 3} = one_div_W .* (d4P{2, 3, 3, 3} - 3 * d3X{2, 3, 3} .* dW{3} - d3X{3, 3, 3} .* dW{2} - 3 * d2X{2, 3} .* d2W{3, 3} - 3 * d2X{3, 3} .* d2W{2, 3} - dX{2} .* d3W{3, 3, 3} - 3 * dX{3} .* d3W{2, 3, 3} - X .* d4W{2, 3, 3, 3});
d4X{3, 3, 3, 3} = one_div_W .* (d4P{3, 3, 3, 3} - 4 * d3X{3, 3, 3} .* dW{3} - 6 * d2X{3, 3} .* d2W{3, 3} - 4 * dX{3} .* d3W{3, 3, 3} - X .* d4W{3, 3, 3, 3});

d4X{2, 1, 1, 1} = d4X{1, 1, 1, 2};
d4X{3, 1, 1, 1} = d4X{1, 1, 1, 3};
d4X{1, 2, 1, 1} = d4X{1, 1, 1, 2};
d4X{2, 2, 1, 1} = d4X{1, 1, 2, 2};
d4X{3, 2, 1, 1} = d4X{1, 1, 2, 3};
d4X{1, 3, 1, 1} = d4X{1, 1, 1, 3};
d4X{2, 3, 1, 1} = d4X{1, 1, 2, 3};
d4X{3, 3, 1, 1} = d4X{1, 1, 3, 3};
d4X{1, 1, 2, 1} = d4X{1, 1, 1, 2};
d4X{2, 1, 2, 1} = d4X{1, 1, 2, 2};
d4X{3, 1, 2, 1} = d4X{1, 1, 2, 3};
d4X{1, 2, 2, 1} = d4X{1, 1, 2, 2};
d4X{2, 2, 2, 1} = d4X{1, 2, 2, 2};
d4X{3, 2, 2, 1} = d4X{1, 2, 2, 3};
d4X{1, 3, 2, 1} = d4X{1, 1, 2, 3};
d4X{2, 3, 2, 1} = d4X{1, 2, 2, 3};
d4X{3, 3, 2, 1} = d4X{1, 2, 3, 3};
d4X{1, 1, 3, 1} = d4X{1, 1, 1, 3};
d4X{2, 1, 3, 1} = d4X{1, 1, 2, 3};
d4X{3, 1, 3, 1} = d4X{1, 1, 3, 3};
d4X{1, 2, 3, 1} = d4X{1, 1, 2, 3};
d4X{2, 2, 3, 1} = d4X{1, 2, 2, 3};
d4X{3, 2, 3, 1} = d4X{1, 2, 3, 3};
d4X{1, 3, 3, 1} = d4X{1, 1, 3, 3};
d4X{2, 3, 3, 1} = d4X{1, 2, 3, 3};
d4X{3, 3, 3, 1} = d4X{1, 3, 3, 3};
d4X{2, 1, 1, 2} = d4X{1, 1, 2, 2};
d4X{3, 1, 1, 2} = d4X{1, 1, 2, 3};
d4X{1, 2, 1, 2} = d4X{1, 1, 2, 2};
d4X{2, 2, 1, 2} = d4X{1, 2, 2, 2};
d4X{3, 2, 1, 2} = d4X{1, 2, 2, 3};
d4X{1, 3, 1, 2} = d4X{1, 1, 2, 3};
d4X{2, 3, 1, 2} = d4X{1, 2, 2, 3};
d4X{3, 3, 1, 2} = d4X{1, 2, 3, 3};
d4X{2, 1, 2, 2} = d4X{1, 2, 2, 2};
d4X{3, 1, 2, 2} = d4X{1, 2, 2, 3};
d4X{3, 2, 2, 2} = d4X{2, 2, 2, 3};
d4X{1, 3, 2, 2} = d4X{1, 2, 2, 3};
d4X{2, 3, 2, 2} = d4X{2, 2, 2, 3};
d4X{3, 3, 2, 2} = d4X{2, 2, 3, 3};
d4X{1, 1, 3, 2} = d4X{1, 1, 2, 3};
d4X{2, 1, 3, 2} = d4X{1, 2, 2, 3};
d4X{3, 1, 3, 2} = d4X{1, 2, 3, 3};
d4X{1, 2, 3, 2} = d4X{1, 2, 2, 3};
d4X{2, 2, 3, 2} = d4X{2, 2, 2, 3};
d4X{3, 2, 3, 2} = d4X{2, 2, 3, 3};
d4X{1, 3, 3, 2} = d4X{1, 2, 3, 3};
d4X{2, 3, 3, 2} = d4X{2, 2, 3, 3};
d4X{3, 3, 3, 2} = d4X{2, 3, 3, 3};
d4X{2, 1, 1, 3} = d4X{1, 1, 2, 3};
d4X{3, 1, 1, 3} = d4X{1, 1, 3, 3};
d4X{1, 2, 1, 3} = d4X{1, 1, 2, 3};
d4X{2, 2, 1, 3} = d4X{1, 2, 2, 3};
d4X{3, 2, 1, 3} = d4X{1, 2, 3, 3};
d4X{1, 3, 1, 3} = d4X{1, 1, 3, 3};
d4X{2, 3, 1, 3} = d4X{1, 2, 3, 3};
d4X{3, 3, 1, 3} = d4X{1, 3, 3, 3};
d4X{2, 1, 2, 3} = d4X{1, 2, 2, 3};
d4X{3, 1, 2, 3} = d4X{1, 2, 3, 3};
d4X{3, 2, 2, 3} = d4X{2, 2, 3, 3};
d4X{1, 3, 2, 3} = d4X{1, 2, 3, 3};
d4X{2, 3, 2, 3} = d4X{2, 2, 3, 3};
d4X{3, 3, 2, 3} = d4X{2, 3, 3, 3};
d4X{2, 1, 3, 3} = d4X{1, 2, 3, 3};
d4X{3, 1, 3, 3} = d4X{1, 3, 3, 3};
d4X{3, 2, 3, 3} = d4X{2, 3, 3, 3};


end

function rational = is_rational(nurbs)
if (size(nurbs.coefs, 1) < 4)
    rational = false;
else
    tolerance = 1.0e-15;
    rational = any(abs(nurbs.coefs(4, :) - 1) > tolerance);
end
end

%!demo
%! crv = nrbtestcrv;
%! nrbplot(crv,48);
%! title('First derivatives along a test curve.');
%! 
%! tt = linspace(0.0,1.0,9);
%! 
%! dcrv = nrbderiv(crv);
%! 
%! [p1, dp] = nrbdeval(crv,dcrv,tt);
%! 
%! p2 = vecnormalize(dp);
%! 
%! hold on;
%! plot(p1(1,:),p1(2,:),'ro');
%! h = quiver(p1(1,:),p1(2,:),p2(1,:),p2(2,:),0);
%! set(h,'Color','black');
%! hold off;

%!demo
%! srf = nrbtestsrf;
%! p = nrbeval(srf,{linspace(0.0,1.0,20) linspace(0.0,1.0,20)});
%! h = surf(squeeze(p(1,:,:)),squeeze(p(2,:,:)),squeeze(p(3,:,:)));
%! set(h,'FaceColor','blue','EdgeColor','blue');
%! title('First derivatives over a test surface.');
%!
%! npts = 5;
%! tt = linspace(0.0,1.0,npts);
%! dsrf = nrbderiv(srf);
%! 
%! [p1, dp] = nrbdeval(srf, dsrf, {tt, tt});
%! 
%! up2 = vecnormalize(dp{1});
%! vp2 = vecnormalize(dp{2});
%! 
%! hold on;
%! plot3(p1(1,:),p1(2,:),p1(3,:),'ro');
%! h1 = quiver3(p1(1,:),p1(2,:),p1(3,:),up2(1,:),up2(2,:),up2(3,:));
%! h2 = quiver3(p1(1,:),p1(2,:),p1(3,:),vp2(1,:),vp2(2,:),vp2(3,:));
%! set(h1,'Color','black');
%! set(h2,'Color','black');
%! 
%! hold off;

%!test
%! knots{1} = [0 0 0 1 1 1];
%! knots{2} = [0 0 0 .5 1 1 1];
%! knots{3} = [0 0 0 0 1 1 1 1];
%! cx = [0 0.5 1]; nx = length(cx);
%! cy = [0 0.25 0.75 1]; ny = length(cy);
%! cz = [0 1/3 2/3 1]; nz = length(cz);
%! coefs(1,:,:,:) = repmat(reshape(cx,nx,1,1),[1 ny nz]);
%! coefs(2,:,:,:) = repmat(reshape(cy,1,ny,1),[nx 1 nz]);
%! coefs(3,:,:,:) = repmat(reshape(cz,1,1,nz),[nx ny 1]);
%! coefs(4,:,:,:) = 1;
%! nurbs = nrbmak(coefs, knots);
%! x = rand(5,1); y = rand(5,1); z = rand(5,1);
%! tt = [x y z]';
%! ders = nrbderiv(nurbs);
%! [points,jac] = nrbdeval(nurbs,ders,tt);
%! assert(points,tt,1e-10)
%! assert(jac{1}(1,:,:),ones(size(jac{1}(1,:,:))),1e-12)
%! assert(jac{2}(2,:,:),ones(size(jac{2}(2,:,:))),1e-12)
%! assert(jac{3}(3,:,:),ones(size(jac{3}(3,:,:))),1e-12)
%! 
%!test
%! knots{1} = [0 0 0 1 1 1];
%! knots{2} = [0 0 0 0 1 1 1 1];
%! knots{3} = [0 0 0 1 1 1];
%! cx = [0 0 1]; nx = length(cx);
%! cy = [0 0 0 1]; ny = length(cy);
%! cz = [0 0.5 1]; nz = length(cz);
%! coefs(1,:,:,:) = repmat(reshape(cx,nx,1,1),[1 ny nz]);
%! coefs(2,:,:,:) = repmat(reshape(cy,1,ny,1),[nx 1 nz]);
%! coefs(3,:,:,:) = repmat(reshape(cz,1,1,nz),[nx ny 1]);
%! coefs(4,:,:,:) = 1;
%! coefs = coefs([2 1 3 4],:,:,:);
%! nurbs = nrbmak(coefs, knots);
%! x = rand(5,1); y = rand(5,1); z = rand(5,1);
%! tt = [x y z]';
%! dnurbs = nrbderiv(nurbs);
%! [points, jac] = nrbdeval(nurbs,dnurbs,tt);
%! assert(points,[y.^3 x.^2 z]',1e-10);
%! assert(jac{2}(1,:,:),3*y'.^2,1e-12)
%! assert(jac{1}(2,:,:),2*x',1e-12)
%! assert(jac{3}(3,:,:),ones(size(z')),1e-12)
