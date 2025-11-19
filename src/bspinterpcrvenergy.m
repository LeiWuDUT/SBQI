function [ctrlpts, crvenergy] = bspinterpcrvenergy(Q, u, p, knts)
    %
    % bspinterpcrvenergy: B-Spline interpolation with a given knot vector by minimizing curve energy.
    %
    % Calling Sequence:
    %
    %   [ctrlpts, crvenergy] = bspinterpcrvenergy(Q, u, p, knts)
    %
    %    INPUT:
    %
    %      Q      - points to be interpolated in the form [x_coord; y_coord; z_coord].
    %      u      - parameters associated with the points Q.
    %      p      - degree of the interpolation curve.
    %      knts   - knot vector of the interpolation curve.
    %
    %    OUTPUT:
    %
    %      ctrlpts - control points of the B-Spline curve.
    %      crvenergy - the curve energy of the B-spline curve.
    %
    %    Copyright (C) 2025 Lei Wu, Dalian University of Technology
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
    %%


    n = size (Q, 2); % the number of the interpolation points Q
    ndim = size(Q, 1); % dimension of the interpolation points Q
    nc = length(knts)-p-1; % the number of control points

    % interpolation process
    H = 2*calHessianMat(nc, p, knts, 0, 1, 2); % hessian matrix for smoothness purpose, times 2 or not does not matter (only influences the Lagrange multiplier)
    A = zeros (n, nc); % matrix that contains all basis function evaluated at u
    A(1,1) = 1;
    A(n,nc) = 1;
    for ii=2:n-1
    span = findspan (nc-1, p, u(ii), knts);
    A(ii,span-p+1:span+1) = basisfun (span, u(ii), p, knts);
    end
    coefmtx = [H, A'; A, zeros(n, n)];
    rhsvec  = [zeros(nc, ndim); Q'];

    pnts = coefmtx \ rhsvec;
    ctrlpts = pnts(1:nc, :)';

    % curve energy
    crvenergy = 0;
    for ii = 1:size(ctrlpts, 1) % dimension of ctrlpts
        crvenergy = crvenergy + ctrlpts(ii, :)*H*ctrlpts(ii, :)';
    end
    crvenergy = 0.5*crvenergy; % using hessian matrix that times 2
end