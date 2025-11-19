function hessianMatrix = calHessianMat(numCtrlPoints, degree, knotVector, paramBegin, paramEnd, order)
    % Compute the Hessian Matrix when interpolating using curve energy method 
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

    numDiscretized = 1000;                                  % Divide the parameter domain [0, 1] into 1000 parts
    uParamArray = linspace(paramBegin, paramEnd, numDiscretized);  % Discrete parameter sequence
    uParamInterval = 1 / (numDiscretized - 1);
    indexArray = zeros(1, numDiscretized);                  % B-spline basis function indices for each discrete parameter
    hessianMatrix = zeros(numCtrlPoints, numCtrlPoints);    % Symmetric matrix with a fixed bandwidth of 2*(degree+1)
    tempDersMtx = zeros(numDiscretized, numCtrlPoints);     % Each row is the derivative value at a fixed parameter, each column is the derivative value at a fixed index
    
    % Calculate the derivatives of b-spline basis
    for idxI = 1 : length(uParamArray)
        indexArray(idxI) = findspan(numCtrlPoints - 1, degree, uParamArray(idxI), knotVector);
        tempDers = basisfunder(indexArray(idxI), degree, uParamArray(idxI), knotVector, order);
        tempDers = squeeze(tempDers(1, :, :));
        tempDersMtx(idxI, indexArray(idxI) - degree + 1 : indexArray(idxI) + 1) = tempDers(order + 1, :);
    end
    
    % Assemble the Hessian matrix using the derivatives
    for idxI = 0 : degree - 1
        for idxJ = 0 : idxI
            hessianMatrix(idxI + 1, idxJ + 1) = sum(tempDersMtx(:, idxI + 1) .* tempDersMtx(:, idxJ + 1)) * uParamInterval;
        end
    end

    for idxI = degree : numCtrlPoints - 1
        for idxJ = idxI - degree : idxI
            hessianMatrix(idxI + 1, idxJ + 1) = sum(tempDersMtx(:, idxI + 1) .* tempDersMtx(:, idxJ + 1)) * uParamInterval;
        end
    end

    for idxI = 1 : numCtrlPoints
        for idxJ = 1 : numCtrlPoints
            hessianMatrix(idxI, idxJ) = hessianMatrix(idxJ, idxI);
        end
    end
end