function uParameters = centripetalParameterization(pointsMatrix)
    % Calculate the centripetal parameterization of a set of points.
    % The function takes as input a matrix of points and returns a row vector of parameter values.
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

    numOfPoints = size(pointsMatrix, 1);

    % Preallocate array for square root of distances
    sqrtDistBetweenPoints = zeros(numOfPoints-1, 1);
    
    % Calculate the square root of the distance between each pair of neighboring points
    for ii = 1 : numOfPoints - 1
        sqrtDistBetweenPoints(ii) = sqrt(norm(pointsMatrix(ii, :) - pointsMatrix(ii + 1, :)));
    end
    
    totalDistance = sum(sqrtDistBetweenPoints);

    % Preallocate uParameter and set the first value to 0
    uParameters = zeros(1, numOfPoints);
    
    % Calculate the cumulative distance of each point from the first point as a proportion of the total distance
    uParameters(2:end) = cumsum(sqrtDistBetweenPoints) / totalDistance;
    
    % Ensure the last parameter value is 1
    uParameters(end) = 1; 
end
