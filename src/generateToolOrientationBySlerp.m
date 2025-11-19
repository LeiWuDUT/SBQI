function interpolatedToolOrientation = generateToolOrientationBySlerp(uParamsArray, Keyframes, indKeyframes)
    % Generate tool orientation using quaternion spherical linear interpolation (Slerp)
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
    
    % Preallocate quaternion array
    numOfKeyframes = length(indKeyframes);
    keyframeQuaternions = repmat(quaternion([1, 0, 0, 0]), numOfKeyframes, 1);
    referenceVector = [0, 0, 1];
    
    % Convert tool orientation vectors to quaternion objects
    for idx = 1:numOfKeyframes
        keyframeQuaternions(idx) = convertToolOrientationToQuaternion(Keyframes(idx, :));
    end
    
    % Preallocate output matrix
    numCCPoints = length(uParamsArray);
    interpolatedToolOrientation = zeros(numCCPoints, 3);
    
    % Calculate tool orientation by Slerp interpolation
    for idx = 1:numOfKeyframes - 1
        fQuaternion = keyframeQuaternions(idx);
        sQuaternion = keyframeQuaternions(idx + 1);
        startIndex = indKeyframes(idx);
        endIndex = indKeyframes(idx + 1);
        for jdx = startIndex:endIndex
            % We use the parameterization of the CC points for interpolation 
            % to better capture the tool orientation changes.
            uParam = (uParamsArray(jdx) - uParamsArray(startIndex)) / (uParamsArray(endIndex) - uParamsArray(startIndex));
            quaternionInterpolated = slerp(fQuaternion, sQuaternion, uParam);
            interpolatedToolOrientation(jdx, :) = rotatepoint(quaternionInterpolated, referenceVector);
        end
    end
end


% Helper Function: Convert Tool Orientation to Quaternion
function quat = convertToolOrientationToQuaternion(toolOrientation, referenceToolOrientation)
    % Convert tool orientation vector to quaternion representation
    
    % Default reference tool orientation
    if nargin < 2
        referenceToolOrientation = [0, 0, 1];
    end
    referenceToolOrientation = referenceToolOrientation / norm(referenceToolOrientation);
    
    % Handle the case where tool orientation is the same as reference
    if isequal(referenceToolOrientation, toolOrientation)
        quat = quaternion(1, 0, 0, 0);
        return;
    elseif isequal(referenceToolOrientation, -toolOrientation)
        quat = quaternion(0, 1, 0, 0);
        return;
    end
    
    % Calculate the rotation axis and angle
    rotationAxis = cross(referenceToolOrientation, toolOrientation);
    rotationAxis = rotationAxis / norm(rotationAxis);
    rotationAngle = acos(dot(toolOrientation, referenceToolOrientation));
    
    % Create the quaternion
    quat = quaternion([cos(rotationAngle / 2), sin(rotationAngle / 2) * rotationAxis]);
end
