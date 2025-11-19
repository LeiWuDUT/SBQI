function interpolatedToolOrientation = generateToolOrientationBySquad(paramArray, keyframes, indKeyframes)
    % Generate tool orientation using quaternion squad interpolation
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

    % Convert tool orientation vectors to quaternion objects
    numOfKeyframes = length(indKeyframes); % Number of keyframes
    keyframeQuaternions = repmat(quaternion([1, 0, 0, 0]), numOfKeyframes, 1);
    for idx = 1:numOfKeyframes
        keyframeQuaternions(idx) = convertToolOrientationToQuaternion(keyframes(idx, :));
    end

    % Get control points for squad interpolation
    controlPointsSi = getControlPointSi(keyframeQuaternions);

    % Preallocate output array
    numCCPoints = length(paramArray);
    interpolatedToolOrientation = zeros(numCCPoints, 3);
    referenceVector = [0, 0, 1];
    
    % Perform squad interpolation between keyframes
    for idx = 1:numOfKeyframes - 1
        % Current and next keyframe indices and their corresponding quaternions
        fIter = indKeyframes(idx);
        sIter = indKeyframes(idx + 1);
        for i = fIter:sIter - 1
            tParameter = (paramArray(i) - paramArray(fIter)) / (paramArray(sIter) - paramArray(fIter));
            interQuaternion = quaternionSquad(keyframeQuaternions(idx), keyframeQuaternions(idx + 1), ...
                                              controlPointsSi(idx), controlPointsSi(idx + 1), tParameter);
            interpolatedToolOrientation(i, :) = rotatepoint(interQuaternion, referenceVector);
        end
    end

    % Store the last keyframe's tool orientation
    interpolatedToolOrientation(indKeyframes(end), :) = keyframes(end, :);
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

% Helper Function: Get Control Points for Squad Interpolation
function Si = getControlPointSi(quaternionArray)
    % Calculate control points for squad interpolation, given a quaternion array
    
    % Number of quaternions
    numOfQuat = length(quaternionArray);
    
    % Preallocate control points array
    Si = repmat(quaternion([1, 0, 0, 0]), numOfQuat, 1);
    Si(1) = quaternionArray(1);
    Si(end) = quaternionArray(end); % The first and last points are the same as the interpolation points
    
    % Calculate control points for intermediate quaternions
    for ii = 2:numOfQuat-1
        qi = quaternionArray(ii);
        qi_m1 = quaternionArray(ii-1);
        qi_a1 = quaternionArray(ii+1);
        
%         % Ensure the shortest path for quaternion multiplication
%         if dot(qi, qi_m1) < 0
%             qi_m1 = -qi_m1;
%         end
%         if dot(qi, qi_a1) < 0
%             qi_a1 = -qi_a1;
%         end
        
        % Calculate the control point using quaternion logarithm and exponentiation
        qi_conj = conj(qi);
        m0 = qi_conj * qi_m1;
        m1 = qi_conj * qi_a1;
        
        m0_log = log(m0);
        m1_log = log(m1);
        
        m_log_sum = m0_log + m1_log;
        k = -0.25 * m_log_sum;
        k_exp = exp(k);
        
        Si(ii) = qi * k_exp;
    end
end

% Helper Function: Quaternion Squad Interpolation
function quatSquad = quaternionSquad(qi, qi_1, si, si_1, tparam)
    % Perform squad interpolation on multiple keyframes

    % Perform spherical linear interpolation (slerp)
    k1 = slerp(qi, qi_1, tparam);
    k2 = slerp(si, si_1, tparam);

    % Perform squad interpolation
    quatSquad = slerp(k1, k2, 2 * tparam * (1 - tparam));
end

