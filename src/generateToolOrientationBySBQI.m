function interpolatedToolOrientation = generateToolOrientationBySBQI(uParamsArray, keyframes, indKeyframes)
    % Generate tool orientation using quaternion SBQI interpolation
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
    
    % Convert quaternion objects to rotation vectors
    paramsKeyframes = uParamsArray(indKeyframes);
    rotationVectorsKeyframes = zeros(numOfKeyframes, 3);
    for ii = 1 : numOfKeyframes
        % [quatAngle, quatAxis] = extractQuatAngleAxis(keyframeQuaternions(ii));
        % rotationVectorsKeyframes(ii, :) = quatAngle*quatAxis;
        rotationVectorsKeyframes(ii, :) = rotvec(keyframeQuaternions(ii));
        % built-in function with the same functionality
    end
    
    % Interpolate the rotation vectors using quintic B-spline
    numOfCCPoints = length(uParamsArray);
    degree = 5;
    % user-defined minmum and maximum values for the number of the control points
    nMinCtrlpts = floor(numOfCCPoints/2);
    nMaxCtrlpts = floor(numOfCCPoints*4/5);
    stepIter = 2;
    numCtrlpts = max([numOfKeyframes, degree+1, nMinCtrlpts]);
    iter = 1;
    maxIter = 100;
    crvEnergyFormer = -1;
    fTol = 1e-5;

    while iter < maxIter
        kntVec = assignKnotVectorUsingAVGKTP(numCtrlpts, uParamsArray, degree);
        [ctrlpts, crvEnergy] = bspinterpcrvenergy(rotationVectorsKeyframes', paramsKeyframes, degree, kntVec);
        % check termination
        if crvEnergyFormer ~= -1
            exitFlag = (abs(crvEnergyFormer - crvEnergy) < fTol*max(1, crvEnergyFormer)) || ...
                       (numCtrlpts >= nMaxCtrlpts);
            if exitFlag
                break;
            end
        end
        crvEnergyFormer = crvEnergy;
        iter = iter + 1;
        numCtrlpts = numCtrlpts + stepIter; 
    end
    bscrv = nrbmak([ctrlpts; ones(1,size(ctrlpts, 2))], kntVec);
    
    % Calculate intermediate rotation vectors and tool orientations
    interpolatedToolOrientation  = zeros(numOfCCPoints, 3);
    pntsOnCrv = nrbeval(bscrv, uParamsArray); % intermediate rotation vectors
    referenceVector = [0 0 1];
    
    for ii = 1 : numOfCCPoints
        % angleInterpolatedQuat = norm(pntsOnCrv(:, ii));
        % axisInterpolatedQuat  = pntsOnCrv(:, ii) / angleInterpolatedQuat;
        % interpolatedQuat_ = [cos(angleInterpolatedQuat/2), ...
        %                      sin(angleInterpolatedQuat/2)*axisInterpolatedQuat'];
        % interpolatedQuat  = quaternion(interpolatedQuat_);
        interpolatedQuat = quaternion(pntsOnCrv(:, ii)', 'rotvec');
        % built-in function with the same functionality
        
        interpolatedToolOrientation(ii, :) = rotatepoint(interpolatedQuat, referenceVector);
    end
end

% Helper function: Extract the rotation angle and axis of the quaternion
function [quatAngle, quatAxis] = extractQuatAngleAxis(q)
    % Convert quaternion into four-element vector
    q_compact = compact(q);
    
    % Extract the imaginary components
    x = q_compact(2);
    y = q_compact(3);
    z = q_compact(4);

    % Compute the norm of the imaginary components
    norm_imag = norm([x, y, z]);

    % Compute the angle using atan2
    angle_half = atan(norm_imag/q_compact(1));

    % Complete the angle by multiplying by 2
    quatAngle = 2 * angle_half;
    
    % Extract the axis
    quatAxis = [x, y, z] / norm_imag;
end

% Helper function: Assign knots using the technique of AVGKTP.
function knts = assignKnotVectorUsingAVGKTP(nCtrlpts, uParamArray, degree)
    % CalKnotVectorUsingAVGKTP: Assign knots using the technique of AVGKTP.
    % refer to the Nurbs Book
    nParamArray = length(uParamArray);
    knts        = zeros(1, nCtrlpts + degree + 1);
    const       = nParamArray / (nCtrlpts-degree);
    for i = 1:nCtrlpts-degree-1
        coefs  = floor(const*i);
        alfa = const*i-coefs;
        knts(i+1+degree) = (1-alfa)*uParamArray(coefs) + alfa*uParamArray(coefs+1);
    end
    for i = 1:degree+1
       knts(i+nCtrlpts) = 1;
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
    

