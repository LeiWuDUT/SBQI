% Main function for quaternion interpolation methods: Slerp, Squad, and SBQI
%
%    Copyright (C) 2025 Lei Wu, Dalian University of Technology
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program. If not, see <http://www.gnu.org/licenses/>.

%% Initialization
clear; clc; close all;
addpath(genpath(pwd));

%% Load workpiece and keyframes
% Choose one of the following cases:
% load('caseTwoKeyframes.mat'); % Case 0: TwoKeyframes
load('caseImpeller.mat');     % Case 1: Impeller    
% load('caseSurfObs.mat');      % Case 2: Surface with obstacles
% load('caseBlade.mat');        % Case 3: Blade


% Set visualization parameters
lenToolOrien = 50; % Arrow length for tool orientations

%% Parameterization of CC curve
uParamsArray = centripetalParameterization(CCPoints); % centripetal parameterization

%% Quaternion interpolation
interToolOrientationSlerp = generateToolOrientationBySlerp(uParamsArray, keyframes, indKeyframes);
interToolOrientationSquad = generateToolOrientationBySquad(uParamsArray, keyframes, indKeyframes);
interToolOrientationSBQI  = generateToolOrientationBySBQI(uParamsArray, keyframes, indKeyframes);

%% Plot results using modularized function
if 1
    plotToolOrientationInterpolationWithSlider( ...
        CCPoints, ...
        fv, ...
        interToolOrientationSlerp, ...
        interToolOrientationSquad, ...
        interToolOrientationSBQI, ...
        keyframes, ...
        indKeyframes, ...
        lenToolOrien ...
    );
end

%% End of main function