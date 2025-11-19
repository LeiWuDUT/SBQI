function plotToolOrientationInterpolationWithSlider(CCPoints, ...
                                                     fv, ...
                                                     interToolOrientationSlerp, ...
                                                     interToolOrientationSquad, ...
                                                     interToolOrientationSBQI, ...
                                                     keyframes, ...
                                                     indKeyframes, ...
                                                     lenToolOrien)
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
% Plot tool orientation interpolation with two panels, checkboxes and slider
% Inputs:
%   CCPoints: Center curve points (Nx3)
%   fv: Infos of the surface, including faces & vertices, for visualization
%   interToolOrientationSlerp: Interpolated orientations (Slerp) (Nx3)
%   interToolOrientationSquad: Interpolated orientations (Squad) (Nx3)
%   interToolOrientationSBQI: Interpolated orientations (SBQI) (Nx3)
%   keyframes: Keyframe orientations (Mx3)
%   indKeyframes: Indices of keyframes in CCPoints
%   lenToolOrien: Arrow length scaling factor

%% Create figure
hFig = figure('Name', 'Tool Orientation Interpolation', ...
              'Color', 'w', 'Position', [100, 100, 1200, 600]);

%% Create panels
hPanel1 = uipanel('Parent', hFig, 'Title', 'Orientation on Sphere', ...
                  'FontSize', 12, 'Position', [0.02 0.2 0.48 0.78]);
hPanel2 = uipanel('Parent', hFig, 'Title', 'Spatial Path', ...
                  'FontSize', 12, 'Position', [0.5 0.2 0.48 0.78]);

%% Create axes
hAx1 = axes('Parent', hPanel1);
hold(hAx1, 'on');
axis(hAx1, 'equal');
grid(hAx1, 'on');
view(hAx1, [130, 32]);
xlabel(hAx1, 'i');
ylabel(hAx1, 'j');
zlabel(hAx1, 'k');
title(hAx1, 'Orientation on Unit Sphere');
set(hAx1, 'FontSize', 12, 'FontName', 'Times New Roman');

hAx2 = axes('Parent', hPanel2);
hold(hAx2, 'on');
axis(hAx2, 'equal');
grid(hAx2, 'on');
view(hAx2, [130, 32]);
xlabel(hAx2, 'X');
ylabel(hAx2, 'Y');
zlabel(hAx2, 'Z');
title(hAx2, 'Spatial Path');
set(hAx2, 'FontSize', 12, 'FontName', 'Times New Roman');

%% Draw unit sphere on hAx1
drawUnitSphere(hAx1);

%% Plot keyframes, surface, and CC curve
hKeyframes = quiver3(hAx2, CCPoints(indKeyframes,1), ...
                           CCPoints(indKeyframes,2), ...
                           CCPoints(indKeyframes,3), ...
                           (lenToolOrien+5)*keyframes(:,1), ...
                           (lenToolOrien+5)*keyframes(:,2), ...
                           (lenToolOrien+5)*keyframes(:,3), ...
                           'k-', 'AutoScale', 'off', 'LineWidth', 1.5);
                       
hSurf = patch(hAx2, fv,  'FaceColor',       [120 250 255]/255, ...
                   'EdgeColor',       'none',        ...
                   'FaceLighting',    'gouraud',     ...
                   'AmbientStrength', 0.8);
camlight(hAx2);           % Add a light source
material(hSurf, 'shiny'); % Set the material property
view(hAx2, [130,32]);     % view angle

plot3(hAx2, CCPoints(:, 1), CCPoints(:, 2), CCPoints(:, 3), 'r', 'LineWidth', 2)

%% Initialize plots for Slerp, Squad, SBQI
hPlots = struct();

[hPlots.Slerp.arrow, hPlots.Slerp.line] = plotInterpolation(hAx1, hAx2, CCPoints, interToolOrientationSlerp, lenToolOrien, 'b');
[hPlots.Squad.arrow, hPlots.Squad.line] = plotInterpolation(hAx1, hAx2, CCPoints, interToolOrientationSquad, lenToolOrien, 'r');
[hPlots.SBQI.arrow,  hPlots.SBQI.line]  = plotInterpolation(hAx1, hAx2, CCPoints, interToolOrientationSBQI, lenToolOrien, 'g');
legend([hKeyframes hPlots.Slerp.arrow hPlots.Squad.arrow hPlots.SBQI.arrow],{'Keyframes', 'Slerp','Squad', 'SBQI'})
legend([hPlots.Slerp.line hPlots.Squad.line hPlots.SBQI.line],{'Slerp','Squad', 'SBQI'})


%% Create control panel (checkboxes and slider)
hControlPanel = uipanel('Parent', hFig, 'Title', 'Controls', ...
                        'FontSize', 12, 'Position', [0.02 0.02 0.96 0.16]);

% Checkbox for Slerp
hCheckboxSlerp = uicontrol('Parent', hControlPanel, 'Style', 'checkbox', ...
    'String', 'Slerp', 'Units', 'normalized', 'Position', [0.05 0.6 0.1 0.3], ...
    'Value', 1, 'FontSize', 12, ...
    'Callback', @(src, ~) toggleVisibility(src, hPlots.Slerp));

% Checkbox for Squad
hCheckboxSquad = uicontrol('Parent', hControlPanel, 'Style', 'checkbox', ...
    'String', 'Squad', 'Units', 'normalized', 'Position', [0.15 0.6 0.1 0.3], ...
    'Value', 1, 'FontSize', 12, ...
    'Callback', @(src, ~) toggleVisibility(src, hPlots.Squad));

% Checkbox for SBQI
hCheckboxSBQI = uicontrol('Parent', hControlPanel, 'Style', 'checkbox', ...
    'String', 'SBQI', 'Units', 'normalized', 'Position', [0.25 0.6 0.1 0.3], ...
    'Value', 1, 'FontSize', 12, ...
    'Callback', @(src, ~) toggleVisibility(src, hPlots.SBQI));

% Slider to control number of points
numPoints = size(interToolOrientationSlerp, 1);
hSlider = uicontrol('Parent', hControlPanel, 'Style', 'slider', ...
    'Units', 'normalized', 'Position', [0.4 0.4 0.5 0.4], ...
    'Min', 2, 'Max', numPoints, 'Value', numPoints, ...
    'SliderStep', [1/(numPoints-2) , 5/(numPoints-2)], ...
    'Callback', @(src, ~) updatePlots(round(src.Value), hPlots, CCPoints, ...
                                      interToolOrientationSlerp, ...
                                      interToolOrientationSquad, ...
                                      interToolOrientationSBQI, ...
                                      lenToolOrien));

% Slider label
uicontrol('Parent', hControlPanel, 'Style', 'text', ...
    'Units', 'normalized', 'Position', [0.4 0.05 0.5 0.3], ...
    'String', 'Select number of points', ...
    'FontSize', 12);

end

%% ===== Helper Functions =====

function drawUnitSphere(hAx)
% Draw a unit sphere on the given axes
    [uSphere, vSphere, wSphere] = sphere(56);
    mesh(hAx, uSphere, vSphere, wSphere);
%     mesh(hAx, uSphere, vSphere, wSphere, ...
%          'FaceAlpha', 0.3, 'EdgeAlpha', 0.5, ...
%          'EdgeColor', [0.7 0.7 0.7]);
end

function [hArrow, hLine] = plotInterpolation(hAx1, hAx2, CCPoints, orientations, lenToolOrien, color)
% Plot interpolation arrows and spatial path

    hArrow = quiver3(hAx2, CCPoints(:,1), CCPoints(:,2), CCPoints(:,3), ...
                           lenToolOrien*orientations(:,1), ...
                           lenToolOrien*orientations(:,2), ...
                           lenToolOrien*orientations(:,3), ...
                           '-', 'Color', color, ...
                           'AutoScale', 'off', 'LineWidth', 0.5);
                       
    hLine  = plot3(hAx1, orientations(:, 1), ...  
                         orientations(:, 2), ...  
                         orientations(:, 3), ...  
                         'Color',color, 'LineWidth', 2, 'LineStyle', '-');  
end

function toggleVisibility(src, hPlot)
% Toggle visibility of interpolation plots
    if src.Value
        set(hPlot.arrow, 'Visible', 'on');
        set(hPlot.line, 'Visible', 'on');
    else
        set(hPlot.arrow, 'Visible', 'off');
        set(hPlot.line, 'Visible', 'off');
    end
end

function updatePlots(nPoints, hPlots, CCPoints, slerp, squad, sbqi, lenToolOrien)
% Update plots based on selected number of interpolation points

    % Update Slerp
    updateSinglePlot(hPlots.Slerp, CCPoints(1:nPoints,:), slerp(1:nPoints,:), lenToolOrien);

    % Update Squad
    updateSinglePlot(hPlots.Squad, CCPoints(1:nPoints,:), squad(1:nPoints,:), lenToolOrien);

    % Update SBQI
    updateSinglePlot(hPlots.SBQI, CCPoints(1:nPoints,:), sbqi(1:nPoints,:), lenToolOrien);
end

function updateSinglePlot(hPlot, CCPointsSubset, orientationsSubset, lenToolOrien)
% Update single interpolation plot (arrow + line)
    set(hPlot.arrow, 'XData', CCPointsSubset(:,1), ...
                     'YData', CCPointsSubset(:,2), ...
                     'ZData', CCPointsSubset(:,3), ...
                     'UData', lenToolOrien * orientationsSubset(:,1), ...
                     'VData', lenToolOrien * orientationsSubset(:,2), ...
                     'WData', lenToolOrien * orientationsSubset(:,3));
    
    set(hPlot.line, 'XData', orientationsSubset(:,1), ...
                    'YData', orientationsSubset(:,2), ...
                    'ZData', orientationsSubset(:,3));
end
