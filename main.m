% main.m - Launcher for UHS_CCT MATLAB GUI
% -----------------------------------------
% This script launches the UHS_CCT.mlapp MATLAB GUI.
%
% Ensure that 'MRST' is available in the `external/` folder.
%
% Author: [Emrah SARI]
% License: GPL v3
% -----------------------------------------

clc;        % Clear the command window
clear;      % Clear workspace variables
close all;  % Close any open figure windows

% Add custom function directory to the top of MATLAB's search path
addpath('functions');
addpath('external\MRST\');

disp('Launching UHS_CCT MATLAB GUI...');
try
    % Launch the UHS_CCT App
    app = UHS_CCT;
catch ME
    % Display error if the app fails to launch
    warning('Failed to launch UHS_CCT.mlapp: %s', ME.message);
end
