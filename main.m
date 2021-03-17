% MAIN analysis script
% This scripts runs the analysis and plots/saves figures.
clear('all');

% Change to the repo path and add the src directory
PREV_DIRECTORY = pwd;
cd(fileparts(which(matlab.desktop.editor.getActiveFilename)));
addpath('src');

% Run the scripts
analysis_formatData;
analysis_classification;
analysis_statistics;

figure_classifier;
figure_obInjection;
figure_tiling;
figure_pcInjection;
figure_mapping;

% Return to previous directory
cd(PREV_DIRECTORY);