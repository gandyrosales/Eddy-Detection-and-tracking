%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takeyoshi Nagai@UMassD 12/9/2010 -- Applied in the California Current System
% see Nagai et al. 2015. https://doi.org/10.1002/2015JC010889
% In this modified version, we apply this algorithm in the Peru-Chile EBUS
% in Rosales-Quintana et al
%
% Tracking function of detected eddies
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% Set paths and outputs name
pathin = '/Volumes/18TB_GR/Elements2/ROMS_eddyDetection/Sub_EddyTrackerOriginal/Detected_Eddy_plus_Puddy';
basein = 'Output_name';

% years and months to be use
sy=7; 
ey=8;
sm=1;
em=12;
% tracking starts
CCSEddyTrack(pathin,basein,sy,ey,sm,em,plotflag)
