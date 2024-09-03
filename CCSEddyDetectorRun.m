 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takeyoshi Nagai@UMassD 12/9/2010 -- Applied in the California Current System
% see Nagai et al. 2015. https://doi.org/10.1002/2015JC010889
% In this modified version, we apply this algorithm in the Peru-Chile EBUS
% in Rosales-Quintana et al
%
% Detection of eddies
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc


pathin =  'your/output/data';
basein='name_output';
sy=7; %year initial
ey=8; % year end
sm=1; % month initial
em=12; % month end
depi=-0; % for surface
depe=-50; % for surface
Isoi= 25.5;
Isoe = 26.5;
crow= -2e-11;
crar=[pi*10000^2 pi*200000^2]; 
plotflag=10;
OW_smooth = 10;
tnum_critic = 7;
% detection starts
CCSEddyDetector(pathin,basein,sy,ey,sm,em,depi,depe,Isoi,Isoe,crar,crow,OW_smooth,tnum_critic,plotflag)








