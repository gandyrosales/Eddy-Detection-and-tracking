 % For subsurfces eddies using rooms outputs   
clear all
close all
clc


pathin =  'your/output/data/folder';
basein='name_output';
% use any year needed, e.g. Y7 to Y8
sy=7; %year initial
ey=8; % year end
sm=1; % month initial
em=12; % month end

% layers criteria
depi=-0; % for surface
depe=-50; % for surface
Isoi= 25.5; %  for subsurface
Isoe = 26.5; % for subsurface
crow= -2e-11; % criteria for OW
crar=[pi*10000^2 pi*200000^2]; 
plotflag=10;
OW_smooth = 10;
tnum_critic = 7;

CCSEddyDetector(pathin,basein,sy,ey,sm,em,depi,depe,Isoi,Isoe,crar,crow,OW_smooth,tnum_critic,plotflag)








