function [Vort,OW,S,Rabs]=CPRNCScompCurl(us,vs,dx,dy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takeyoshi Nagai@UMassD 12/9/2010 -- Applied in the California Current System
% see Nagai et al. 2015. https://doi.org/10.1002/2015JC010889
% In this modified version, we apply this algorithm in the Peru-Chile EBUS
% in Rosales-Quintana et al -- 2024, september.
%
% To compute relative vorticity
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ux,uy]=gradient(ur);
[vx,vy]=gradient(vr);
uy=uy./dy;
vx=vx./dx;
ux=ux./dx;
vy=vy./dy;
Vort=vx-uy;

S2= (0.25.*(vx+uy).^2-ux.*vy);
S2(S2<0)=NaN;
S=sqrt(S2);
OW=4.*ux.^2+4.*vx.*uy;
% for absolute valocity
Rabs= sqrt(us.^2 + vs.^2);