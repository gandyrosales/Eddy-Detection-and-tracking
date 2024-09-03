function [maskd,maskr,long,latg,dx,dy,h,f,distr]=CCSETHdistmask

pathgridn = ['/path_output/'];
fn = 'name_output.nc';

ncgrd=netcdf.open([pathgridn fn],'nc_nowrite');
long=getvarnetcdfTN(ncgrd,'lon_rho');
latg=getvarnetcdfTN(ncgrd,'lat_rho');
dx=1./getvarnetcdfTN(ncgrd,'pm'); % in meter
dy=1./getvarnetcdfTN(ncgrd,'pn'); % in meter
maskr=getvarnetcdfTN(ncgrd,'mask_rho');
f=getvarnetcdfTN(ncgrd,'f');
h=getvarnetcdfTN(ncgrd,'h');
netcdf.close(ncgrd)

maskd = maskr;

distc1 = fliplr(cumsum(fliplr(dx.*maskd),2));
distc2 = cumsum(dx.*(maskd-1),2);
distr = distc1+distc2;
