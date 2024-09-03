function CCSEddyDetector(pathin,basein,sy,ey,sm,em,depi,depe,Isoi,Isoe,crar,crow,OW_smooth,tnum_critic,plotflag)

em0=em;
% em=12;
% input file check
flagfile=0;
for iy=sy:1:ey
    if iy==ey
        em=em0;
    end
    for im=sm:1:em
        fn=sprintf([basein 'Y%02dM%02d.nc.1'],iy,im);
        fn=fullfile(pathin,fn);
	dh=dir(fn);
	if size(dh,1)==0
	    disp(['no such file: ' fn])
	    flagfile=-1;
	end
    end
end
if flagfile==-1
    error('there are missing files')
end


%%% here any data bc it is just taking grid data info
pathgridn = ['/path_outputs/'];

fng = 'data.nc';
ncz=netcdf.open(fullfile(pathgridn,fng),'nc_nowrite');
[zw,dzw,err] = getZ_DZTN(ncz,1,1);
[zr,dzr,err] = getZ_DZTN(ncz,0,1);
zu=rho2u_3d(zr);
zv=rho2v_3d(zr);
dzwu=rho2u_3d(dzw);
dzwv=rho2v_3d(dzw);
[maskd,maskr,lon,lat,dx,dy,h,f,distr]=CCSETHdistmask;

[xi,yi]=meshgrid(1:size(maskr,2),1:size(maskr,1));
[dimlens, dimnames]=ncinqdimlen(ncz,'temp');
recn=dimlens(end); % days
counts=dimlens; counts(end)=1;
countsu=counts;
countsv=counts;
countstemp=counts;
countsu(1)=countsu(1)-1;
countsv(2)=countsv(2)-1;
netcdf.close(ncz)

%-------------------------------------------------
% starting the main loop
%-------------------------------------------------
% em0=em;
% em=12;

for iy=sy:1:ey
  j=0;  
    for im=sm:1:em
        
        irec=0;
        for irec= 1:1:recn

            j=j+1;
            fn=sprintf([basein 'Y%02dM%02d.nc.1'],iy,im);
            fn=fullfile(pathin,fn);
            ncin=netcdf.open(fn,'nc_nowrite');
            disp(['Processing Year:' num2str(iy) ' Month:' num2str(im) ' Day:' num2str(irec)] )


            % for opening ncdf
            u=getvarnetcdfTN(ncin,'u',[0 0 0 irec-1],[countsu]);
            v=getvarnetcdfTN(ncin,'v',[0 0 0 irec-1],[countsv]);
            temp=getvarnetcdfTN(ncin,'temp',[0 0 0 irec-1],[countstemp]);
            sal=getvarnetcdfTN(ncin,'salt',[0 0 0 irec-1],[countstemp]);
            sigma3d = sigmat(temp,sal);


            % getting velocities with the same size as for sigma
            ur=0.5.*(u(:,:,1:end-1)+u(:,:,2:end));
            ur=cat(3,u(:,:,1),ur,u(:,:,end));
            vr=0.5.*(v(:,1:end-1,:)+v(:,2:end,:));
            vr=cat(2,v(:,1,:),vr,v(:,end,:));

            % Create a mask matrix with NaNs in the boundary regions
            grid_cells = 12;
            mask = ones(size(ur));
            mask(:,1:grid_cells,:) = NaN; % Mask the top rows
            mask(:,end-grid_cells+1:end, :) = NaN; % Mask the bottom rows
            mask(:,:, 1:grid_cells) = NaN; % Mask the left columns
            mask(:,:, end-grid_cells+1:end) = NaN; % Mask the right columns

            %%%%%%%% masked
            urx = ur.*mask;
            vrx = vr.*mask;
            sigma3dx = sigma3d.*mask;

            % for surface eddies
            indx1 = double(zr<depi & zr>depe);
            usurf = sum(urx.*dzr.*indx1,1)./sum(dzr.*indx1,1);
            usurf = squeeze(usurf);
            vsurf = sum(vrx.*dzr.*indx1,1)./sum(dzr.*indx1,1);
            vsurf = squeeze(vsurf);

            % for puddies averaged velocitites within isopycnals 25.5 to 26.5
            indx2 = double(sigma3dx>=Isoi & sigma3dx<=Isoe);
            uIso = sum(urx.*dzr.*indx2,1)./sum(dzr.*indx2,1);
            uIso = squeeze(uIso);
            vIso = sum(vrx.*dzr.*indx2,1)./sum(dzr.*indx2,1);
            vIso = squeeze(vIso);


            % detecting eddies
            [data(j).eddy] = DetectEddyOW(usurf,vsurf,usub,vsub,lon,lat,xi,yi,dx,dy,f,maskr,urx,vrx,sigma3dx,crar,crow,OW_smooth,tnum_critic,plotflag);

        end

    end
    outfn=sprintf(['Eddy_' basein 'EddyDetect.mat'],iy);
    netcdf.close(ncin)
    disp([outfn ' : now saving...'])
    save(outfn,'data','xi','yi','lon','lat','dx','dy')
    disp([outfn ' is saved'])

end






