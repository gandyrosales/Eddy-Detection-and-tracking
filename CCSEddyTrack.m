function CCSEddyTrack(pathin,basein,sy,ey,sm,em)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takeyoshi Nagai@UMassD 12/9/2010 -- Applied in the California Current System
% see Nagai et al. 2015. https://doi.org/10.1002/2015JC010889
% In this modified version, we apply this algorithm in the Peru-Chile EBUS
% in Rosales-Quintana et al
%
% Tracking function of detected eddies
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------
% starting the main loop
%-------------------------------------------------
for iy=sy:1:ey
        fn=sprintf([basein 'detected_output_name.mat'],iy);
        fn=fullfile(pathin,fn);
    	load(fn);


    for im=sm:1:em
         disp(['Processing Year:' num2str(iy) 'month: ' num2str(im)]);
       
    	recn=size(data,2);
    	if iy==sy
    	    peddy=data(1).eddy;
    	    recstr=2;
    	else
    	    recstr=1;
    	end
    	for irec=recstr:1:recn
    	    eddy=data(irec).eddy;
    	    eddy=CCSEddyTracker(peddy,eddy,dx,dy);
    	    peddy=eddy;
    	    data(irec).eddy=eddy;
    	end
    	outfn=fn;
        disp(' now saving...')
    	save(outfn,'data','-append')
    end
end


