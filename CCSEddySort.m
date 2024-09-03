function Eddies=CCSEddySort(data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takeyoshi Nagai@UMassD 12/9/2010 -- Applied in the California Current System
% see Nagai et al. 2015. https://doi.org/10.1002/2015JC010889
% In this modified version, we apply this algorithm in the Peru-Chile EBUS
% in Rosales-Quintana et al -- 2024, september.
%
% Tracking function of detected eddies
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort eddy dataset by ednum stored in
% data.eddy.ednum
% resultant dataset struct Eddies will have
%  Eddies(ednum).data(i).x
%  Eddies(ednum).data(i).y
%  Eddies(ednum).data(i).ind2d
%  Eddies(ednum).data(i).zeta
%  Eddies(ednum).data(i).time
%
% so that one can follow an eddy history easier.
% It is recommanded that the data is cat-ed to have
% entire record.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nums=size(data,2);
% because the first index of data, eddy
% does not have ednum, it should be added.
data1=data(1);
ned1=size(data1.eddy,2);
for i1=1:1:ned1
    data(1).eddy(i1).ednum=i1;
end
% The last entry of data is taken and max 
% number of eddy is read.
data2=data(end);
ned2=size(data2.eddy,2);
maxed=-999;
for i2=1:1:ned2
    maxed=max([maxed data2.eddy(i2).ednum]);
end
% prepare the time index for all the eddy as a vector, I.
I = zeros(maxed,1);

% Start main loop for each snapshot
for i=1:1:nums
    eddy=data(i).eddy;
    ned=size(eddy,2);
    time=data(i).time;
    % extracting eddy data from each eddy
    for ied=1:1:ned
        ednum=eddy(ied).ednum;
        x=eddy(ied).x;
        y=eddy(ied).y;
        cx=eddy(ied).cx;
        cy=eddy(ied).cy;
        clon=eddy(ied).clon;
        clat=eddy(ied).clat;
        lon=eddy(ied).lon;
        lat=eddy(ied).lat;
        Curlz_sub=eddy(ied).Curlz_sub;
        radious=eddy(ied).radious;
        Rvabs=eddy(ied).Rvabs;
        ind2d=eddy(ied).ind2d;

    	% updating I
    	I(ednum)=I(ednum)+1;
    	% Store the data in Eddies
    	Eddies(ednum).data(I(ednum)).x=x;
    	Eddies(ednum).data(I(ednum)).y=y;
    	Eddies(ednum).data(I(ednum)).cx=cx;
        Eddies(ednum).data(I(ednum)).cy=cy;
        Eddies(ednum).data(I(ednum)).clon=clon;
        Eddies(ednum).data(I(ednum)).clat=clat;
        Eddies(ednum).data(I(ednum)).lon=lon;
        Eddies(ednum).data(I(ednum)).lat=lat;
        Eddies(ednum).data(I(ednum)).Curlz_sub=Curlz_sub;
        Eddies(ednum).data(I(ednum)).Rvabs=Rvabs;
        Eddies(ednum).data(I(ednum)).radious=radious;
        Eddies(ednum).data(I(ednum)).ind2d=ind2d;
        Eddies(ednum).data(I(ednum)).time=time;
    end

end

