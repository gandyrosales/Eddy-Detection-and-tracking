function eddy=CCSEddyTracker(peddy,eddy,dx,dy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takeyoshi Nagai@UMassD 12/9/2010 -- Applied in the California Current System
% see Nagai et al. 2015. https://doi.org/10.1002/2015JC010889
% In this modified version, we apply this algorithm in the Peru-Chile EBUS.
% % in Rosales-Quintana et al -- 2024, september.
% Tracking function of detected eddies
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%-------------------------------------------------
% starting the main loop
%-------------------------------------------------
Ar=dx.*dy;
numped=size(peddy,2);
numed=size(eddy,2);
flagfirst=0;
if ~isfield(peddy,'ednum')
    flagfirst=1;
end
% looking for the same eddy
maxednum=-999;
idassigned=[];
for iped=1:1:numped
    if flagfirst
        peddy(iped).ednum=iped;
    end
    ind2dp=peddy(iped).ind2d;
    cxp=peddy(iped).cx;
    cyp=peddy(iped).cy;
    pednum=peddy(iped).ednum; % eddy number
    for ied=1:1:numed
        ind2d=eddy(ied).ind2d;
        cx=eddy(ied).cx;
        cy=eddy(ied).cy;
        inds=ind2d(ismember(ind2d,ind2dp));
        numshare(iped,ied)=sum(Ar(inds));
        dr(iped,ied)=sqrt((cxp-cx)^2+(cyp-cy)^2);
    end
    maxnum=max(numshare(iped,:));
    mindist=min(dr(iped,:));
    if maxnum~=0
        % new eddy number is defined by the maximum area sharing eddy.
        % But there is a posibility that multiple eddies in a new time step
        % have the same area of the grid. This will be taken care by
        % using distance closest eddy criteria below.
        eddyid=find(maxnum==numshare(iped,:));

        if length(eddyid)>1
            disp('multiple eddies')
            eddyid2=find(min(dr(iped,eddyid))==dr(iped,eddyid));
            eddyid=eddyid(eddyid2);
        end
        eddyid = eddyid(1);
        eddy(eddyid).ednum=pednum;
        idassigned=cat(2,idassigned,eddyid);
        maxednum=max([maxednum pednum]);
    end

    % taking care eddies not assigned any number yet.
    eddyorgids=[1:1:numed];
    eddyrest=eddyorgids(~ismember([1:1:numed],idassigned));
    for iedrest=eddyrest
        maxednum=maxednum+1;
        eddy(iedrest).ednum=maxednum;
    end



