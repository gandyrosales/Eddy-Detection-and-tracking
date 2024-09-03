function eddy = DetectEddyOW(usurf,vsurf,usub,vsub,loni,lati,xi,yi,dx,dy,maskr,OW_smooth,crow,crar,plotflag)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takeyoshi Nagai@UMassD 12/9/2010 -- Applied in the California Current System
% see Nagai et al. 2015. https://doi.org/10.1002/2015JC010889
% In this modified version, we apply this algorithm in the Peru-Chile EBUS
% % in Rosales-Quintana et al -- 2024, september.
%
% Detection of eddies
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Ar=dx.*dy;

% conmputing OW and vorticity
[Vort_sub,OW_sub,~,Rabs] = CCScompCurl(usub,vsub,dx,dy);
[Vort_surf,OW_surf,~,Rabs] = CCScompCurl(usurf,vsurf,dx,dy);

% OW and vorticity subsuface
OW_sub_smooth = movmean(OW_sub,OW_smooth,1,'omitnan');
OW_sub_smooth = movmean(OW_sub_smooth,OW_smooth,2,'omitnan');
Vort_sub_smooth = movmean(Vort_sub,OW_smooth,1,'omitnan');
Vort_sub_smooth = movmean(Vort_sub_smooth,OW_smooth,2,'omitnan');

% OW and vorticity surface
OW_surf_smooth = movmean(OW_surf,OW_smooth,1,'omitnan');
OW_surf_smooth = movmean(OW_surf_smooth,OW_smooth,2,'omitnan');
Vort_surf_smooth = movmean(Vort_surf,OW_smooth,1,'omitnan');
Vort_surf_smooth = movmean(Vort_surf_smooth,OW_smooth,2,'omitnan');


% contour at critical value of OW subsurface target
[c] = contourc(xi(1,:),yi(:,1),double(OW_sub_smooth).*maskr,[-abs(crow) -abs(crow)]);
tnum=size(c,2);
num=c(2,1);
nk=1;

%-------------------------------------------------------------
%-------------------------------------------------------------
% extract all closed contour
i=0;
flag=1;
while flag==1
    x=c(1,nk+1:nk+num);
    y=c(2,nk+1:nk+num);
    if x(:,1) == x(:,end) & y(:,1) == y(:,end) % when contour is closed
        i=i+1;
        eddyall(i).x=x;
        eddyall(i).y=y;
        eddyall(i).lon=interp2(xi,yi,lon,eddyall(i).x,eddyall(i).y);
        eddyall(i).lat=interp2(xi,yi,lat,eddyall(i).x,eddyall(i).y);

    end
    nk=nk+num+1;
    if nk>tnum
        flag=-1;
    else
        num=c(2,nk);
    end
end

%-------------------------------------------------------------
%-------------------------------------------------------------
% for plotting results
if plotflag==1
    set(gcf, 'position',[100 1300 900 900])
    colormap(colormapcat('bck','bk','b','w','r','rk','rkm'))
    pcolor(loni,lati,OW_sub_smooth.*maskr)
    shading flat
    colorbar
    caxis([-2e-10 2e-10])
    ylim([-20 -10])
    yticks([-20:5:-10])
    xlim([-85 -70])
    xticks([-85:5:-70])
    set(gca,'ydir','normal')
    drawnow
    hold on
end

%-------------------------------------------------------------
% screening eddies
%-------------------------------------------------------------


ieok=0;
for ie= 1:1:size(eddyall,2)
    if mod(ie,100)==0
        disp(['screening eddy #' num2str(ie) ' of ' num2str(size(eddyall,2))])
    end
    XX=eddyall(ie).x(1:end-1);
    YY=eddyall(ie).y(1:end-1);

    if length(XX)>tnum_critic && length(YY)>tnum_critic
        tri = delaunayn([XX(:) YY(:)]);
        warning off
        T = tsearchn([XX(:) YY(:)], tri, [xi(:) yi(:)]); % K is the total index
        ind2d = find(~isnan(T));
        Arsum=sum(Ar(ind2d)); % area within eddiess
        radius = sqrt(Arsum/pi)./1000; % radious in Km

        % threshold
        eqrud = sqrt(Arsum./pi);
        lone=eddyall(ie).lon;
        late=eddyall(ie).lat;
        dr = sum(compdist(lone(1:end-1),late(1:end-1),lone(2:end),late(2:end)));
        eqrud2 = dr./2./pi;
        ratioRudi = eqrud2./eqrud;

        Curlz_surf=mean(Vort_surf_smooth(ind2d));
        Curlz_sub=mean(Vort_sub_smooth(ind2d));

        % criteria for Puddies
        if crar(1) <=Arsum && Arsum <= crar(2) ...
                && (abs(Curlz_sub)>abs(Curlz_surf))...
                && (0.95 <= ratioRudi)  && (ratioRudi <= 1.2)
            ieok=ieok+1;
            % only puddies information will be saved
            % compute eddy center variables
            cix=mean(ix(ind2d));
            ciy=mean(iy(ind2d));
            cx=mean(xi(ind2d));
            cy=mean(yi(ind2d));
            clon=mean(loni(ind2d));
            clat=mean(lati(ind2d));
            Curlz_sub=mean(Vort_sub_smooth(ind2d)); % 
            Rv= Rabs(ind2d);

            % storing eddy information to a struct named "eddy"
            eddy(ieok).ix=eddyall(ie).ix;
            eddy(ieok).iy=eddyall(ie).iy;
            eddy(ieok).x=eddyall(ie).x;
            eddy(ieok).y=eddyall(ie).y;
            eddy(ieok).lon=eddyall(ie).lon;
            eddy(ieok).lat=eddyall(ie).lat;
            eddy(ieok).cix=cix;
            eddy(ieok).ciy=ciy;
            eddy(ieok).cx=cx;
            eddy(ieok).cy=cy;
            eddy(ieok).clon=clon;
            eddy(ieok).clat=clat;
            eddy(ieok).ind2d=ind2d;
            eddy(ieok).Curlz_sub=Curlz_sub;
            eddy(ieok).radius=radius;
            eddy(ieok).Rv=Rv;

        end

    end
end


%-------------------------------------------------------------
%-------------------------------------------------------------
% for plotting results
if plotflag==1
    set(gcf, 'position',[300 1200  1000 800])
    pcolor(lon-360,lat,OW_sub_smooth.*maskr)
    colormap(colormapcat('bck','bk','b','w','r','rk','rkm'))
    shading flat
    colorbar
    caxis([-1e-10 1e-10])
    hold on
    for i=1:1:size(eddy,2)
        plot(eddy(i).lon,eddy(i).lat,'y-','LineWidth',2);
        plot(eddy(i).clon,eddy(i).clat,'cx','markersize',15,'LineWidth',3);
        grid on
        drawnow
    end
end



