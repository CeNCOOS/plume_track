% grab_HF_data_datepick
%
% script to talk to CORDC to get data for specific dates and length of time
%%
% clear all
% close all

% User inputs
% use 2km data off of SF Bay, 6km beyond SF Bay
res=input('Enter the resolution (2 or 6 km): ');
startdate=input('Enter the start date ([yr,mo,dy,hr,mn,sc]): ');
runlength = input('Enter the number of days to capture: ');

% URL to the data
if res==2
    url='http://hfrnet-tds.ucsd.edu/thredds/dodsC/HFR/USWC/2km/hourly/RTV/HFRADAR_US_West_Coast_2km_Resolution_Hourly_RTV_best.ncd';
end
if res==6
    url='http://hfrnet-tds.ucsd.edu/thredds/dodsC/HFR/USWC/6km/hourly/RTV/HFRADAR_US_West_Coast_6km_Resolution_Hourly_RTV_best.ncd';
end
%urlalt='http://sdf.ndbc.noaa.gov:8080/thredds/dodsC/hfradar_uswc_6km';
%
%url='c:\plume_software\SantaClaraPlumeTraj_example\RTV_HFRADAR_US_West_Coast_2km_Resolution_Hourly_RTV_best.nc';
% need our subset
% max and min latitude 32.5 42
%maxlat=42;
x1=input('Enter north and south boundary: []');
maxlat=x1(1);
minlat=x1(2);
x1=input('Enter west and east boundary: []');
minlon=x1(1);
maxlon=x1(2);
%maxlat=40;
%minlat=37;
%%minlat=32.5;
%% longitude 234 243
%minlon=234-360;
%maxlon=238.8-360;
% get 25 hours of time Note time returned is seconds since 1970
offset=datenum(1970,1,1);
%offset=0;
% current data is too sparse
time1=datenum(startdate);
% time1=now-0/24;
%time1=now-7; % try 1 week back....
% this actually goes 2 weeks back
time2=time1+runlength;
%time2=time1-(25/24);
%time1=((now-offset)*24*60*60)-(5*60*60); % start 5 hours back
%time2=time1-(25*60*60); % go back 25 hours from time1
avgtime=(time1+time2)/2;
%try
    %disp('trying primary')
    ncid=netcdf.open(url);
    talt=0;
%catch
%    disp('trying alternate')
%    ncid=netcdf.open(urlalt);
%    talt=1;
    %keyboard;
%end
% get variables id's and latitude, longitude, and time variables
lonid=netcdf.inqVarID(ncid,'lon');
latid=netcdf.inqVarID(ncid,'lat');
timeid=netcdf.inqVarID(ncid,'time');
toffid=netcdf.inqVarID(ncid,'time_offset');
longitude=netcdf.getVar(ncid,lonid);
latitude=netcdf.getVar(ncid,latid);
time=netcdf.getVar(ncid,timeid);
% The following just returns and array of zeros.
%time_offset=netcdf.getVar(ncid,toffid); %#ok<NASGU>
time_units=netcdf.getAtt(ncid,timeid,'units');
% this will work until 2100 
% parse out the parts we need compute the offset which could change in the
% future
idst=strfind(time_units,'20');
year=str2double(time_units(idst:idst+3));
month=str2double(time_units(idst+5:idst+6));
day=str2double(time_units(idst+8:idst+9));
hour=str2double(time_units(idst+11:idst+12));
minu=str2double(time_units(idst+14:idst+15));
secs=str2double(time_units(idst+17:idst+18));
time_offset=datenum(year,month,day,hour,minu,secs);
% This is the current offset value
%time_offset=datenum(2011,10,01,0,0,0);
time=time/24+time_offset;
%timeo=time;
uid=netcdf.inqVarID(ncid,'u');
vid=netcdf.inqVarID(ncid,'v');
%keyboard;
indlat=find((latitude >= minlat)&(latitude <= maxlat));
indlon=find((longitude>= minlon)&(longitude <= maxlon));
indtime=find((time > time1)&(time <= time2));
%keyboard;
time=time(indtime);
latitude=latitude(indlat);
longitude=longitude(indlon);

start=[indlon(1)-1 indlat(1)-1 indtime(1)-1];

count=[length(indlon) length(indlat) length(indtime)];
% The following won't work since the arrays of u and v are too big to load
% for my computer.
%u=netcdf.getVar(ncid,uid);
%v=netcdf.getVar(ncid,vid);
% This gets a subset of the data
u=netcdf.getVar(ncid,uid,start,count);
v=netcdf.getVar(ncid,vid,start,count);
% get the variable data and convert to cm/s
u=u.*100;
v=v.*100;
netcdf.close(ncid);
[mlon,mlat]=meshgrid(longitude,latitude);
% now we need to flatten arrays...
% need to remove grid points over land...
% need only points over water...
[lx,ly]=size(mlon);
xarr=[];
yarr=[];
%
% load ad simple coast line and outside of bound and sf-bay 
load simple_coast;
xb(411:675)=[];
yb(411:675)=[];
xb(263:308)=[];
yb(263:308)=[];

ll=yb > 40;
yb(ll)=[];
xb(ll)=[];

zind=xb < -122;
nx=xb;
ny=yb;
nx=nx(zind);
ny=ny(zind);
x1=[-123.6 -123 -122.6 -122.4 -122.36];
y1=[39 38.5 38 37.5 37.4];
% simple polygon to blank over land data.
p1=polyfit(y1(1:3),x1(1:3),1);
p2=polyfit(y1(3:end),x1(3:end),1);
lex=length(nx);
mx=[];
my=[];
icount=1;
for i=1:lex
    if ny(i) > 39
        mx(icount)=nx(i); %#ok<SAGROW>
        my(icount)=ny(i); %#ok<SAGROW>
        icount=icount+1;
    else
        if ny(i) < 37.4
            mx(icount)=nx(i); %#ok<SAGROW>
            my(icount)=ny(i); %#ok<SAGROW>
            icount=icount+1;
        else
            if ny(i) > 38
                xjnk=polyval(p1,ny(i));
                if nx(i) < xjnk
                    mx(icount)=nx(i); %#ok<SAGROW>
                    my(icount)=ny(i); %#ok<SAGROW>
                    icount=icount+1;
                end
            else
                if ((ny(i) > 37.76)&&(ny(i) < 37.77))                    
                    if nx(i) < -122.48
                        if nx(i) > -123.08
                            mx(icount)=nx(i); %#ok<SAGROW>
                            my(icount)=ny(i); %#ok<SAGROW>
                            icount=icount+1;
                        end
                    end
                else
                    xjnk=polyval(p2,ny(i));
                    if nx(i) < xjnk
                        mx(icount)=nx(i); %#ok<SAGROW>
                        my(icount)=ny(i); %#ok<SAGROW>
                        icount=icount+1;
                    end
                end
            end
        end
    end
end
ind=find((mx > -122.5)&(my > 39));
mx(ind)=[];
my(ind)=[];
nx=mx;
ny=my;
% reshuffle axis to what we want.
% This sets the axes to x,y,time
u=permute(u,[2 1 3]);
v=permute(v,[2 1 3]);
uarr=[];
varr=[];
zboundx=[];
zboundy=[];
% find the nearest coast boundaries grid points.
for i=1:lx
    dx=abs(ny-mlat(i,1));
    [jnk,kind]=min(dx);
    ind=find(mlon(i,:) < nx(kind));
    if isempty(ind)==0
        %[i length(mlat(i,1)) length(mlon(i,ind(end)))]
%        mlon(i,ind(end))
        zboundy=[zboundy; mlat(i,1)]; %#ok<AGROW>
        zboundx=[zboundx; mlon(i,ind(end))]; %#ok<AGROW>
    end
end
%zboundx(56)=zboundx(55);
%zboundx(58:59)=zboundx(60);
%
kx=[];
ky=[];
% create grid data flattened.
for i=1:lx
    dx=abs(zboundy-mlat(i,1));
    [jnk,kind]=min(dx);
    jind=find(mlon(i,:) < zboundx(kind));
    tmpx=mlon(i,jind);
    tmpy=mlat(i,jind);
    tmpu=squeeze(u(i,jind,:));
    tmpv=squeeze(v(i,jind,:));
    [lux,luy]=size(tmpu);
    if luy > 1
        xarr=[xarr;tmpx']; %#ok<AGROW>
        yarr=[yarr;tmpy']; %#ok<AGROW>
        uarr=[uarr;tmpu]; %#ok<AGROW>
        varr=[varr;tmpv]; %#ok<AGROW>
    end
end

lex=size(uarr,2);
xarr=double(xarr);
yarr=double(yarr);
uarr=double(uarr);
varr=double(varr);
xuse=[]; %#ok<NASGU>
yuse=[]; %#ok<NASGU>
uuse=[];
vuse=[];
uerr=2;
% remove data north and south of grid we want.
ll=yarr > 38.4;
yarr(ll)=[]; xarr(ll)=[];uarr(ll,:)=[];varr(ll,:)=[];
ll=yarr < 37;
yarr(ll)=[]; xarr(ll)=[];uarr(ll,:)=[];varr(ll,:)=[];
ll=xarr > -122.44;
yarr(ll)=[]; xarr(ll)=[];uarr(ll,:)=[];varr(ll,:)=[];
ll=xarr < -123.4;
yarr(ll)=[]; xarr(ll)=[];uarr(ll,:)=[];varr(ll,:)=[];
% we need a filled "grid", so that is what is done below.
for i=1:lex
    ind=find(isnan(uarr(:,i)));
    tmpt=time;
    utmp=uarr(:,i);
    vtmp=varr(:,i);
    xtmp=xarr;
    ytmp=yarr;
    if isempty(ind)==0
        utmp(ind)=[];
        vtmp(ind)=[];
        xtmp(ind)=[];
        ytmp(ind)=[];        
        if length(tmpt) < 2
            % do nothing
        else
          ufill=griddata(xtmp,ytmp,utmp,xarr,yarr);
          vfill=griddata(xtmp,ytmp,vtmp,xarr,yarr);
          uuse=[uuse ufill]; %#ok<AGROW>
          vuse=[vuse vfill]; %#ok<AGROW>
          % we have partial grid that was filled
          %quiver(xarr,yarr,ufill,vfill);
          %drawnow;
          disp('partial');
        end    
    else
        % We have a full grid
        % need to flatten the grid to append?  Already flat...
        uuse=[uuse uarr(:,i)]; %#ok<AGROW>
        vuse=[vuse varr(:,i)]; %#ok<AGROW>
        disp('Full');
    end
end
xuse=xarr;
yuse=yarr;
uarr=uuse;
varr=vuse;
p=double([xarr yarr]);
% define boundary of data points
% The boundary function returns the vector point of indices representing a
% 2-D boundary around the points. 
k=boundary(p,1);
%keyboard;
%set up for 6km totals
if res==6
    k([1:14 48:end])=[];
    k([13 14])=[];
end

% remove open boundary of data set
% below is set up for 2km totals
if res==2
    k([1:45 142:end])=[];
% test to remove some points that are causeing boundary problems
    k([62 65 70 71])=[];
end
% k is now llist below
llist=k;
clear i;
for j=1:length(llist)
    gx_=xarr(llist(j));
    gy_=yarr(llist(j));
    [a,b]=lonlat2km(gx_,gy_,xb,yb);
    c=abs(a+i*b);
    [jnk,ii]=min(c);
    blist(j)=ii; %#ok<SAGROW>
end
% what if we set blist = llist? No because xb and yb is much bigger...
% need to find nearest neighbor?
blist=blist';
ind=find(diff(blist)==0);
blist(ind)=[];
llist(ind)=[];
% original code
[a, b] = lonlat2km(xb(1:end-2), yb(1:end-2), xb(3:end), yb(3:end)); 
thb = atan2(b, a);
thb = [NaN; thb; NaN]; %add NaNs
% new code
%[a,b]=lonlat2km(xb(1:end-1),yb(1:end-1),xb(2:end),yb(2:end));
%thb=atan2(b,a);
%thb=[thb; NaN];
% original code
[a, b] = lonlat2km(xb(2:end-1), yb(2:end-1), xb(3:end), yb(3:end)); 
c = abs(a + i*b);
baxis = cumsum(c); %baxis has 2 less elements than phb xb etc
% new code
%[a,b]=lonlat2km(xb(1:end-1),yb(1:end-1),xb(2:end),yb(2:end));
%c=abs(a+i*b);
%baxis=cumsum(c);
dcf = 0.3;
cv=1e-5*3600;
npts=50;
nlife=24*3;
xc=-122.5726;
yc=37.7836;
gx=xarr;
gy=yarr;
x_=[];
y_=[];
[ns,nt]=size(uarr);
for l=1:length(time)
    t0=clock;
    u_ = uarr(llist,l); v_ = varr(llist, l);    
    mg = abs(u_ + i*v_)*dcf;
    ph = atan2(v_, u_);
    phb = thb(blist);
    mgb = mg.*cos(ph - phb);
    ub = mgb.*cos(phb); vb = mgb.*sin(phb); 
    emg = interp1(baxis(blist), mgb, baxis);
    emg = [NaN; emg; NaN]; %#ok<AGROW>
    % currents are along coast? this appears to be what this does?
    eub = emg.*cos(thb); evb = emg.*sin(thb);
    u_ = uarr(:, l);
    v_ = varr(:, l);        
    ival = find(~isnan(u_));
    jval = find(~isnan(eub));
    % Add (append) new particles (locations) at the source
    x_ = [x_; xc*ones(npts,1)]; y_ = [y_; yc*ones(npts,1)];          %#ok<AGROW>
    % Remove particles beyond nlife age
    if l > nlife
        x_(1:npts) = [];
        y_(1:npts) = [];
    end
    gtmpx=gx(ival);
    gtmpy=gy(ival);
    utmp=u_(ival);
    vtmp=v_(ival);
    ig=find(isnan(gtmpx)==1);
    gtmpx(ig)=[];
    gtmpy(ig)=[];
    utmp(ig)=[];
    vtmp(ig)=[];
    testx=[gtmpx; xb(jval)];
    testy=[gtmpy; yb(jval)];
    testu=[utmp; eub(jval)];
    testv=[vtmp; evb(jval)];
    ind=find(isnan(testu)==1);
    testx(ind)=[];
    testy(ind)=[];
    testu(ind)=[];
    testv(ind)=[];
    % what about advecting data?
    u_ = griddata(testx,testy,testu, x_, y_); 
    v_ = griddata(testx,testy,testv, x_, y_); 
    
    ntotalpt = length(x_);
    th = 2*pi.*rand([ntotalpt 1]);
    un_ = u_ + uerr*cos(th);
    vn_ = v_ + uerr*sin(th);
    mgn_ = abs(un_ + i*vn_); 
    phn = atan2(vn_, un_); 

    dx = un_*cv;
    dy = vn_*cv;
    [xn_, yn_] = km2lonlat(x_, y_, dx, dy); jj = 1; 
    
    for j = 1: ntotalpt
        [cx_, cy_] = polyxpoly(xb, yb, [x_(j) xn_(j)], [y_(j) yn_(j)], 'unique');
        if isempty(cx_)
            continue
        end
        % If an intersection is found, re-compute displacement in
        % along-shore direction
        cx_ = cx_(1);
        cy_ = cy_(1);
        [a, b] = lonlat2km(cx_, cy_, xb, yb);
        c = abs(a + i*b);
        ii = find(c == min(c));
        ii = ii(1);
        dx_ = eub(ii)*cv;
        dy_ = evb(ii)*cv;
        [xn_(j), yn_(j)] = km2lonlat(x_(j), y_(j), dx_, dy_); 
    end
    % Save each timestep to a cell array
    xf{l} = xn_; yf{l} = yn_;  %#ok<SAGROW>
    
    % Update the particle locations for the next iteration
    x_ = xn_; y_ = yn_;
    disp(sprintf('%d of %d, totalpt = %d : time = %4.3f', l, length(time), ntotalpt, etime(clock, t0))); %#ok<DSPS>
end
if numel(xf) <= nlife
    n=numel(xf);
else
    n=nlife;
end
clr=flipud(jet(nlife+1));
load sfcoast.mat;
% This is where we need to get a good set of plots.
for k=2:n
    if k==2
        plot(xb,yb,'k');
        hold on
        plot(zboundx,zboundy,'k.-');
    end  
    disp(['step ',num2str(k),' of ',num2str(n)]);
    for j=k-1:k
%    for j=2:k
%         if j==2
%             plot(xb,yb,'k');
%             hold on
%             plot(zboundx,zboundy,'k.-');
%             %hold on
%         end
    %plot(xf{end}(1+npts*(j-1):npts*j),...
    %    yf{end}(1+npts*(j-1):npts*j),...
    %    '.','color',clr(n-j+1,:));
        x=[];
        y=[];
        x=[x xf{k}(1:50)]; %#ok<AGROW>
        y=[y yf{k}(1:50)]; %#ok<AGROW>
        x=[x xf{k-1}(1:50)]; %#ok<AGROW>
        y=[y yf{k-1}(1:50)]; %#ok<AGROW>
%        for l=k:-1:j           
%            x=[x xf{l}([1+50*(j-2):50*(j-1)])];
%            y=[y yf{l}([1+50*(j-2):50*(j-1)])];
%        end
%        x1=xf{j-1}(1:50);
%        y1=yf{j-1}(1:50);
%        x2=xf{j}(1:50);
%        y2=yf{j}(1:50);
%        x=[x1 x2]';
%        y=[y1 y2]';
         x=x';
         y=y';
         plot(x,y,'.-','color',clr(n-k+1,:));
         set(gca,'xlim',[-123.2 -122.4]);
         set(gca,'ylim',[36.5 38.5]);
         axis equal
%         plot(x,y,'.-','color',clr(n-j+1,:));
        %quiver(xarr,yarr,uarr(:,k),varr(:,k),'m');
%    plot(xf{j}(1:50),...
%        yf{j}(1:50),...
%        '.','color',clr(n-j+1,:));
%    plot(xf{j}(1+npts*(j-1):npts*j),...
%        yf{j}(1+npts*(j-1):npts*j),...
%        '.','color',clr(n-j+1,:));
        hold on
        drawnow;
        %title(datestr(time(k)));
    end
    %if k < n
    %pause;
    %    clf;
    %end
end
%keyboard;
% ulat=unique(yuse);
% lu=length(ulat);
% indx=[];
% indy=[];
% for i=1:lu,
%     ind=find(yuse==ulat(i));
%     [myval,lind]=max(xuse(ind));
%     ky=[ky;ulat(i)]; %#ok<AGROW>
%     kx=[kx;myval]; %#ok<AGROW>
%     indx=[indx;ind(lind)]; %#ok<AGROW>
% end
% [t,xc,yc,npts,nlife,xf,yf]=SFBrmwalk(time,xuse,yuse,uuse,vuse,2);
    

