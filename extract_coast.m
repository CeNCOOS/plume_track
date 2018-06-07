% extract coast data for region
load coastOR2Mex
maxlat=input('Enter the maximum latitudinal extent to get: ');
minlat=input('Enter the minimum latitudinal extent to get: ');
inde=find(OR2Mex(:,2) < maxlat);
stop=min(inde);
%stop=inde(end);
inds=find(OR2Mex(:,2) > minlat);
start=max(inds);
%ind=find((OR2Mex(:,2) >= minlat)&(OR2Mex(:,2) < maxlat));
if start > stop,
    tmp=start;
    start=stop;
    stop=tmp;
end
xb=OR2Mex(start:stop,1);
yb=OR2Mex(start:stop,2);
xb=xb';
yb=yb';
fname=input('Enter the coastline file name to save to: ','s');
eval(['save ',fname,'.mat xb yb']);