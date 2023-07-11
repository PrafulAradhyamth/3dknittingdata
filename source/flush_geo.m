function fig = flush_geo(fig)

val = -inf;    
% while abs(val) > 1
%     try
%         val = input('Specify the time value between -1 and 1: ');
%         validateattributes(val,{'numeric'},{'>=',-1,'<=',1});
%     catch
%         disp('Invalid input!')
%     end
% end


disp('Draw new geodesic. Flushing existing user data.')

trist = get(fig,'UserData');
trist.clicked = 0;
paths = [];
for i=2:length(trist.phandle)
    set(trist.phandle{i},'markerfacecolor','b') 
end
for i=3:length(trist.lhandle)
    set(trist.lhandle{i},'color','b')
    usr = trist.lhandle{i}.UserData;
    x = trist.lhandle{i}.XData;
    y = trist.lhandle{i}.YData;
    z = trist.lhandle{i}.ZData;
    id = usr(:,1);
    type = usr(:,2);
    partpath = [x(:),y(:),z(:),id(:),type(:)];
    paths = cat(1,paths,partpath);
end
trpaths = trist.paths;
ismid = trist.ismid;
trpaths = cat(1,trpaths,paths);
ismid = cat(1,ismid,trist.drawmid);
tvals = cat(1,trist.timevals,val);
trist.val = val;
trist.ismid = ismid;
trist.paths = trpaths;
trist.timevals = tvals;
trist.phandle = {NaN};
trist.lhandle = {NaN,NaN};
trist.lendpts = [NaN,NaN;NaN,NaN];
set(fig,'UserData',trist);