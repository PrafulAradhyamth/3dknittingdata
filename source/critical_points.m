function crit = critical_points(all_res, all_idx, TR_start, w)
    crit = [];
    figure;
    trisurf(TR_start,'facecolor','none','edgecolor',[0.7,0.7,0.7]); 
    axis equal
    hold on
    wales = unique(all_idx(:,1));
    for i=1:length(wales)
        wi = all_idx(all_idx(:,1)==wales(i),:); % all points in a wale
        % get their coordinates
        pts = zeros(size(wi,1),3);
        for j=1:size(wi,1)
            pts(j,:) = all_res{wi(j,3)}{wi(j,4)}(:,wi(j,2))'; 
        end
        plot3(pts(:,1),pts(:,2),pts(:,3),'r')
        if size(wi,1)<3
            continue
        end
        % compute angles, smooth them and get peaks
        angles = mid_angles(1:size(wi,1), pts);
        framelen = round(size(wi,1)/6);
        peak_width = round(framelen/2);
        if rem(framelen,2) == 0
            framelen = framelen + 1;
        end
        order = 3;
        if order > framelen-1 %|| peak_width == 0
            continue
        end
        [~,loc] = findpeaks(sgolayfilt(angles,order,framelen),'MinPeakWidth',peak_width);
        crit = cat(1,crit,pts(loc+1,:));
        %scatter3(pts(loc+1,1),pts(loc+1,2),pts(loc+1,3),'ro','filled')
    end
%     cluster points
    idx = dbscan(crit,5*w,round(0.02*size(crit,1)));
    
    uidx = unique(idx);
    c = jet(length(uidx));
    %c = zeros(length(unique(idx)),1);
%     scatter3(crit(:,1),crit(:,2),crit(:,3),20,jet(length(unique(idx))),'filled')
    for k=1:length(uidx)
        p = crit(idx==uidx(k),:);
        col = repmat(c(k,:),size(p,1),1);
        scatter3(p(:,1),p(:,2),p(:,3),30,col,'filled')
    end
    
    
    