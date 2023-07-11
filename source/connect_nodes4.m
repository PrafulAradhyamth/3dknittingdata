function [all_cont,onedge,wvec, all_time, all_pts] = connect_nodes4(all_idx,all_res,w,ltime)
% connect nodes starting from the mid
% input:
% all_idx: sampled points structure
% all_res: coordinates of all sampled points on wales
% w: stitch width
% ltime: time values of all_res
% output:
% all_cont: matrix with point indices, indicating the structure of knitting
% courses
% onedge: boolean whether a point is on edge or not
% wvec: vector of all signed distances from middle line
% all_time: all time values for all points
% all_pts: all point coordinates
all_pts = zeros(size(all_idx,1),3);
all_time = zeros(size(all_idx,1),1);
onedge = false(size(all_idx,1),1);
for i=1:size(all_idx,1)
    i0 = all_idx(i,:);
    all_pts(i,:) = all_res{i0(3)}{i0(4)}(:,i0(2))'; 
    all_time(i) = ltime{i0(3)}{i0(4)}(i0(2));
    if ( i0(2)==1 || i0(2)==size(all_res{i0(3)}{i0(4)},2) ) && abs(all_time(i)) < 1-1e-3
        onedge(i) = true;
    end
end

figure;
scatter3(all_pts(:,1),all_pts(:,2),all_pts(:,3),'r.')
axis equal
hold on
for i=1:length(all_res)
    for j=1:length(all_res{i})
        plot3(all_res{i}{j}(1,:),all_res{i}{j}(2,:),all_res{i}{j}(3,:),'r')
    end
end

maxwale = max(all_idx(:,1));
minwale = min(all_idx(:,1));
wvec = minwale:maxwale;
all_cont = [];
all_local = [];
all_visit = [];

notvisited = true(size(all_idx,1),1);


while any(notvisited)
    % find the smallest (by abs value) non-visited wale and the smallest node there 
    nvidx = find(notvisited);
    nonvisitedwale = all_idx(nvidx,1);
    [~,minwaleidx] = min(abs(nonvisitedwale)); % always finds first which is also 
    % the one with the smallest time, i.e. the smallest component
    minabswale = nonvisitedwale(minwaleidx);
    st_idx = nvidx(minwaleidx);
    notvisited(st_idx) = false;
    current_cont = zeros(1,length(wvec));
    local_visit = false(1,length(wvec));
    local_visit(wvec==minabswale) = true;
    current_cont(wvec==minabswale) = st_idx;
    local_course = zeros(1,length(wvec));
    local_course(wvec==minabswale) = all_idx(st_idx,5);

    % build contour to the right
    [notvisited,current_cont,local_course,local_visit] = ...
        build_cont(minabswale,maxwale,1,all_idx,all_pts,st_idx,...
        notvisited,all_cont,local_course,all_local,current_cont,wvec,...
        all_time,local_visit,onedge,w);

    st_idx = nvidx(minwaleidx);

    % build contour to the left
    [notvisited,current_cont,local_course,local_visit] = ...
        build_cont(minabswale,minwale,-1,all_idx,all_pts,st_idx,...
        notvisited,all_cont,local_course,all_local,current_cont,wvec,...
        all_time,local_visit,onedge,w);
    
    % placement of contours immediately after they are built
    % find all visited nodes in all_cont which lay entirely on the same 
    % wales as new nodes in current cont
    if sum(local_visit)<2
        continue
    end
    if ~isempty(all_cont)    
        lcidx = find(local_course>0);
        lc = local_course(local_course>0);
        lastprev = zeros(length(lc),1);
        for i=1:length(lc)
            lpi = find(all_local(:,lcidx(i))>0 & all_local(:,lcidx(i))<lc(i),1,'last'); 
            if ~isempty(lpi)
                lastprev(i) = lpi;
            end
        end
        all_cont = cat(1, all_cont, current_cont)
        all_local = cat(1,all_local,local_course)
        all_visit = cat(1,all_visit,local_visit)
        pl = max(lastprev)+1;
%         pl = C(si(si==size(all_temp,1)));
        all_cont = [all_cont(1:pl-1,:); all_cont(end,:); all_cont(pl:end-1,:)]
        all_local = [all_local(1:pl-1,:); all_local(end,:); all_local(pl:end-1,:)]
        all_visit = [all_visit(1:pl-1,:); all_visit(end,:); all_visit(pl:end-1,:)]
    else
        all_cont = current_cont
        all_local = local_course
        all_visit = local_visit
    end
    plot3(all_pts(current_cont(local_visit),1),...
        all_pts(current_cont(local_visit),2),...
        all_pts(current_cont(local_visit),3),'b')
    drawnow;
%     pause

   
end


% figure; spy(flipud(all_visit))

% create graph
% wales
% G = graph;
% for j=1:size(all_cont,2)
%     nz = unique(all_cont(all_visit(:,j),j),'stable');    
%     weights = vecnorm(all_pts(nz(1:end-1),:) - all_pts(nz(2:end),:),2,2);
%     EdgeTable = table(weights,false(length(nz)-1,1),true(length(nz)-1,1),...
%         'VariableNames',{'Weights','Course','Wale'});
%     G = addedge(G,nz(1:end-1),nz(2:end),EdgeTable);
% end
% 
% for i=1:size(all_cont,1)
%     nz = unique(all_cont(i,all_visit(i,:)),'stable');
%     weights = vecnorm(all_pts(nz(1:end-1),:) - all_pts(nz(2:end),:),2,2);
%     EdgeTable = table(weights,true(length(nz)-1,1),false(length(nz)-1,1),...
%         'VariableNames',{'Weights','Course','Wale'});
%     G = addedge(G,nz(1:end-1),nz(2:end),EdgeTable);
% end
% pad all_cont
for j=1:size(all_cont,2)
    for i=1:size(all_cont,1)
        idx = find(all_cont(:,j)>0); 
        if all_cont(i,j)==0
            [~,idx2] = min(abs(idx-i));
            all_cont(i,j) = all_cont(idx(idx2),j); %#ok<AGROW>
        end
    end
end

end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [notvisited,current_cont,local_course,local_visit] = ...
        build_cont(minabswale,endwale,step,all_idx,all_pts,st_idx,...
        notvisited,all_cont,local_course,all_local,current_cont,wvec,...
        all_time,local_visit,onedge,w)
%     st_idx1 = st_idx;
    edge_reached = false;
    for i=minabswale:step:endwale-step 
        nextwale = find(all_idx(:,1) == i+step);
        
        % prefer non visited ones unless there are none
        notvis = nextwale(notvisited(nextwale));
%         currentwale = find(all_idx(:,1) == i);
        nextwpts = all_pts(nextwale,:);
        
        stpt = all_pts(st_idx,:);
        d1 = pdist2(nextwpts,stpt);
        [d1s,si] = sort(d1);
        ei = find(d1s<2*w,1,'last');
        if isempty(ei)
            break
        end
%         ei = min(length(si),4);
        closest_candidate_t = all_time(nextwale(si(1:ei)));
%         stpt_time = all_time(st_idx1);
        stpt_time = all_time(st_idx);
%         stpt_time = mean(all_time(current_cont(current_cont>0)));
        [~,i1] = min(abs(closest_candidate_t-stpt_time));
        idx1 = si(i1);    
%         idx1 = si(1);
        nextpt = nextwale(idx1);
        
        % alternative; just the closest
%         nextpt = nextwale(si(1));
%         if onedge(nextpt) || d1(idx1) > 4*w
%             edge_reached = true;
%         end
        current_cont(wvec==i+step) = nextpt;
        local_course(wvec==i+step) = all_idx(nextpt,5);
        local_visit(wvec==i+step) = true;
        if ~notvisited(nextpt) || (onedge(nextpt))% && ~onedge(current_cont(wvec==minabswale))) %|| d1(idx1) > 3*w
            notvisited(nextpt) = false;
            break
        end
        
        notvisited(nextpt) = false;

        st_idx = nextpt;
    end    
   
end



