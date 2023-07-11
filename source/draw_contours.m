function [ls, res_temp, local_tri,ltime] = draw_contours(lset,tri_idx,c,h,level_time)
% resample, order and draw contours with respect to input scalar field
% input: 
% lset: 3D coordinates for a certain level
% tri_idx: index of triangles which contain edges with lset
% c: point color for plot
% h: stitch height
% level_time: time values of lset

% output:
% ls: reordered but not resampled coordinates
% res_temp: resampled coordinates
% local_tri: local triangle index in the ordering of triangles
% ltime: time values of res_temp

[~,~,ti3] = unique(tri_idx);
ti4 = reshape(ti3,numel(ti3)/2,2);
G = graph(ti4(:,1),ti4(:,2));
concom = conncomp(G,'outputform','cell');
valb = zeros(size(ti4,1),1);
for j=1:size(ti4,1)
    valb(j) = findedge(G,ti4(j,1),ti4(j,2)); 
end
res_temp = [];
local_tri = [];
ls = [];
t1 = [];
ltime = [];
for j=1:length(concom)
    [ls_reorder,~] = contour_forming2(concom{j},G,valb);
    ls_ordered = lset(:,ls_reorder);
    plot3(ls_ordered(1,:),ls_ordered(2,:),ls_ordered(3,:),c); axis equal; hold on;
    drawnow
    pause(0.0001)
    
    if h~=0
        t_ordered = level_time(ls_reorder);
        tidx = tri_idx(ls_reorder,:);
        if t_ordered(1) > t_ordered(end)
            ls_ordered = fliplr(ls_ordered);
            t_ordered = flipud(t_ordered);
            tidx = flipud(tidx);
            ls_reorder = flipud(ls_reorder); 
        end
    % compute cumulative sum of distances from a starting point
        B = diff(ls_ordered,1,2); 
        ls_temp = sqrt(diag(B'*B));
        ls_temp_cs = [0;cumsum(ls_temp)];
        cont_len = ls_temp_cs(end)-1e-8; 
        k0 = round(cont_len/h); 
        hopt = cont_len/k0;
        [resample_pts,lt,res_time] = isocontour_sampling_v03(ls_ordered,ls_temp_cs,hopt,tidx,t_ordered);
        resample_pts = cat(2,ls_ordered(:,1),resample_pts);
        res_time = cat(1,t_ordered(1),res_time);
%         if size(resample_pts,2)<2
%             continue
%         end
        lt = cat(1, max(tri_idx(ls_reorder(1),:)), lt);
%         scatter3(resample_pts(1,:),resample_pts(2,:),resample_pts(3,:),'ro','filled') 
%         scatter3(resample_pts(1,1),resample_pts(2,1),resample_pts(3,1),'ko','filled')
        if isempty(res_temp) || sum(ismember(resample_pts',cell2mat(res_temp(:)')','rows'))==0
                res_temp = cat(1,res_temp,{resample_pts});
                local_tri = cat(1,local_tri,{lt});
                ltime = cat(1,ltime,{res_time});
                ls = cat(1,ls,{ls_ordered});
                t1 = cat(1,t1,t_ordered(1));
            
        end
    end
end
[~,idx] = sort(t1);
ls = ls(idx);
if h~=0
    res_temp = res_temp(idx);
    local_tri = local_tri(idx);
    ltime = ltime(idx);
end

