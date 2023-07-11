function all_idx = reorder_wales(all_res,lt,TR,w)
% reorder wale components according to their position and assign them
% poisitive or negative dmid

% input:
% all_res: all resampled points of all wales
% lt: cell of local triangle indices in the ordering of triangles for all
% wales
% TR: triangulation
% w: stitch width

% output: 
% all_idx: a matrix containing info about all wales points:
%   global wale number, i.e. signed distance from the middle line
%   local node order
%   wale in all_res (wale index, middle line is 1)
%   index of component in all_res
%   local node order within wale with several components
addpath(genpath('..\Extras\geodesic_matlab-master\matlab'))
mesh = geodesic_new_mesh(TR.Points,TR.ConnectivityList);     
alg = geodesic_new_algorithm(mesh, 'exact');

ns = faceNormal(TR);

all_res2 = cell2mat(all_res{1}');
lt2 = cell2mat(lt{1});

all_idx = [];
for i=1:length(all_res{1})
    n1 = size(all_res{1}{i},2);
    sz = size(all_idx,1);
    vec = [zeros(n1,1),(1:n1)',ones(n1,1),i*ones(n1,1),(sz+1:sz+n1)']; 
    all_idx = cat(1,all_idx,vec);
end

% all_res2 = all_res{1}{1};
% lt2 = lt{1}{1};
% n1 = size(all_res{1}{1},2);
% all_idx = [zeros(n1,1),(1:n1)',ones(n1,1),ones(n1,1),(1:n1)'];
% all_idx: [wale number (global), local node order, wale in all_res,
% component in all_res, local node order within wale with several components]
for i=2:length(all_res)
    count1 = 0;
    count2 = 0;
    for j=1:length(all_res{i})
%         if i==1 %&& j==1 % special case handled before the loop
%             continue
%         end
        n1 = size(all_res{i}{j},2);
        destB = cartesianToBarycentric(TR,lt2,all_res2');
        dest_vtx = TR.ConnectivityList(lt2,:);
        Nvtx = size(all_res{i}{j},2);
%         src2 = cell(Nvtx-1,1);
        cl = zeros(Nvtx-1,1);
        dest = zeros(Nvtx-1,1);
        for k=1:Nvtx
%             src2{k} = geodesic_create_surface_point('face',lt{i}{j}(k),all_res{i}{j}(:,k)');
            src = {geodesic_create_surface_point('face',lt{i}{j}(k),all_res{i}{j}(:,k)')};
            geodesic_propagate(alg, src,[],4*w);
            [~, dcurr] = geodesic_distance_and_source(alg);
            destdtri = dcurr(dest_vtx);
            destd = sum(destdtri.*destB,2);
            [dest(k),cl(k)] = min(destd);
        end
        [~,s] = min(dest);
        if length(dest)>1 && s==length(dest)
            [~,s] = min(dest(1:end-1));
        end
        clidx = cl(s);
%         geodesic_propagate(alg, src2,[],4*w);
%         [source_id, dcurr] = geodesic_distance_and_source(alg);
%         destdtri = dcurr(dest_vtx);
%         destd = sum(destdtri.*destB,2);
%         [~,clidx] = min(destd);
%         s = min(source_id(dest_vtx(clidx,:))); 

% decide whether a wale should be prepended or appended in the list of
% wales
% from the middle line, there are two closest wales: -w and w
% we decide on which side are we by computing a cross product (see videos)
% when further from the middle line, we also have two candidates for the
% ordering. We prepend the one with the negative signed distance from the
% middle and append the one with the positive signed distance.
        p2 = all_res2(:,clidx)';       
        if i==2
            p0 = all_res{i}{j}(:,s)';
            p1 = all_res{i}{j}(:,s+1)';
            v3 = ns(lt2(clidx),:);
            v1 = p1-p0;
            v2 = p2-p0;
            v4 = cross(v1,v2);
            v4 = v4/norm(v4);
            crit = v3*v4'<0;
        else
            sz = cellfun(@(x)size(x,2),all_res{i-1});
            
            cmp = find(clidx<=cumsum(sz),1,'first');
            
            row = find(ismember(all_idx(:,2:4),[clidx-sum(sz(1:cmp-1)),i-1,cmp],'rows'));
            wl = all_idx(row(1),1);
            crit = wl<0;
        end
        
%         if norm(v3-v4) > 1e-3 % true if on negative side
        if crit
            numl = sum(all_idx(:,1)==min(all_idx(:,1))-1+sign(count1));
            vec = [(min(all_idx(:,1))-1+sign(count1))*ones(n1,1),...
                (1:n1)',i*ones(n1,1),j*ones(n1,1),(numl+1:numl+n1)'];
            all_idx = cat(1, vec, all_idx);
            count1 = count1 + 1;
        else
            numl = sum(all_idx(:,1)==max(all_idx(:,1))+1-sign(count2));
            vec = [(max(all_idx(:,1))+1-sign(count2))*ones(n1,1),...
                (1:n1)',i*ones(n1,1),j*ones(n1,1),(numl+1:numl+n1)'];
            all_idx = cat(1, all_idx, vec);
            count2 = count2 + 1;
        end
        
    end
    all_res2 = cell2mat(all_res{i}');
    lt2 = cell2mat(lt{i});
end

% visualize points on wales
figure;
trisurf(TR,'facecolor','none','edgecolor',[0.7,0.7,0.7]); 
axis equal
hold on
pts = zeros(length(all_idx),3);
c = zeros(length(all_idx),1);
for i=1:length(all_idx)
    rw = all_idx(i,:);
    pts(i,:) = all_res{rw(3)}{rw(4)}(:,rw(2))';
    c(i) = rw(1);
end
scatter3(pts(:,1),pts(:,2),pts(:,3),30,c,'filled')
