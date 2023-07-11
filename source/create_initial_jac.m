function jac = create_initial_jac(all_cont,onedge,wvec, idx_first, idx_last)
% based on Liu et al.
% form the initial 2D pattern
% input: 
% all_cont: point structure indicating the points on courses
% onedge: boolean vector indicating if a point is on edge
% wvec: vector of signed distances from the middle
% idx_first: indices of points on the first contour
% idx_last: indices of points on the last contour
% output:
% jacquard structure as matrix
jac = zeros(size(all_cont)-1);
mid = find(wvec==0);
% initial jac: 0=no stitch, 1=regular stitch, 2=left edge, 3=right edge
for j=1:size(all_cont,2)-1
    for i=1:size(all_cont,1)-1
        ptidx = [all_cont(i,j),all_cont(i+1,j),all_cont(i,j+1),all_cont(i+1,j+1)];
        b = [all_cont(i,j) == all_cont(i+1,j), all_cont(i,j+1) == all_cont(i+1,j+1)];
        if b(1) && b(2)
            jac(i,j) = 0;
%         elseif (~b(1) && b(2)) || (b(1) && ~b(2))
%             jac(i,j) = 1;
%             patch(all_pts(ptidx,1),all_pts(ptidx,2),all_pts(ptidx,3),'m')
        else
            if any(onedge(ptidx))
                jac(i,j) = 2;
            else
                jac(i,j) = 1;
            end
        end

    end
end
for i=1:size(jac,1)
    e = find(jac(i,:)>0);
    if length(e)==1
        if e>mid
            jac(i,e-1) = 1;
        else
            jac(i,e+1) = 1;
        end
    end
end
for i=1:size(jac,1)
    e = find(jac(i,:)==2);
    mn = find(jac(i,:)>0,1,'first');
    mx = find(jac(i,:)>0,1,'last');
    for j=1:length(e)
        if e(j)~=mn && e(j)~=mx
            jac(i,e(j)) = 1;
        end
        if e(j)==mn && e(j)>mid
            jac(i,e(j)) = 1;
        end
        if e(j)==mx && e(j)<mid
            jac(i,e(j)) = 1;
        end
        if e(j)==mx && e(j)>mid
            jac(i,e(j)) = 3;
        end
    end
    
end

% fill tiny holes (should look in 3D)
for i=1:size(jac,1)
    for j=2:size(jac,2)-1
        if jac(i,j)==0 && jac(i,j-1)>0 && length(unique(jac(i,j+1:end)))>1 
            % it is a hole
            jac(i,j) = 1; 
        end
    end
end

% correct first and last row (should include the whole boundary)
[~,~,ib] = intersect(idx_first, all_cont(1,:));
jac(1,1:ib(1)-1) = 0;
jac(1,ib(1):ib(end)) = 1;
jac(1,ib(end)+1:end) = 0;
jac(1,ib(1)) = 2;
jac(1,ib(end)) = 3;

[~,~,ib] = intersect(idx_last, all_cont(end,:));
jac(end,1:ib(1)-1) = 0;
jac(end,ib(1):ib(end)) = 1;
jac(end,ib(end)+1:end) = 0;
jac(end,ib(1)) = 2;
jac(end,ib(end)) = 3;

figure;
[r,c] = find(jac==1);
scatter(c,r,'b.'); axis equal; hold on
[r,c] = find(jac==2);
scatter(c,r,'r.'); 
[r,c] = find(jac==3);
scatter(c,r,'m.');
grid on