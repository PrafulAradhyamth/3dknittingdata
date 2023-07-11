function [resample_pts,level_field_res,res_time] = isocontour_sampling_v03(level_set,ls_cs,w,tri_idx,time)
    % resample contours
    % input: 
    % level_set: level set itself; to resample
    % ls_cs: cumsum of distances in a level set
    % w: sampling rate (stitch width or height)
    % tri_idx: triangle inidces corresponding to ordered nodes in level_set
    % time: time values of level_set
    
    % output:
    % resample_pts: resampled coordinates
    % level_field_res: triangle indices that contain resampled points
    % res_time: time values of resample_pts
w_vec = w:w:max(ls_cs);
resample_pts = zeros(3,length(w_vec));
level_field_res = zeros(length(w_vec),1);
res_time = zeros(length(w_vec),1);
% find a last point which still has lower distance from the start than the 
% specified distance. The next point is already further away. The point
% with the desired distance is somewhere between the two.
for i=1:length(w_vec)
    d_low = find(ls_cs < w_vec(i),1,'last');
    d_up = d_low+1;

    p1 = level_set(:,d_low);
    p2 = level_set(:,d_up);
    t0 = (w_vec(i)-ls_cs(d_low))/(ls_cs(d_up)-ls_cs(d_low));
    p0 = (1-t0)*p1 + t0*p2;
    resample_pts(:,i) = p0;
    
    time1 = time(d_low);
    time2 = time(d_up);
    res_time(i) = (1-t0)*time1 + t0*time2;
    tri1 = tri_idx(d_low,:);
    tri2 = tri_idx(d_up,:);
    triins = intersect(tri1,tri2);
    if ~isempty(triins)
        level_field_res(i) = max(triins); % local triangle
    else
        level_field_res(i) = 0;
    end
    
end