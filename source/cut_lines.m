function [TR_start,edges_unique,seamidx2,lb,ub,vint,seamtm,mididx2] = cut_lines(TR_start,seamidx,edges_unique,lb,ub,vint,mididx)
    % cut at seam, i.e. double the points at the seam and separate them
    % topologically
    % input:
    % TR_start: triangulation, i.e. mesh
    % seamidx: indices of points on the seam
    % edges_unique: edges in TR_start
    % lb: indices of the points on the lower boundary of the knitting time
    % ub: indices of the points on the upper boundary of the knitting time
    % vint: time field, i.e. time values for all points in the mesh
    % mididx: indices of points on the middle line
    
    % output:
    % TR_start: updated triangulation
    % edges_uniqeu: updated edges
    % seamidx2: updated indices of seam points
    % lb, ub, vint: updated lower, upper boundary and time values, resp.
    % seamtm: time values of seam points
    % mididx2: updated indices of middle line points
    
    seamedges = ismember(edges_unique,seamidx);
    seamedg = sum(seamedges,2)==2;
    [Fs,I] = cut_edges(TR_start.ConnectivityList,edges_unique(seamedg,:));
    Vs = TR_start.Points(I,:);
    try
        FF = flip_ears(Vs,Fs);
    catch
        FF = Fs;
    end
    TR_start = triangulation(FF,Vs);
    TR_start = remedy_mesh(TR_start);% not working properly? remap?
    edges_unique = TR_start.edges();
    % remap nodes on the seam 
    seamidx2 = find(ismember(I,seamidx));
    ub = find(ismember(I,ub));
    lb = find(ismember(I,lb));
    vint = vint(I);
    seamtm = vint(seamidx2);    
%     mididx2 = find(ismember(I,mididx));
    mididx2 = zeros(length(mididx),1);
    for i=1:length(mididx)
        mididx2(i) = find(I==mididx(i));
    end
%     [lia,locb] = ismember(I,mididx);
%     locb = locb(locb>0);
%     liaidx = find(lia);
%     mididx2 = liaidx(locb);