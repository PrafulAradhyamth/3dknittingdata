function [jac,gleft,gright] = postprocess2D(jac,wvec)
    % correct jacquard to consider the yarn carrier direction
    % input: 
    % jac: jacquard matrix
    % wvec: vector of signed distances from the middle
    % output:
    % jac: updated jacquard as cell array
    % gleft: indices of short rows on the left
    % gright: indices of short rows on the right
mid = find(wvec==0);

%% entferne leere Zeichen - Bugfix by DITF - 2022-10-24
%jacOhneLeereReihen = [];
%for i=1:size(jac,1)    
%    jacZuTesten = find(jac(i,:)>0);
%    if not(isempty(jacZuTesten))
%        jacOhneLeereReihen(end+1, :) = jac(i,:)
%    end    
%end
%jac = jacOhneLeereReihen;

%% convert to cell
jacmat = jac;
jac = cell(size(jacmat,1),1);
gleft = [];
gright = [];

for i=1:size(jacmat,1)
    jac{i} = find(jacmat(i,:)>0);
    if jacmat(i,jac{i}(1))~=2
        gleft = cat(1,gleft,i);
    end
    if jacmat(i,jac{i}(end))~=3
        gright = cat(1,gright,i);
    end
end
N = length(jac);
vec = (1:N)';
norm_edge_left = setdiff(vec,gleft);
norm_edge_right = setdiff(vec,gright);
jacmins = cellfun(@min,jac);
jacmaxs = cellfun(@max,jac);
% workaround for empty norm_edge_left/_right
if isempty(norm_edge_left)
    norm_edge_left = vec; 
    gleft = [];
end
if isempty(norm_edge_right)
    norm_edge_right = vec; 
    gright = [];
end
%% correct increases
% special case 1 on the right side!!!!
df = jac{norm_edge_right(1)}(end) - jac{norm_edge_right(2)}(end);
if df>0
    jac{norm_edge_right(1)} = jac{norm_edge_right(1)}(1:end-df);
end
if df<0
    newmax = jacmaxs(norm_edge_right(1))-df;
    jac{norm_edge_right(1)} = [jac{norm_edge_right(1)},jacmaxs(norm_edge_right(1))+1:newmax];
end

% left
for i=1:2:length(norm_edge_left)-2
    df = jac{norm_edge_left(i)}(1) - jac{norm_edge_left(i+2)}(1);
    if df > 1
        newmin = jacmins(norm_edge_left(i))-df+1;
        jac{norm_edge_left(i)} = [newmin:jacmins(norm_edge_left(i))-1,jac{norm_edge_left(i)}];
        jac{norm_edge_left(i)+1} = [newmin:jacmins(norm_edge_left(i)+1)-1,jac{norm_edge_left(i)+1}];
    end
end

% right
for i=2:2:length(norm_edge_right)-2
    df = jac{norm_edge_right(i+2)}(end) - jac{norm_edge_right(i)}(end);
    if df > 1
        newmax = jacmaxs(norm_edge_right(i))+df-1;
        jac{norm_edge_right(i)} = [jac{norm_edge_right(i)},jacmaxs(norm_edge_right(i))+1:newmax];
        jac{norm_edge_right(i)+1} = [jac{norm_edge_right(i)+1},jacmaxs(norm_edge_right(i)+1)+1:newmax];
    end
end
jacmins = cellfun(@min,jac);
jacmaxs = cellfun(@max,jac);
%% correct decreases
% left
for i=1:length(norm_edge_left)-1
    df = jac{norm_edge_left(i+1)}(1) - jac{norm_edge_left(i)}(1);
    if df > 1 % extend next one
        newmin = jacmins(norm_edge_left(i+1))-df+1;
        jac{norm_edge_left(i+1)} = [newmin:jacmins(norm_edge_left(i+1))-1,jac{norm_edge_left(i+1)}];
    end
end

% right
for i=1:length(norm_edge_right)-1
    df = jac{norm_edge_right(i)}(end) - jac{norm_edge_right(i+1)}(end);
    if df > 1 % extend next one
        newmax = jacmaxs(norm_edge_right(i+1))+df-1;
        jac{norm_edge_right(i+1)} = [jac{norm_edge_right(i+1)},jacmaxs(norm_edge_right(i+1))+1:newmax];    
    end
end
jacmins = cellfun(@min,jac);
jacmaxs = cellfun(@max,jac);
%% extend to min 3 stitches per row
jac_len = cellfun(@length,jac);

too_short_left = find(jacmins<mid & jac_len<3);
too_short_right = find(jacmaxs>mid & jac_len<3);
for i=1:length(too_short_left)
    rowi = too_short_left(i);
    if jac_len(rowi)>3; continue; end
    if rem(rowi,2)==1
        jac{rowi} = jacmins(rowi):jacmins(rowi)+2;
        jac{rowi-1} = jacmins(rowi-1):jacmins(rowi)+2;
    else
        jac{rowi} = jacmins(rowi):jacmins(rowi)+2;
        jac{rowi+1} = jacmins(rowi+1):jacmins(rowi)+2;
    end
end
clear i
for i=1:length(too_short_right)
    rowi = too_short_right(i);
    if jac_len(rowi)>3; continue; end
    if rem(rowi,2)==1
        jac{rowi} = jacmaxs(rowi)-2:jacmaxs(rowi);
        jac{rowi+1} = jacmaxs(rowi)-2:jacmaxs(rowi+1);
    else
        jac{rowi} = jacmaxs(rowi)-2:jacmaxs(rowi);
        jac{rowi-1} = jacmaxs(rowi)-2:jacmaxs(rowi-1);
    end
end

return