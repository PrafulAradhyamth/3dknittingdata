function jac_str = generate_ascii_json(jac,gleft,gright,stf,trline_idx,wave_lines,wave_info,jac_patch2)
% generate ASCII code out of 2D Jacquard coordinates
% input: 
% jac: jacquard as cell array
% gleft: short rows on the left
% gright: short rows on the right
% everything else is obsolete
% output:
% jac_str: jacquard as string
dont_knit = 'N';%'.';
default_color = 'A';
stoll_color = 'I';
edge_right = 'A';%'+';
edge_left = 'A';%'B';
dec_right = 'G';
dec_left = 'H';
shortrowedge_right = 'O';%W
shortrowedge_left = 'W';%O
separation_yarn = 'Z';
within_form = '.';%'K';
strobelung = 'M'; % strobel seam
trline = 'L';
wave = 'Q';
clrs = {'A','Y','T','*','I'};


stoll_logo_length = 0;
placeholders = cell(length(jac),1);
starting_line = 2100;
mins = cellfun(@min,jac); maxs = cellfun(@max,jac);
fullrows_left = setdiff((1:length(jac))',[gleft;wave_lines]);
fullrows_right = setdiff((1:length(jac))',[gright;wave_lines]);
jac_str = cell(length(jac),1);

aftergapidx_left = fullrows_left(find(diff(fullrows_left)>1)+1);
aftergapidx_right = fullrows_right(find(diff(fullrows_right)>1)+1);

map_coord2idx = length(min(mins)-1:max(maxs)+1)-max(maxs)-1;
mins_aftergap = mins(aftergapidx_left) + map_coord2idx;
maxs_aftergap = maxs(aftergapidx_right) + map_coord2idx;


for i=1:length(jac)
    jac_str{i} = repmat(dont_knit,1,length(min(mins)-1:max(maxs)+1));
    Q_idx = find(ismember(min(mins)-1:max(maxs)+1,jac{i}));
    if ismember(i,wave_lines)
        jac_str{i}(Q_idx) = wave;
%     elseif ismember(i,trline_idx)
%         jac_str{i}(Q_idx) = trline;
%         jac_str{i}(Q_idx) = default_color;
    else
%         jac_str{i}(Q_idx) = default_color;
%     end
        c2 = unique(jac_patch2{i});
        for j = 1:length(c2)
            regidx = jac_patch2{i}==c2(j);
            jac_str{i}(Q_idx(regidx)) = clrs{c2(j)};
        end
    end
    
    placeholders{i} = [num2str(starting_line+i-1),' ',repmat(dont_knit,1,20)];
end

% Q1 = cellfun(@(x)union(strfind(x,default_color),strfind(x,trline)),jac_str,'UniformOutput',false);
Q1 = cellfun(@(x)find(x~=dont_knit),jac_str,'UniformOutput',0);
Q2 = cellfun(@(x)strfind(x,wave),jac_str,'UniformOutput',false);
clear i
for i=1:length(fullrows_left)
    fulli = fullrows_left(i);
    jac_str{fulli}(Q1{fulli}(1)) = edge_left;
    if ismember(fulli,wave_lines); continue; end
    if i<length(fullrows_left)
        fulli1 = fullrows_left(i+1);
        if mins(fulli1)-1 == mins(fulli)
            jac_str{fulli}(Q1{fulli}(1)) = dec_left;%-1
        end
    end
%     if stf(fulli)&& rem(fulli,32)==0 && length(Q1{fulli})>5
%          jac_str{fulli}(Q1{fulli}(3:4)) = strobelung;
%     end
end
clear i
for i=1:length(gleft)
    jac_str{gleft(i)}(Q1{gleft(i)}(1)) = shortrowedge_left; 
    nextaftergap = find(aftergapidx_left>gleft(i),1,'first');
    firstQidx = Q1{gleft(i)}(1);
    jac_str{gleft(i)}(mins_aftergap(nextaftergap):firstQidx-1) = within_form;
end
clear i
for i=1:length(fullrows_right)
    fulli = fullrows_right(i);    
    if ismember(fulli,wave_lines); continue; end
    jac_str{fulli}(Q1{fulli}(end)) = edge_right;
    if i<length(fullrows_right)
        fulli1 = fullrows_right(i+1);
        if maxs(fulli1)+1 == maxs(fulli)
            jac_str{fulli}(Q1{fulli}(end)) = dec_right;%-1
        end
    end
%     if stf(fulli)&& rem(fulli,32)==0 && length(Q1{fulli})>5
%          jac_str{fulli}(Q1{fulli}(end-3:end-2)) = strobelung;
%     end
end
clear i
for i=1:length(gright)
    jac_str{gright(i)}(Q1{gright(i)}(end)) = shortrowedge_right;
    nextaftergap = find(aftergapidx_right>gright(i),1,'first');
    lastQidx = Q1{gright(i)}(end);
    jac_str{gright(i)}(lastQidx+1:maxs_aftergap(nextaftergap)) = within_form;
end

clear i
for i=1:length(wave_lines)
    nextaftergap = find(aftergapidx_left>wave_lines(i),1,'first');
    firstQidx = Q2{wave_lines(i)}(1);
    jac_str{wave_lines(i)}(mins_aftergap(nextaftergap):firstQidx-1) = within_form; 
    
    nextaftergap = find(aftergapidx_right>wave_lines(i),1,'first');
    lastQidx = Q2{wave_lines(i)}(end);
    jac_str{wave_lines(i)}(lastQidx+1:maxs_aftergap(nextaftergap)) = within_form;
    % some boundaries have 'K'. Correct them
%     jac_str{wave_lines(i)}(firstQidx) = wave;
%     jac_str{wave_lines(i)}(lastQidx) = wave;
end

idxtochange = Q1{end};
xtrarows = cell(stoll_logo_length+1,1);
for i=1:length(xtrarows)
    xtrarows{i} = repmat(dont_knit,1,length(min(mins)-1:max(maxs)+1));
    xtrarows{i}(idxtochange) = stoll_color;
    placeholders{i+length(jac)} = [num2str(starting_line+i-1),' ',repmat(dont_knit,1,20)];
end
xtrarows{1}(idxtochange) = separation_yarn;
jac_str = [jac_str; xtrarows];
% jac_str = flipud(jac_str);   
placeholders{end}(1:4) = '9999';
jac_str = [placeholders,flipud(jac_str)];
    


