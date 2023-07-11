function [rgbs,symbs,ops] = importCAjson(filename,jsonname)
% import colors for symbols
% prolly boo boo
% input:
% filename: specification file name
% jsonname: name of the json file containing symbols for knitting ops
% output:
% rgbs: rgb values of colors
% smybs: knitting ops symbols
% ops: names of knitting ops
data = loadjson(jsonname);
rel_IDs_all = cellfun(@(x)x.color,data.codes,'UniformOutput',false);
symbs_all = cellfun(@(x)x.sign,data.codes,'UniformOutput',false);
ops_all = cellfun(@(x)x.op,data.codes,'UniformOutput',false);
rel_IDs = cell2mat(rel_IDs_all(~cellfun(@isempty,rel_IDs_all)));
symbs = cell2mat(symbs_all(~cellfun(@isempty,symbs_all)));
ops = ops_all(~cellfun(@isempty,symbs_all));
rgbs = zeros(length(rel_IDs),1,3);
% rgbs2 = zeros(length(rel_IDs),3);
s = xml2struct(filename);
for i=1:length(rel_IDs)
    clr = s.ycTable.yc{rel_IDs(i)}.Attributes.rgb;
%     rgbs2(i,:) = [hex2dec(clr(1:2)),hex2dec(clr(3:4)),hex2dec(clr(5:6))];
    rgbs(i,:,1) = hex2dec(clr(1:2));
    rgbs(i,:,2) = hex2dec(clr(3:4));
    rgbs(i,:,3) = hex2dec(clr(5:6));
end
rgbs = uint8(rgbs);