function [rgbs,symbs] = importCA(filename)

rel_IDs = [1,2,3,4,5,6,7,8,9,10,11,12,13,18,20,21,22,23,28,53];
symbs = '.AYT*I+BGHOWZKMEQLaP';
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