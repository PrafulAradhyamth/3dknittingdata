function write_jacquard_file(jac,fn)
    % write jacquard into file
fileID = fopen(fn,'w');
for i=1:length(jac)
    fprintf(fileID,strcat([jac{i}],'\n'));
end
fclose(fileID);