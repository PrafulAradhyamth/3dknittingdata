function clear_interpreter_paths

% Extras
eikonal = fullfile(pwd,strcat(['..',filesep,'Extras']),'Eikonal_2D');
geom3d = fullfile(pwd,strcat(['..',filesep,'Extras']),'geom3d-2017.05.05');
grtheory = fullfile(pwd,strcat(['..',filesep,'Extras']),'GrTheory');
meshutil = fullfile(pwd,strcat(['..',filesep],'Extras'),'gptoolbox-master');
fren = fullfile(pwd,strcat(['..',filesep],'Extras'),'frenet_robust');
piecelin = fullfile(pwd,strcat(['..',filesep],'Extras'),'SLMtools');
tb_graph = fullfile(pwd,strcat(['..',filesep],'Extras'),'toolbox_graph');
jsonlab = fullfile(pwd,strcat(['..',filesep],'Extras'),'jsonlab-master');
% opcode = fullfile(pwd,strcat(['..',filesep],'Extras'),'opcodemesh');
% geodesic = fullfile(pwd,strcat(['..',filesep],'Extras'),'geodesic_matlab-master');
% matmesh = fullfile(pwd,strcat(['..',filesep],'Extras'),'matlabmesh');
rmpath(genpath(eikonal))
rmpath(genpath(geom3d))
rmpath(genpath(grtheory))
rmpath(genpath(meshutil))
rmpath(genpath(fren))
rmpath(genpath(piecelin))
rmpath(genpath(tb_graph))
rmpath(genpath(jsonlab))
% rmpath(genpath(opcode))
% rmpath(genpath(geodesic))
% addpath(genpath(matmesh))


