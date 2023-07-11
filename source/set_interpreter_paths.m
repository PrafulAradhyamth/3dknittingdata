function set_interpreter_paths

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
addpath(genpath(eikonal))
addpath(genpath(geom3d))
addpath(genpath(grtheory))
addpath(genpath(meshutil))
addpath(genpath(fren))
addpath(genpath(piecelin))
addpath(genpath(tb_graph))
addpath(genpath(jsonlab))
% addpath(genpath(opcode))
% addpath(genpath(geodesic))
% addpath(genpath(matmesh))


