close all; clear; clc
fclose('all');
dr = pwd;

global geodesic_library;                
geodesic_library = 'geodesic_release';
%% import file

[filename,pathobj] = uigetfile('*.obj','Select a .obj file',dr);
if filename == 0 
    disp('No file selected. Nothing imported.')
    clear_interpreter_paths
    fclose('all');
    return 
end

[V,F,UV,TF] = readOBJ(strcat([pathobj,filename]));

%% correct errors in the geometry
if isempty(UV)
    UV = zeros(size(V,1),2);
    TF = F;
end
FTF = [F(:),TF(:)];
FTFsort = sortrows(FTF);
uFTF = unique(FTFsort,'rows');
UVs = UV(uFTF(:,2),:);
%[C,CFF] = connected_components(F);
[SV,SVI,SVJ] = remove_duplicate_vertices(V,1e-5);
SF = SVJ(F);
[~,CFF] = connected_components(SF);
UV_orig = UVs(SVI,:);
TR_start = triangulation(SF,SV);

%% Point click

hg = clickA3DPointDS(SV,SF);
uiwait(hg);
lns = hg.UserData.paths;
ismid = hg.UserData.ismid;
vals = hg.UserData.timevals;
reset(hg)

%% insert points into the triangulation and update triangulation
edges_unique = TR_start.edges();
ln = cell(length(lns), 1);
inserted_vertex_ID = cell(length(lns), 1);
for i=1:length(lns)
    [ln{i},TR_start,inserted_vertex_ID{i}] = insert_seam(lns{i},TR_start,edges_unique);
    edges_unique = TR_start.edges();
end


