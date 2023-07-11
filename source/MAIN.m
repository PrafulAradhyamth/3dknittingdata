%% main interpreter file:

% Date: 25 Oct 2019
% Author: Dominik Surc

% generate a knitting pattern from .obj triangular mesh files
% input: obj file, stitch parameters (must be set in this file)
% output: knitting pattern in jac and bmp format

% Funktionen von Dominik angepasst oder erstellt DS im Namen

close all; clear; clc
fclose('all');
dr = pwd;

%% stitch parameters
% possible parameters (width, height) for the following knit constructions:
% rechts-rechts (RR), Halbschlauch (HS), italienisches Patent (IP), 
% rechts-links (RL), Halbschlauch 2 (HS2)
% ask Uwe for details :)
W = [0.96, 0.88, 0.84, 0.76, 10/6, 50/6, 1, 10, 40];
H = [0.42, 0.3, 0.85, 0.44, 10/13, 50/13, 0.5, 5, 20];
KNIT_STR = {'RR', 'HS', 'IP', 'RL', 'HS2', 'Test1', 'Test2','Test3'};
KNIT = 9; % change here!
h = H(KNIT); %Maschenhöhe
w = W(KNIT); %Maschenweite

%% set paths, mainly files from folder "extras"
set_interpreter_paths

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

[V,F] = readOBJ(strcat([pathobj,filename]));

% perform some sanity checks. A mesh shouldn't have any duplicate nodes and
% have no free hanging triangles
[SV,SVI,SVJ] = remove_duplicate_vertices(V,1e-5);
SF = SVJ(F);
[~,CFF] = connected_components(SF);

% generate initial triangulation
TR_start = triangulation(SF,SV);

%% GUI 
% see readme for usage

% order of lines: upper bnd, lower bnd, seam, one or several middle lines

% Oeffnet die Figure und erlaubt das setzen von Linien (Funktion fcw
hg = clickA3DPointDS(SV,SF);
uiwait(hg);

%% alle Linien, die vom Benutzer markiert wurden
lns = hg.UserData.paths; % clicked lines

%% Elemente innerhalb von lns
% Array mit 5 Elementen
% 1-3 x,y,z
% 4 Index von Knoten, bzw. Kante (hängt von Spalte 5 ab)
% 5 Knoten = 1, Kanten = 2

%% Linie, die als Mittellinie markiert wurde
ismid = hg.UserData.ismid; % boolean vector labeling lines, true for middle line

reset(hg)
disp('moving forward')

%% check if we have a seam
% usrseam - Naht
usrseam = false;
if length(lns)>2 && ~ismid(3); usrseam = true; end

% define edge structure
edges_unique = TR_start.edges();

% insert seam points into the mesh
% midpts - Punkte der Mittellinie
midpts = cell2mat(lns(logical(ismid)));
if usrseam
    % wenn Naht vorhanden, wird Naht und Mittellinie inslines zugewiesen
    inslines = [lns{3};midpts];
    % Länge (Anzahl Punkte) der Naht
    seamlen = size(lns{3},1);
else
    % ohne Naht, wird Mittellinie inslines zugewiesen
    inslines = midpts;    
    seamlen = 0;
end

%% wenn Naht oder Mittellinie vorhanden (wenn inslinies ist nicht leer)
if ~isempty(inslines)
    %% - % Naht wird eingefügt
    [midseam,TR_start,midseamidx] = insert_seam(inslines,TR_start,edges_unique);
    
    edges_unique = TR_start.edges();
    
    mid = midseam(seamlen+1:end,:);    
    mididx = midseamidx(seamlen+1:end);
    seam = midseam(1:seamlen,:);
    seamidx = midseamidx(1:seamlen);
    
    mid = unique(mid,'rows','stable');
    mididx = unique(mididx,'stable');
    seam = unique(seam,'rows','stable');
    seamidx = unique(seamidx,'stable');
else
    mid = [];
    mididx = [];
    seam = [];
    seamidx = [];
end

%% Figure mit eingezeichnter Mittelline / Naht zeichnen
figure; trisurf(TR_start); axis equal; hold on

%% fill holes
% contours of the knitting time are perpendicular to mesh boundary.
% holes (=boundary) in the mesh cause distortions in the knitting time.
% therefore, we temporarely fill them with triangles.

HF = fill_holesDS(TR_start.Points,TR_start.ConnectivityList);
[C,ia,ic] = unique(HF);
Cmat = reshape(ic,size(HF,1),size(HF,2));
HV = TR_start.Points(C,:); % !!! probably only ok to continue if no nodes are removed
[VV,FF,I] = remove_degenerate_faces_DS(HV,Cmat,'Epsilon',1e-6);
FI = I(FF);
CF = C(FI);
ConnectivityListWithFilledHoles = [TR_start.ConnectivityList;CF];

%% get laplace interpolation on filled mesh
% laplace interpolation computes the knitting time,
% we use (only here!) the filled triangle structure
lap = cotmatrix_embedded(TR_start.Points,ConnectivityListWithFilledHoles);

if unique(lns{1}(:,5))~=1
    error('Edge in upper boundary!')
end
if unique(lns{2}(:,5))~=1
    error('Edge in lower boundary!')
end

% define boundary values for laplace interpolation
% lower boundary gets value -1
% upper boundary gets value 1
% can be freely chosen
ub = unique(lns{1}(:,4),'stable');
lb = unique(lns{2}(:,4),'stable');
bn = [ub;lb];
f2 = ones(length(bn),1);
f2(1:size(ub,1)) = -1;

% Strickzeiten aller Knoten werden in vint gespeichert
vint = mesh_lap_int(lap,bn,f2);
%%%kann weg %vint_orig = vint;

%% plot time field
% Es wird im Bereich der Stickzeit über den ganzen Körper, in gleiche
% Abschnitte mit 0.05 gebildet


% visualize knitting time with time contours
levelT2 = min(vint):0.05:max(vint);
level_time = cell(length(levelT2),1);
figure;
% optionally: compute skeleton if needed (not needed at the moment)
skel = zeros(length(levelT2),3);
for i=1:length(levelT2)
    % compute (unordered) level sets
    [level_time{i},triidx] = compute_level_set4_KP(TR_start,edges_unique,levelT2(i),vint,[]);
    res_temp = draw_contours(level_time{i},triidx,'b',0,[]);   
    %skel(i,:) = mean(level_time{i},2)';
    level_time{i} = pc_rotate(level_time{i});
end
%plot the graph to see the skeleton 
%figure;
%plot3(skel(:,1),skel(:,2),skel(:,3))
%% cut through seam
% separate point at seam
if ~isempty(seam)
    % Wenn eine Naht da ist, werden die Punkte verdoppelt und das Mesh an
    % der Stelle der Naht durchgeschnitten
    [TR_start,edges_unique,seamidx2,lb,ub,vint,seamtm,mididx2] = ...
        cut_lines(TR_start,seamidx,edges_unique,lb,ub,vint,mididx);
    dmid = dist_from_mid(TR_start,mididx2);
else
    seamidx2 = [];
    seamtm = [];
    dmid = dist_from_mid(TR_start,mididx);
end
% dmid ist der Abstand zur Mittelline über das Mesh

%% plot distance from seam (or mid)
% compute and plot distance field, starting from the middle line
levelT2 = min(dmid):w:max(dmid);
level_dmid = cell(length(levelT2),1);
scalar_field = {vint};
level_field = cell(length(levelT2),length(scalar_field));
all_res = cell(length(levelT2),1);
lt = cell(length(levelT2),1);
ltime = cell(length(levelT2),1);
for i=1:length(levelT2)
    [level_dmid{i},triidx,level_field(i,:)] = compute_level_set4_KP(TR_start,edges_unique,levelT2(i),dmid,scalar_field);    
    level_time = level_field{i,1};
    [~,all_res{i},lt{i},ltime{i}] = draw_contours(level_dmid{i},triidx,'r',2*h,level_time);
end

% reorder wales so that they occur in the correct order
all_idx = reorder_wales(all_res,lt,TR_start,w);
% build knitting courses
[all_cont,onedge,wvec, all_time, all_pts] = connect_nodes4(all_idx,all_res,w,ltime);
idx_first = sort(find(all_time<-1+1e-3));
idx_last = sort(find(all_time>1-1e-3));
% create the initial 2D structure
jac = create_initial_jac(all_cont,onedge,wvec,idx_first, idx_last);

% correct the initial jacquard, partly repeated in postprocess2D
jac = stitch_placement(jac);

[jac2,gleft,gright] = postprocess2D(jac,wvec);

jac_patch2 = cell(length(jac2),1);

for k=1:length(jac2)
    jaccolors = ones(length(jac2{k}),1);
     jac_patch2{k} = jaccolors;
     %scatter(jac2{k},k*ones(length(jac2{k}),1),5,jaccolors,'o','filled')
     %hold on; grid on
     %drawnow
end

%% generate ascii and a bitmap
[rgbs,trueops,ops] = importCAjson('ISPO2019_Farbpalette.scpx','knit_encoding_s2k.json');
jac_str = generate_ascii_json(jac2,gleft,gright,[],[],[],[],jac_patch2);
jac_flat = cell(length(jac_str),1);
for i=1:length(jac_str)
    jac_flat{i} = strcat(jac_str{i,:}); 
end

jac_ascii = jac_flat;
jacim = ones(length(jac_ascii), length(jac_ascii{1}), 3); % change here for background

colorsmat = zeros(length(trueops),3);
for i = 1:length(trueops)
    colorsmat(i,:) = reshape(rgbs(i,:,:),1,3); 
end

colorsmat = colorsmat./255;

for k = 1:length(jac_ascii)
    for k2 = 1:length(trueops)
        jacim(k,((jac_ascii{k})==trueops(k2)),1) = colorsmat(k2,1);
        jacim(k,((jac_ascii{k})==trueops(k2)),2) = colorsmat(k2,2); 
        jacim(k,((jac_ascii{k})==trueops(k2)),3) = colorsmat(k2,3);
    end 
end

figure 
image(jacim)
title('Final Color-Coded Pattern')
axis equal

%% save
jacpic = strcat([filename(1:end-4),'_',date,'_',KNIT_STR{KNIT},'.bmp']);
jacjac = strcat([filename(1:end-4),'_',date,'_',KNIT_STR{KNIT},'.jac']);
fnpic = fullfile(pwd,strcat(['..',filesep,'Results']),jacpic);
fnjac = fullfile(pwd,strcat(['..',filesep,'Results']),jacjac);
imwrite(jacim, fnpic)
write_jacquard_file(jac_flat,fnjac);

if ~isdeployed
    clear_interpreter_paths
end
fclose('all');
return
