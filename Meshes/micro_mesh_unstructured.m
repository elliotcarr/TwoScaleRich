% ------------------------------------------------------------------------------
% micro_mesh_unstructured.m
%
% Written by Dr Elliot Carr (2013-2014)
% Ecole Centrale Paris and Queensland University of Technology
% 
% This code is part of TwoScalRich.
%
% This MATLAB script generates the microscopic mesh.
%
% ------------------------------------------------------------------------------

commandwindow
clc
close all
clear all

addpath('..') % To use image_based_meshing.m

%-------------------------------------------------------------------------------
% User changeable parameters 

% Meshing algorithm (MeshAdapt, Delaunay or Frontal)
mesh_algorithm = 'Frontal';

% Micro-cell
image  = '../micro_cell.png';

% Cell dimensions
M_elements_x = 20;
M_elements_y = 20;
m_width = 1.0/M_elements_x;
m_height = 1.0/M_elements_y;

% Mesh refinement parameter
r = 0.07;

% Mesh names
mesh_nameA = ['microA_unstructured_',num2str(M_elements_x),'by',...
    num2str(M_elements_y)];
mesh_nameB = ['microB_unstructured_',num2str(M_elements_x),'by',...
    num2str(M_elements_y)];

% Mesh to calculate effective conductivity (see Keff/run.m)
mesh_nameKeff = 'mesh_Keff';

% 'A' for sub-domain A (connected part) only
% 'AB' for full cell
domain_Keff = 'A'; 

% Path to gmsh
gmsh_path = '/Volumes/gmsh-2.7.1-MacOSX/Gmsh.app/Contents/MacOS/gmsh ';

% Visualise mesh (true or false)
visualise = true;
%-------------------------------------------------------------------------------

switch mesh_algorithm
    case 'MeshAdapt'
        ma = 1;
    case 'Delaunay'
        ma = 5;
    case 'Frontal'
        ma = 6;
    otherwise
        error(['Mesh algorithm must be one of either ''MeshAdapt'',' ...
            '''Frontal'' or ''Delaunay''.']);
end

% Generate gmsh .geo file
a = 1;
b = 1;
domain = 'B';
fid = fopen('mesh.geo', 'w');
image_based_meshing(image,fid,r,a,b,domain,1,1);
fclose(fid);

% Call gmsh and mesh the geometry
fprintf('%% Meshing microscopic domain\n');
system([gmsh_path,'mesh.geo -2']);

% Read gmsh .msh file
fid = fopen('mesh.msh', 'r');
while 1
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end
    switch tline
        case '$MeshFormat'
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            MF=sscanf(tline,'%g %g %g')';
            tline = fgetl(fid);
            if ~ischar(tline) || ~strcmp(tline,'$EndMeshFormat'), break, end
        case '$Nodes'
            tline = fgetl(fid);
            no_nodes = sscanf(tline,'%g');
            
            nodes = zeros(no_nodes,2);
            
            for i=1:no_nodes
                tline = fgetl(fid);
                type = sscanf(tline,'%g %g %g %g');
                nodes(i,1:2) = type(2:3);
            end
            
            tline = fgetl(fid);
            if ~ischar(tline) || ~strcmp(tline,'$EndNodes')
                break
            end
            
        case '$Elements'
            tline = fgetl(fid);
            if ~ischar(tline)
                break
            end
            no_elements = 0;
            no_objects  = sscanf(tline,'%g'); elements = [];
            boundary_nodes = []; no_boundary_edges = 0;
            for i=1:no_objects
                tline = fgetl(fid);
                type=sscanf(tline,'%g %g %g %g %g %g %g %g')';
                if type(2)==2
                    no_elements = no_elements+1;
                    elements(no_elements,1:3) = type(6:8);
                end
                if type(4) == 99
                    no_boundary_edges = no_boundary_edges + 1;
                    boundary_nodes = union(boundary_nodes,type(6:7));
                    boundary_edges(no_boundary_edges,1:3) = [type(6:7),-1];
                end
            end
            tline = fgetl(fid);
            if ~ischar(tline) || ~strcmp(tline,'$EndElements')
                break
            end
        otherwise
            disp('Unknown type encountered...')
            break
    end
    
end

fclose(fid);

nodes = [nodes(:,1)*m_width,nodes(:,2)*m_height];

%-------------------------------------------------------------------------------
% Generate .msh file read in by code (see micro_mesh_properties.cpp)
fid_mesh = fopen([mesh_nameB,'.msh'], 'w');
fprintf(fid_mesh,'%2.16f\n',m_width);
fprintf(fid_mesh,'%2.16f\n',m_height);
fprintf(fid_mesh,'3\n');

fprintf(fid_mesh,'%i\n',no_nodes);
for i = 1:no_nodes
    fprintf(fid_mesh,'%2.16f %2.16f\n',nodes(i,1:2));
end
fprintf(fid_mesh,'%i\n',no_boundary_edges);
for i = 1:no_boundary_edges
    fprintf(fid_mesh,'%i %i %i\n',boundary_edges(i,1:3));
end
fprintf(fid_mesh,'%i\n',no_elements);
for i = 1:no_elements
    fprintf(fid_mesh,'%i %i %i\n',elements(i,1:3));
end
fclose(fid_mesh);
%-------------------------------------------------------------------------------

% Plot mesh
if visualise
    figure; 
    triplot(elements,nodes(:,1),nodes(:,2),'k');
    view(2)
    hold on    
    % Plot boundary edges
    for i = 1:no_boundary_edges
        line([nodes(boundary_edges(i,1),1),nodes(boundary_edges(i,2),1)],...
            [nodes(boundary_edges(i,1),2),nodes(boundary_edges(i,2),2)],...
            'Color','r','LineWidth',2);
    end
    plot(nodes(boundary_nodes,1),nodes(boundary_nodes,2),'b.','MarkerSize',12)    
    axis([0,m_width,0,m_height]);
    set(gca,'DataAspectRatio',[1,1,1])
    set(gcf,'Color','white')
    drawnow    
end

% Save to use for plotting
save(mesh_nameB, 'elements', 'nodes','m_width','m_height','boundary_nodes')

%-------------------------------------------------------------------------------
% Mesh over sub-domain A (see visualise_twoscale.m)

% Generate gmsh .geo file
fid = fopen('mesh.geo', 'w');
a = 1;
b = 1;
domain = 'A';
image_based_meshing(image,fid,r,a,b,domain,1,1);
fclose(fid);

% Call gmsh and mesh the geometry
fprintf('%% Meshing microscopic domain\n');
system([gmsh_path,'mesh.geo -2']);

% Read gmsh .msh file
fid = fopen('mesh.msh', 'r');
while 1
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end
    switch tline
        case '$MeshFormat'
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            MF=sscanf(tline,'%g %g %g')';
            tline = fgetl(fid);
            if ~ischar(tline) || ~strcmp(tline,'$EndMeshFormat'), break, end
        case '$Nodes'
            tline = fgetl(fid);
            no_nodes = sscanf(tline,'%g');
            
            nodes = zeros(no_nodes,2);
            boundary_nodes = zeros(no_nodes,1);
            
            for i=1:no_nodes
                tline = fgetl(fid);
                type = sscanf(tline,'%g %g %g %g');
                nodes(i,1:2) = type(2:3);
            end
            
            tline = fgetl(fid);
            if ~ischar(tline) || ~strcmp(tline,'$EndNodes')
                break
            end
            
        case '$Elements'
            tline = fgetl(fid);
            if ~ischar(tline)
                break
            end
            no_elements = 0;
            no_objects  = sscanf(tline,'%g'); elements = [];
            for i=1:no_objects
                tline = fgetl(fid);
                type=sscanf(tline,'%g %g %g %g %g %g %g %g')';
                if type(2)==2
                    no_elements=no_elements+1;
                    elements(no_elements,1:3) = type(6:8);
                end
            end
            tline = fgetl(fid);
            if ~ischar(tline) || ~strcmp(tline,'$EndElements')
                break
            end
        otherwise
            disp('Unknown type encountered...')
            break
    end
    
end

fclose(fid);
nodes = [nodes(:,1)*m_width,nodes(:,2)*m_height];

% Plot mesh
if visualise
    figure;
    triplot(elements,nodes(:,1),nodes(:,2),'k');
    view(2)      
    axis([0,m_width,0,m_height]);
    set(gca,'DataAspectRatio',[1,1,1])
    set(gcf,'Color','white')
    drawnow    
end

% Save mesh (see visualise_twoscale.m)
save(mesh_nameA,'elements','nodes','m_width','m_height')

%-------------------------------------------------------------------------------
% Mesh for effective conductivity (see Keff/run.m)

% Generate gmsh .geo file
fid = fopen('mesh.geo', 'w');
a = 1;
b = 1;
image_based_meshing(image,fid,r,a,b,domain_Keff,1,1);
fclose(fid);

% Call gmsh and mesh the geometry
fprintf('%% Meshing microscopic domain\n');
system([gmsh_path,'mesh.geo -2']);

% Read gmsh .msh file
fid = fopen('mesh.msh', 'r');
while 1
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end
    switch tline
        case '$MeshFormat'
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            MF=sscanf(tline,'%g %g %g')';
            tline = fgetl(fid);
            if ~ischar(tline) || ~strcmp(tline,'$EndMeshFormat'), break, end
        case '$Nodes'
            tline = fgetl(fid);
            no_nodes = sscanf(tline,'%g');
            
            nodes = zeros(no_nodes,2);
            boundary_nodes = zeros(no_nodes,1);
            
            for i=1:no_nodes
                tline = fgetl(fid);
                type = sscanf(tline,'%g %g %g %g');
                nodes(i,1:2) = type(2:3);
            end
            
            tline = fgetl(fid);
            if ~ischar(tline) || ~strcmp(tline,'$EndNodes')
                break
            end
            
        case '$Elements'
            tline = fgetl(fid);
            if ~ischar(tline)
                break
            end
            no_elements = 0;
            no_objects  = sscanf(tline,'%g'); elements = []; soil = [];
            for i=1:no_objects
                tline = fgetl(fid);
                type=sscanf(tline,'%g %g %g %g %g %g %g %g')';
                if type(2)==2
                    no_elements=no_elements+1;
                    elements(no_elements,1:3) = type(6:8);
                    if ismember(type(4),1)
                        soil(no_elements) = 1;
                    elseif ismember(type(4),2)
                        soil(no_elements) = 2;
                    end
                end
            end
            tline = fgetl(fid);
            if ~ischar(tline) || ~strcmp(tline,'$EndElements')
                break
            end
        otherwise
            disp('Unknown type encountered...')
            break
    end
    
end

fclose(fid);
nodes = [nodes(:,1)*m_width,nodes(:,2)*m_height];

% Plot mesh
if visualise
    figure;
    colormap([0.8*[ones(32,1),ones(32,1),ones(32,1)]; ...
        0.6*[ones(32,1),ones(32,1),ones(32,1)]])
    p1 = patch('Faces',elements,'Vertices',nodes(:,1:2),'FaceColor','flat',...
        'EdgeColor','k');
    set(p1,'FaceVertexCData',soil');
    caxis([1,2])
    view(2) 
    hold on
    triplot(elements,nodes(:,1),nodes(:,2),'k');
    hold off         
    axis([0,m_width,0,m_height]);
    set(gca,'DataAspectRatio',[1,1,1])
    set(gcf,'Color','white')
    drawnow    
end

% Save mesh
save('../Keff/mesh_Keff','elements','nodes','m_width','m_height','soil')

% Remove unused files
delete('mesh.geo','mesh.msh')