% ------------------------------------------------------------------------------
% macro_mesh_structured.m
%
% Written by Dr Elliot Carr (2013-2014)
% Ecole Centrale Paris and Queensland University of Technology
% 
% This code is part of TwoScalRich.
%
% This MATLAB script generates the macroscopic mesh.
%
% ------------------------------------------------------------------------------

commandwindow
clc
close all
clear all

%-------------------------------------------------------------------------------
% User changeable parameters

a = 1.0; % Domain width
b = 1.0; % Domain height
elements_x = 20; % Number of macro elements (x direction)
elements_y = 20; % Number of macro elements (y direction)

% Influx condition applied at north boundary for left < x < right
left = 0.0; 
right = 0.4; 

% Mesh name
mesh_name = ['macro_structured_',num2str(elements_x),'by',num2str(elements_y)];

% Path to gmsh
gmsh_path = '/Volumes/gmsh-2.7.1-MacOSX/Gmsh.app/Contents/MacOS/gmsh ';

% Visualise mesh (true or false)
visualise = true;

%-------------------------------------------------------------------------------

% Generate gmsh .geo file
fid = fopen('mesh.geo', 'w');
fprintf(fid,'a = %g;\n',a);
fprintf(fid,'b = %g;\n',b);
fprintf(fid,'left = %g;\n',left);
fprintf(fid,'right = %g;\n',right);
fprintf(fid,'r = %g;\n',a/(elements_x));
fprintf(fid,'Point(1) = {0,0,0,r};\n');
fprintf(fid,'Point(2) = {a,0,0,r};\n');
fprintf(fid,'Line(1) = {1,2};\n');
fprintf(fid,['Extrude {0,b,0} {Line{1}; Layers{',num2str(elements_y),...
    '}; Recombine; }\n']);
fclose(fid);

% Call gmsh and mesh the geometry
fprintf('%% Meshing macroscopic domain\n');
system([gmsh_path,'mesh.geo -2']);

% Read gmsh .msh file
fid_gmsh = fopen('mesh.msh', 'r');
while 1
    tline = fgetl(fid_gmsh);
    if ~ischar(tline)
        break
    end
    switch tline
        case '$MeshFormat'
            tline = fgetl(fid_gmsh);
            if ~ischar(tline), break, end
            MF=sscanf(tline,'%g %g %g')';
            tline = fgetl(fid_gmsh);
            if ~ischar(tline) || ~strcmp(tline,'$EndMeshFormat'), break, end
        case '$Nodes'
            tline = fgetl(fid_gmsh);
            no_nodes = sscanf(tline,'%g');
            
            nodes = zeros(no_nodes,2);
            boundary_nodes = zeros(no_nodes,5);
            
            for i = 1:no_nodes
                tline = fgetl(fid_gmsh);
                type = sscanf(tline,'%g %g %g %g');
                nodes(i,1:2) = type(2:3);
            end
            
            tline = fgetl(fid_gmsh);
            if ~ischar(tline) || ~strcmp(tline,'$EndNodes')
                break
            end
            
        case '$Elements'
            tline = fgetl(fid_gmsh);
            if ~ischar(tline)
                break
            end
            no_elements = 0;
            no_objects  = sscanf(tline,'%g'); elements = [];
            total = 0;
            for i=1:no_objects
                tline = fgetl(fid_gmsh);
                type=sscanf(tline,'%g %g %g %g %g %g %g %g')';
                if type(2) == 3
                    no_elements=no_elements+1;
                    elements(no_elements,1:4) = type(6:9);
                end
            end
            tline = fgetl(fid_gmsh);
            if ~ischar(tline) || ~strcmp(tline,'$EndElements')
                break
            end
        otherwise
            fprintf('Unknown type encountered... %s\n',tline)
            break;
    end
end

%-------------------------------------------------------------------------------
% Nodes located on influx, north, east, south and west boundaries
influx = [];
north  = [];
east   = [];
south  = [];
west   = [];

for i = 1:no_nodes
    if nodes(i,1) == 0
        west = [west,i];
    end
    if nodes(i,1) == a
        east = [east,i];
    end
    if nodes(i,2) == 0
        south = [south,i];
    end
    if nodes(i,2) == b
        if nodes(i,1) <= left+1.0e-6
            north = [north,i];
        end
        if nodes(i,1) >= left-1.0e-6 && nodes(i,1) <= right+1.0e-6
            influx = [influx,i];
        end
        if nodes(i,1) >= right-1.0e-6
            north = [north,i];
        end
    end
end
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Boundary edge types
% 1: influx, 2: south, 3: east, 4:west, 5: north
boundary_edges = [];
no_boundary_edges = 0;
for i = 1:no_elements    
    icv = elements(i,1:4);    
    for j = 1:4
        jnb = mod(j,4)+1;
        if ismember(icv(j),influx) && ismember(icv(jnb),influx)
            no_boundary_edges = no_boundary_edges + 1;
            boundary_edges(no_boundary_edges,1:3) = [icv(j), icv(jnb), -1];
        end
        if ismember(icv(j),south) && ismember(icv(jnb),south)
            no_boundary_edges = no_boundary_edges + 1;
            boundary_edges(no_boundary_edges,1:3) = [icv(j), icv(jnb), -2];
        end
        if ismember(icv(j),east) && ismember(icv(jnb),east)
            no_boundary_edges = no_boundary_edges + 1;
            boundary_edges(no_boundary_edges,1:3) = [icv(j), icv(jnb), -3];
        end
        if ismember(icv(j),west) && ismember(icv(jnb),west)
            no_boundary_edges = no_boundary_edges + 1;
            boundary_edges(no_boundary_edges,1:3) = [icv(j), icv(jnb), -4];
        end
        if ismember(icv(j),north) && ismember(icv(jnb),north)
            no_boundary_edges = no_boundary_edges + 1;
            boundary_edges(no_boundary_edges,1:3) = [icv(j), icv(jnb), -5];
        end
    end
end
fclose(fid_gmsh);
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Generate .msh file read in by code (see macro_mesh_properties.cpp)
fid_mesh = fopen(['',mesh_name,'.msh'], 'w');
fprintf(fid_mesh,'4\n');
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
    fprintf(fid_mesh,'%i %i %i %i\n',elements(i,1:4));
end
fclose(fid_mesh);
%-------------------------------------------------------------------------------

% Plot mesh
colors = {'r','g','b','c','m'};
if visualise
    trisurf(elements(:,1:4),nodes(:,1),nodes(:,2),zeros(no_nodes,1),...
        'FaceColor','None','EdgeColor','k');
    for i = 1:no_boundary_edges
        line([nodes(boundary_edges(i,1),1),nodes(boundary_edges(i,2),1)],...
            [nodes(boundary_edges(i,1),2),nodes(boundary_edges(i,2),2)],...
            'Color',colors{-boundary_edges(i,3)},'LineWidth',2);
    end
    view(2)
    grid off
    hold on
    plot(nodes(:,1),nodes(:,2),'k.','MarkerSize',12);
    bnd_nodes = find(boundary_nodes == 1);
    plot(nodes(bnd_nodes,1),nodes(bnd_nodes,2),'r.','MarkerSize',12);
    axis equal
    axis([0-b/10 a+b/10 0-b/10 b+b/10])
    hold off
    drawnow
end

% Element centroids (required in plot_solution.m)
centroid_elements = zeros(no_elements,2);
for k = 1:no_elements
    vert = elements(k,1:4);
    centroid_elements(k,1) = sum(nodes(vert,1))/4;
    centroid_elements(k,2) = sum(nodes(vert,2))/4;
    centroid = centroid_elements(k,:);
end

% Save mesh (see visualise_twoscale.m)
no_verts = 4;
save(mesh_name, 'elements', 'nodes','centroid_elements','no_verts')

% Remove unused files
delete('mesh.geo','mesh.msh')