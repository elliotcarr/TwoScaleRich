function mesh = micro_mesh(mesh_name, visualise)
% ------------------------------------------------------------------------------
% micro_mesh.m
%
% Written by Dr Elliot Carr (2013-2014)
% Ecole Centrale Paris and Queensland University of Technology
% 
% This code is part of TwoScalRich.
%
% This MATLAB file generates generates properties such as CV areas, 
% edge lengths, etc. for solving the periodic cell problem 
% (see effective_conductivity.m, Ffunc1.c, Ffunc2.c and compute_Keff.c)
%
% ------------------------------------------------------------------------------

load(mesh_name)
a = m_width;
b = m_height;
no_elements = size(elements,1);
no_nodes = size(nodes,1);

% North, south, east and west boundary nodes
north_nodes = find((abs(nodes(:,2) - b) < 1e-4) == 1);                                       
east_nodes  = find((abs(nodes(:,1) - a) < 1e-4) == 1);
south_nodes = find(nodes(:,2) == 0);
west_nodes  = find(nodes(:,1) == 0);

variable        = ones(no_nodes,1);
variable_number = zeros(no_nodes,1);                                                                                %

% Only nodes not located on north and east boundaries are treated as unknowns
% variable(i) = 1 if node i is treated as unknown
% variable(i) = 0 if node i is treated as known (given by periodic BC)
no_variables = 0;
for i = 1:no_nodes
    if ismember(i,[north_nodes; east_nodes]) 
        % Known from periodic boundary condition
        variable(i) = 0;
    else
        % Variables are nodes not situated on the north or east boundaries
        no_variables = no_variables + 1; 
        variable_number(i) = no_variables;
    end
end

% Make sure number of nodes on opposing boundaries is equal
if length(north_nodes) ~= length(south_nodes)
    error('Number of nodes on north and south boundaries differ');
elseif length(east_nodes) ~= length(west_nodes)
    error('Number of nodes on east and west boundaries differ');
end

% Corner nodes
north_east_node = intersect(north_nodes,east_nodes);
north_west_node = intersect(north_nodes,west_nodes);
south_west_node = intersect(south_nodes,west_nodes);
south_east_node = intersect(south_nodes,east_nodes);

% Remove corner nodes from north, south, east and west boundary nodes
north_nodes = setdiff(north_nodes,[north_east_node; north_west_node]);
east_nodes  = setdiff(east_nodes, [north_east_node; south_east_node]);
south_nodes = setdiff(south_nodes,[south_east_node; south_west_node]);
west_nodes  = setdiff(west_nodes, [north_west_node; south_west_node]);

%-------------------------------------------------------------------------------
% Find matching nodes on opposing boundaries
partner_north = zeros(no_nodes,1);
partner_east  = zeros(no_nodes,1);

% Find partner node on south boundary matching north boundary node k
for k = 1:length(north_nodes)                                               
    found = false;
    x = nodes(north_nodes(k),1);
    for j = 1:no_nodes                                                      
        if (abs(nodes(j,1) - x) < 1e-6) && ismember(j,south_nodes)
            partner = j;
            nodes(partner,1) = x; % Ensure matching x-coordinates                                  
            found = true;
            break;
        end
    end
    if found
        partner_north(north_nodes(k)) = partner;
        variable(north_nodes(k)) = 0;
        variable_number(north_nodes(k)) = variable_number(partner);
    end
end

% Find partner node on west boundary matching east boundary node k
for k = 1:length(east_nodes)                                                
    found = false;                                                          
    y = nodes(east_nodes(k),2);
    for j = 1:no_nodes                                                      
        if (abs(nodes(j,2) - y) < 1e-6) && ismember(j,west_nodes)
            partner = j;
            nodes(partner,2) = y; % Ensure matching y-coordinates                                          
            found = true;
            break;
        end
    end
    if found
        partner_east(east_nodes(k)) = partner;
        variable(east_nodes(k)) = 0;
        variable_number(east_nodes(k)) = variable_number(partner);
    end
end

% All corner nodes qual to vale at south_west_node
if ~isempty(north_east_node)
    variable(north_east_node)        = 0;
    variable_number(north_east_node) = south_west_node;
end
if ~isempty(north_west_node)
    variable(north_west_node)        = 0;
    variable_number(north_west_node) = south_west_node;
end
if ~isempty(south_east_node)
    variable(south_east_node)        = 0;
    variable_number(south_east_node) = south_west_node;
end
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Generate mesh properties such as CV_areas, edge lengths, etc.

SCV_area     = zeros(no_nodes,3);
CV_area      = zeros(no_nodes,1);
normals      = zeros(no_elements,6);
shape_funcs  = zeros(no_elements,9);
element_area = zeros(no_elements,1);

areaA = 0;

for i = 1:no_elements
    vert = elements(i,:); %rectangle vertices
    xv   = nodes(vert,1); %x-coords
    yv   = nodes(vert,2); %y-coords
    
    % Element centroid
    centroid = [(xv(1) + xv(2) + xv(3)) / 3,(yv(1) + yv(2) + yv(3))/3];
    
    % Element edge midpoint coordinates
    midpts1   = [(xv(2) + xv(1)) / 2, (yv(2) + yv(1))/2];
    midpts2   = [(xv(2) + xv(3)) / 2, (yv(2) + yv(3))/2];
    midpts3   = [(xv(3) + xv(1)) / 2, (yv(3) + yv(1))/2];
    
    % Vectors connecting centroid to nodes
    p1        = centroid - [xv(1),yv(1)];
    p2        = centroid - [xv(2),yv(2)];
    p3        = centroid - [xv(3),yv(3)];
    
    % Vectors connecting midpoints
    q1        = midpts1  - midpts3;
    q2        = midpts2  - midpts1;
    q3        = midpts3  - midpts2;
    
    % Sub-control volume areas
    SCV_area1 = 0.5*abs(p1(1)*q1(2)-p1(2)*q1(1));
    SCV_area2 = 0.5*abs(p2(1)*q2(2)-p2(2)*q2(1));
    SCV_area3 = 0.5*abs(p3(1)*q3(2)-p3(2)*q3(1));    
    SCV_area(i,1:3) = [SCV_area1, SCV_area2, SCV_area3];
    
    % Area of element
    element_area(i) = SCV_area1 + SCV_area2 + SCV_area3;
    
    % Control volume edges
    fv1 = centroid - midpts1;
    fv2 = centroid - midpts2;
    fv3 = centroid - midpts3;
    
    % Control volume edge normals
    a = sign(p1(1)*q1(2)-p1(2)*q1(1));
    normals(i,:)  = -a * [fv1(2), -fv1(1), fv2(2), -fv2(1), fv3(2), -fv3(1)];
    
    % Shape function information
    detA = xv(2)*yv(3)-yv(2)*xv(3) - xv(1)*yv(3) + xv(3)*yv(1) + xv(1)*yv(2) ...
        - xv(2)*yv(1);
    shape_funcs(i,1)  = (yv(2)-yv(3)) / detA;
    shape_funcs(i,2)  = (yv(3)-yv(1)) / detA;
    shape_funcs(i,3)  = (yv(1)-yv(2)) / detA;
    shape_funcs(i,4)  = (xv(3)-xv(2)) / detA;
    shape_funcs(i,5)  = (xv(1)-xv(3)) / detA;
    shape_funcs(i,6)  = (xv(2)-xv(1)) / detA;
    shape_funcs(i,7)  = (xv(2)*yv(3) - yv(2)*xv(3)) / detA;
    shape_funcs(i,8)  = (yv(1)*xv(3) - xv(1)*yv(3)) / detA;
    shape_funcs(i,9)  = (xv(1)*yv(2) - yv(1)*xv(2)) / detA;
    
    % Area of sub-domain A
    if soil(i) == 1
        areaA = areaA + element_area(i);
    end
    
    % Control volume area
    CV_area(vert) = CV_area(vert) + SCV_area(i,1:3)';
    
end

cell_area = m_width * m_height;
epsilonA  = areaA / cell_area;
%-------------------------------------------------------------------------------

% Plot mesh
if visualise
    figure;
    colormap([0.8*[ones(32,1),ones(32,1),ones(32,1)]; ...
        0.6*[ones(32,1),ones(32,1),ones(32,1)]])
    p1 = patch('Faces',elements,'Vertices',nodes(:,1:2),'FaceColor','flat',...
        'EdgeColor','k');
    set(p1,'FaceVertexCData',soil');
    caxis([1,2])
    hold on
    vars = find(variable == 1);
    non_vars = find(variable == 0);
    p2 = plot(nodes(vars,1),nodes(vars,2),'b.','MarkerSize',26);
    p3 = plot(nodes(non_vars,1),nodes(non_vars,2),'r.','MarkerSize',26);
    legend([p2,p3],'Unknowns','Known from periodic boundary condition',...
        'Location','NorthOutside')    
    axis([0,m_width,0,m_height]);
    set(gca,'DataAspectRatio',[1,1,1])
    set(gcf,'Color','white')
    axis equal
    hold off
    drawnow    
end

mesh.elements           = elements;
mesh.nodes              = nodes;
mesh.no_elements        = no_elements;
mesh.no_nodes           = no_nodes;
mesh.no_variables       = no_variables;
mesh.normals            = normals;
mesh.shape_funcs        = shape_funcs;
mesh.variable_number    = variable_number;
mesh.variable           = variable;
mesh.epsilonA           = epsilonA;
mesh.CV_area            = CV_area;
mesh.cell_area          = cell_area;
mesh.soil               = soil;
mesh.element_area       = element_area;

end