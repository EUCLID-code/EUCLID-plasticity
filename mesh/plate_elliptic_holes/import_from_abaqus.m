
%% Extract Nodes and Elements from Abaqus Data File
% using `Import data from file` Matlab functionality

% array `node` contains the positions of all nodes in the physical space
% each row corresponds to one node
node = JobRefinex1;
% array `element` contains information about the element connectivity
% each row corresponds to one element
element = JobRefinex1; 

save('mesh\plate_elliptic_holes\node_plate_elliptic_holes_x1','node')
save('mesh\plate_elliptic_holes\element_plate_elliptic_holes_x1','element')








