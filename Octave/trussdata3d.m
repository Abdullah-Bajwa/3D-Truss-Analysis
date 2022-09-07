% Input file for trusscode3D.m

%Truss properties are stored as 3 different matrices
%materials matrix includes Elastic Modulus as well as Crosss Sectional Area of members
%nXYZ includes the location of all nodes in cartesian coordinates
%eprop defines the individual members by specifying the material and area flag along with the end nodes of each member
%Boundary Conditions are stored in 2 matrices
%bc matrix specifies the support status of nodes
%p matrix specifies the loads

% material props, E, A
%each row corresponds to one combination of area and material
%same row can be associated to multiple members if they have the same material as well as cross sectional area
% [material_flag, elastic_modulus, cross_sectional_area]
materials = [
   15e6 1
   3e7 (pi*((0.7)^2))/4
];

%coordinates of nodes
%Every row corresponds to one node
%list all nodes, order is not relevant as long as eprops is specified properly
% [x, y, z]
nXYZ = [
   0 0 30
   0 0 -30
   0 -30 0
   40 0 0
   ];

% element topology & material
% every row corresponds to one member 
%first two entries of the row are the positions of the end nodes in nXYZ
%starting from 1
%Third Entry is the material and area flag from materials matrix
%again, starting from 1 for first element
% [element, node 1(i), node 2(j), material_flag]
eProp = [
   1 4 1
   2 4 1
   3 4 1
   ];

%every node corresponds to one node
% [node, ux, uy, uz]
% boundary conditions,set value of ux/uy/uz = ( fixed-1, free-0)
bc = [
   1 1 1 1
   2 1 1 1
   3 1 1 1
   4 0 0 0
];

%Every row corresponds to one node
% loads in X,Y and Z directions
% [node, Fx, Fy, Fz]
p = [
   4 0 -5000 0
];

