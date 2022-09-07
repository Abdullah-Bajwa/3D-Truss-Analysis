% graphics_toolkit ("gnuplot");
graphics_toolkit ("fltk");

function [Uxyz, strain, stress, axialforce] = truss3D(inputfile)
% Purpose:
%  Solves the plane truss problem
%	 with data as inputfile
% Input:   inputfile
% Output:
%		nodal displacements
%		element strains
%		element stresses
%		axial force in each element

% read in the input file
run(inputfile);

% the basic parameters:
nElems = size(eProp, 1);
nNodes = size(nXYZ,   1);
sDof   = 3 * nNodes;	  % System Dof

% initialization of all required matrices
ke = zeros(3);	             % element stiffness matrix (3x3)
K  = zeros(sDof);	           % global  stiffness matrix
u  = F = zeros(sDof,1);      % global displacement and load vectors
axialforce = stress = strain = zeros(nElems,1);   % result vectors
error = zeros(nNodes, 3);

% computation of elements' stiffness matrices
for e = 1:nElems
    i = eProp(e,1); j = eProp(e,2);	    % node i,j for element e
    m = eProp(e,3);						    % material for element e
   xi = nXYZ(i,1);  yi = nXYZ(i,2); zi = nXYZ(i,3);		 % (x,y) for node i
   xj = nXYZ(j,1);  yj = nXYZ(j,2); zj = nXYZ(j,3);		 % (x,y) for node j
   le = sqrt( (xj-xi)^2 + (yj-yi)^2 +(zj-zi)^2);  % length of element e
    cx = (xj - xi)/le; cy = (yj - yi)/le; cz = (zj-zi)/le; % direction cosines
   %disp func for troubleshooting
   %disp('directional cosines:');
   %disp(cx);
   %disp(cy);
   %disp(cz);
   ke = materials(m, 2)*materials(m, 1)/le * [
              cx*cx   cx*cy   cx*cz		             % element stiffness sub-matrix
              cx*cy   cy*cy   cy*cz
              cx*cz   cy*cz   cz*cz ];		       % in the global coordinates   
   %disp('stiffness matrix');
   %disp(ke);
   % assembly, by scattering [ke] into [K]
   si = 3*i - 2; sj = 3*j - 2;          % system i,j
	 K(si:si+2, si:si+2) = K(si:si+2, si:si+2) + ke;
	 K(si:si+2, sj:sj+2) = K(si:si+2, sj:sj+2) - ke;
	 K(sj:sj+2, si:si+2) = K(sj:sj+2, si:si+2) - ke;
	 K(sj:sj+2, sj:sj+2) = K(sj:sj+2, sj:sj+2) + ke;
end
   disp('stiffness matrix');
   disp(K);
   
% system load vector in GCS
for n=1:size(p,1)
   F(3*p(n,1)-2) = p(n,2);       % Fx at node n
   F(3*p(n,1)-1) = p(n,3);       % Fy at node n
   F(2*p(n,1)  ) = p(n,4);       % Fz at node n
end

% application of boundary conditions by 'row/column removal' method
bcdof=zeros(sDof,1);    % initialization
for n=1:size(bc,1)
   bcnode = 3*bc(n,1) -2;
   bcdof(bcnode:bcnode+2) = bc(n,2:4);
end

% Reduced forms for both [K] and {F}
Kr = K(~bcdof, ~bcdof);    % remove all rows & columns mentioned in bcdof
Fr = F(~bcdof);            % remove all entries corresponding to bcdof

%disp functions for troubleshooting
%disp(K);
%disp(Kr);
%disp(Fr);
% solution to [K]{u} = {F}, reduced equation
Ur = (Kr\Fr);

% Must add rows removed earlier to apply BCs
j = 1;
for i = 1:size(u,1)
    if bcdof(i)
        u(i) = 0;       % adding entries corresponding to removed rows
    else
        u(i) = Ur(j);   % u is now the restored nodal displacement vector
        j = j+1;        % next  entry
    end
end

% post processing
for e = 1:nElems
    i = eProp(e,1); j = eProp(e,2);         % node i,j of element e
    m = eProp(e,3);	                       % material of element e
    xi = nXYZ(i,1);  yi = nXYZ(i,2); zi = nXYZ(i,3);		 % (x,y) for node i
    xj = nXYZ(j,1);  yj = nXYZ(j,2); zj = nXYZ(j,3);		 % (x,y) for node j
    le = sqrt( (xj-xi)^2 + (yj-yi)^2 +(zj-zi)^2);        % length of element e
    cx = (xj - xi)/le; cy = (yj - yi)/le; cz = (zj-zi)/le; % direction cosines
  
  % calculating strains
   strain(e) = ( (cx * u(j*3-2) + cy * u(j*3-1) + cz * u(j*3)) - ...
               ( (cx * u(i*3-2) + cy * u(i*3-1) + cz * u(i*3) ) ) ) / le;
  
  % calculating stresses
   stress(e) = materials(m,1) * strain(e);      % stress = E * strain
 
  % calculating axial forces
   axialforce(e) = materials(m,2) * stress(e);  % F = stress * x-area
end

% Re-writing nodal displacements as three column matrix [Uix, Uiy]
for i=1:sDof/3
   Uxyz(i,1) = u(3*i-2);     %Uix
   Uxyz(i,2) = u(3*i-1);     % Uiy
   Uxyz(i,3) = u(3*i);     %Uiz
endfor

disp('The stresses are:');
disp(stress);
disp('The strains are:');
disp(strain);
disp('The axial forces are:');
disp(axialforce);
disp('The deflections are:');
disp(Uxyz);



% plotting the results of truss analysis

A = zeros(nNodes);	% adjacenecy matrix
for e = 1:nElems
    i = eProp(e,1); j = eProp(e,2);  % node i,j for element e
    A(i,j) = 1; A(j,i) = 1;
    
end
    %disp(A);
    %disp('done');
    %disp(nXYZ);


dXYZ = nXYZ;	% to store the deformed coordinates
    dXYZ = dXYZ + 10 * Uxyz;    % XY-coord of displaced Node i

hold on;
gplot3(A, nXYZ, "r:");
gplot3(A, dXYZ, "b-o");
hold off;
;


endfunction



