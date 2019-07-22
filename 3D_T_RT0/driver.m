function [time, ngdof] = driver (nx, ny, nz)
% This is a general driver to run the 3D Block RT0 code.
% nx:  number of x elements
% ny:  number of y elements
% nz:  number of z elements
% ----------------------------------------------------------------
% by Theodore Weinberg, November 2016

dim = 3;

% (mean) value of the coefficient
coeffs_mean = 1;

[nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_num_3D_T_P1(nx,ny,nz);
% adjust: facewise data for RT0 element
% ToDo: elems2faces = get_edges(elems2nodes)
% etc.
[elems2faces,faces2nodes] = get_faces(elems2nodes');
elems2faces = elems2faces'; 
faces2nodes = faces2nodes';
nfaces = max(max(elems2faces));
ngdof = nfaces + nelem;
coeffs = 1*ones(dim,nelem);
coeffs = rand(dim,nelem);
%...
tic;
[Kloc,signs] = setup_stiffness_matrix_Darcy_3D_T_RT0(...
             coeffs,nelem,elems2nodes,elems2faces,faces2nodes,nodes2coord);
%...
time(1) = toc;
tic;
Kloc0 = mat2cell(Kloc,size(Kloc,1),size(Kloc,2),ones(size(Kloc,3),1));
% 
% assemble global matrices
elems2faces2 = [elems2faces; nfaces+(1:nelem)];
% K = assemble(Kloc,elems2nodes,nodes2dofsassemble(Kloc0,elems2nodes2,1:ngdof,1:nelem,ngdof),1:nelem,ngdof);
K = assembler(Kloc,elems2faces2,ngdof);
time(2) = toc;

tic
% affine transformations (same 2d and 3d, but only for triangles or tetrahedrons)
[B_K,~,B_K_det] = affine_transformations(nodes2coord',elems2nodes');
% mass matrix assembly for RT0 element 
%[A,Aloc] = mass_matrix_RT0(elems2faces',B_K,B_K_det,signs',coeffs);  
%[B,Bloc] = pv_matrix_RT0(elems2faces',B_K_det,signs');
signs2 = signs_faces(nodes2coord',elems2faces',faces2nodes',B_K);
[MPV] = mpv_matrix_RT0(elems2faces',B_K,B_K_det,signs2,coeffs);
time(3) = toc;
tic;
MPV0 = mat2cell(MPV,size(MPV,1),size(MPV,2),ones(size(MPV,3),1));
MPVA = assemble(MPV0,elems2faces2,1:ngdof,1:nelem,ngdof);
time(4) = toc;
% timeDisp = ['Time difference: ', num2str(time(1)-time(2))];
% disp(timeDisp);
% err = zeros(nelem,1);
% for i = 1:size(Kloc,3)
%    err(i) = norm(Kloc(:,:,i)-MPV(:,:,i));
% end
% disp(sum(err))
% [maxe,I]=max(err);
% norm(full(K-MPVA))

%norm(full(K-K2))

return %end of function