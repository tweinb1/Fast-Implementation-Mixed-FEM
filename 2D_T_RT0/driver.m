function [time, ngdof] = driver(nx, ny)
%
% This is a generic driver for triangles
% nx: number of elements in x direction, nx >= 1
% ny: number of elements in y direction, ny >= 1
% ----------------------------------------------------------------------
% by Bedrich Sousedik and Theodore Weinberg, June 2016.

%tic

dim = 2;           % physical space dimension (only ndim = 2 supported)

% (mean) value of the coefficient
coeffs_mean = 1;

%[nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_test01_2D_T_P1;
[nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_num_2D_T_P1(nx,ny);
%[nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_my_own_2D_T_P1;
%[nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_Matlab_2D_T_P1;

% adjust: facewise data for RT0 element
% note: regard edges in 2d as faces
[elems2faces,faces2nodes] = get_edges(elems2nodes'); 
elems2faces = elems2faces';
nfaces = max(max(elems2faces));
ngdof = nfaces + nelem;
%Come back to this
elengdof = nfaces+nelem;
% test
%nodes2coord = nodes2coord * 1/4;

% plot mesh
figure(1); show_mesh(elems2nodes',nodes2coord'); title('mesh');
% ToDo (perhaps): pdeplot(model) % alternatively use PDE toolbox

% coefficient of the equation
% coeffs = coeffs_mean*ones(nelem,1);
coeffs = 1*ones(dim,nelem);%(1:nelem)';
coeffs = [1 0;0 2]*coeffs;
% coeffs = [1 2; 3 4];
% coeffs = [1 1; 2 2];
coeffs = rand(dim,nelem);


% --- setup stiffness matrix: ------------------------------------------
% 1. standard approach: 
tic
[Kloc,signs,J] = setup_stiffness_matrix_Darcy_2D_T_RT0(...
             coeffs,nelem,elems2nodes,elems2faces,nodes2coord);
time(1) = toc;
tic
Kloc0 = mat2cell(Kloc,size(Kloc,1),size(Kloc,2),ones(size(Kloc,3),1));

% 
% assemble global matrices
elems2faces2 = [elems2faces; nfaces+(1:nelem)];

% nodes2dofs2 = 1:ngdof;
% K = assemble(Kloc,elems2nodes2,nodes2dofs2,1:nelem,ngdof);
K = assembler(Kloc,elems2faces2,elengdof);
time(2) = toc;
tic         
% 2. vectorized approach
% tic 
% affine transformations (same 2d and 3d, but only for triangles or tetrahedrons)
[B_K,~,B_K_det] = affine_transformations(nodes2coord',elems2nodes');
signs         = signs_edges(elems2nodes');
MPV = mpv_matrix_RT0(elems2faces',B_K,B_K_det,signs,coeffs); 
time(3) = toc;
MPV0 = mat2cell(MPV,size(MPV,1),size(MPV,2),ones(size(MPV,3),1));
MPVA = assemble(MPV0,elems2faces2,1:ngdof,1:nelem,ngdof);
time(4) = toc;

%Error testing
%timeDisp = ['Time difference: ', num2str(time(1)-time(2))];
%disp(timeDisp);
%err = zeros(nelem,1);
%for i = 1:size(Kloc,3)
%   err(i) = norm(Kloc(:,:,i)-MPV(:,:,i));
%end
%disp(sum(err))
%[maxe,I]=max(err);

%norm(full(K-MPVA))



%File Creation code

% flow 'out' of the domain (same sources in all elements) :
%bc = [1 3 11 12; 0 0 0 0]; 
%f(13:16) = [1 1 -1 -1];%- ones(4,1);

% f = zeros(elengdof,1);
% f(nfaces+1) = 1;
% f(nfaces+2) = 1;
% f(nfaces+3) = -1;
% f(elengdof) = -1;
% ftest = zeros(nfaces);
% [MPVtest2,ftest2] = apply_bc(MPVA,f,bc);
% fem.Kloc = MPV; % or A, B separately(?)
% fem.coeffs = coeffs;
% %fem.elems2faces2 = elems2faces2;
% fem.nf = nfaces;
% %fem.elengdof = elengdof;
% fem.ngdof = elengdof;
% %fem.elems2faces = elems2faces;
% fem.elems2nodes = elems2nodes;
% %fem.faces2nodes = faces2nodes;
% fem.nodes2coord = nodes2coord;
% %fem.nodes2dofs = nodes2dofs;
% %fem.nodes2dofs = 1:elengdof;
% fem.signs = signs;
% fem.nelem = nelem;
% fem.bc = bc; % boundary conditions
% fem.f = f; % right-hand side
% fem.ndim = dim;
% fem.u = MPVtest2\ftest2;
% fem.faces2nodes = faces2nodes';
% fem.elems2faces = elems2faces;
% fem.elems2dofs = elems2faces2;
% 
% % save
% name = 'test_Teddy'; % this should be some identifier to worry about later 
% fprintf('\nSaving the structure fem into a file %s.mat ',name)
% fem.name = name;
% savefile = ['data_fem',name]; % this may change
% save(savefile,'fem','-V7.3')
% fprintf('... saved.\n\n')


return %end of function
