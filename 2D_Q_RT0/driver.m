function [time, elengdof] = driver(nx, ny)
%
% This is a generic driver for quadrilaterals
% nx: number of elements in x direction, nx >= 1
% ny: number of elements in y direction, ny >= 1
% ----------------------------------------------------------------------
% by Bedrich Sousedik and Theodore Weinberg, June 2016.

%tic

dim = 2;           % physical space dimension (only ndim = 2 supported)

% (mean) value of the coefficient

%Create mesh and show our figure, setup input data, currently random coeffs
coeffs_mean = 1;
[nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_test_Q_P1(nx,ny);
[elems2faces,faces2nodes] = get_edges2(elems2nodes');
elems2faces = elems2faces';
nfaces = max(max(elems2faces));
elengdof = nfaces + nelem;
elems2faces2 = [elems2faces; nfaces+(1:nelem)];
figure(1); show_mesh2(elems2nodes',nodes2coord'); title('mesh');%...


%coeffs = rand(dim,nelem)*10;
coeffs = ones(dim,nelem);

%Non-vectorized
tic
[Kloc,signs] = setup_stiffness_matrix_Darcy_2D_Q_RT0(...
            coeffs,nelem,elems2nodes,elems2faces,faces2nodes,nodes2coord);
time(1) = toc;
tic

%Kloc0 = mat2cell(Kloc,size(Kloc,1),size(Kloc,2),ones(size(Kloc,3),1));
%K = assemble(Kloc0,elems2nodes2,1:elengdof,1:nelem,elengdof);
K = assembler(Kloc,elems2faces2,elengdof);
time(2) = toc;


%Vectorized
tic         
[jac,detj] = getDeterminants(elems2nodes,nodes2coord,dim,nelem);
%TED NOTE:  Modify signs by standard rules
%signs = ones(nelem,4);
signs = signs_edges_Q(elems2nodes');
%signs = [1 1 -1 -1; 1 1 -1 -1; 1 1 -1 -1; 1 1 -1 -1];
%test = mass_matrix_RT0(elems2faces,jac,detj,signs,coeffs)
[MPV] = mpv_matrix_RT0(elems2faces',jac,detj,signs,coeffs);
time(3) = toc;
tic
MPV0 = mat2cell(MPV,size(MPV,1),size(MPV,2),ones(size(MPV,3),1));
MPVA = assemble(MPV0,elems2faces2,1:elengdof,1:nelem,elengdof);
time(4) = toc;
timeDisp = ['Time difference: ', num2str(time(1)-time(2))];
disp(timeDisp);
err = zeros(nelem,1);
for i = 1:size(Kloc,3)
  err(i) = norm(Kloc(:,:,i)-MPV(:,:,i));
end
disp(sum(err))
[maxe,I]=max(err);
disp(time(3));
norm(full(K-MPVA))

Atest = MPVA(1:nfaces,1:nfaces);
Btest = MPVA(nfaces+1:elengdof,nfaces+1:elengdof);

%bc = [1 2 3 5 7 10 11 12; 0 0 0 0 0 0 0 0];
%bc = [2 7; 0 0]
% 'parabolic inflow' :
%bc = [2 7 1 3 11 12; 1/5 1/5 0 0 0 0];

% flow 'out' of the domain (same sources in all elements) :
bc = [1 3 11 12; 0 0 0 0]; 
%f(13:16) = [1 1 -1 -1];%- ones(4,1);

f = zeros(elengdof,1);
f(nfaces+1) = 1;
f(nfaces+2) = 1;
f(nfaces+3) = -1;
f(elengdof) = -1;
ftest = zeros(nfaces);
[MPVtest2,ftest2] = apply_bc(MPVA,f,bc);
fem.Kloc = MPV; % or A, B separately(?)
fem.coeffs = coeffs;
%fem.elems2faces2 = elems2faces2;
fem.nf = nfaces;
%fem.elengdof = elengdof;
fem.ngdof = elengdof;
%fem.elems2faces = elems2faces;
fem.elems2nodes = elems2nodes;
%fem.faces2nodes = faces2nodes;
fem.nodes2coord = nodes2coord;
%fem.nodes2dofs = nodes2dofs;
%fem.nodes2dofs = 1:elengdof;
fem.signs = signs;
fem.nelem = nelem;
fem.bc = bc; % boundary conditions
fem.f = f; % right-hand side
fem.ndim = dim;
fem.u = MPVtest2\ftest2;
fem.faces2nodes = faces2nodes';
fem.elems2faces = elems2faces;
fem.elems2dofs = elems2faces2;

%Special case code where original not invertible
% v1 =zeros(nfaces,1);
% v2 =  [1 1 1 1]';
% augment = zeros(17,17);
% augment(1:16,1:16)= MPVtest2;
% augment(17,1:12) = v1';
% augment(17,13:16) = v2';
% augment(1:12,17) = v1;
% augment(13:16,17) =v2;
% augment(13:16,17) =v2;
% rhs = [ftest2; 0];
% final = augment\rhs;
% fem.klocBC = MPVtest2;
% fem.augmentrhs = rhs;
% fem.final = final;
% fem.augment = augment;


% save
name = 'test_Teddy'; % this should be some identifier to worry about later 
fprintf('\nSaving the structure fem into a file %s.mat ',name)
fem.name = name;
savefile = ['data_fem',name]; % this may change
save(savefile,'fem','-V7.3')
fprintf('... saved.\n\n')




return %end of function
