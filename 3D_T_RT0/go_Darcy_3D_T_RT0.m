function go_Darcy_3D_T_RT0
% function go_Darcy_3D_T_RT0
% Darcy flow in 3D, tetrahedral mesh, RT0 finite element discretization.
%
% for reference, the old naming is:
% nodes2coord (= xyz)
% nodes2dofs (= node_dofs)
% elems2nodes (= ien)
% ----------------------------------------------------------------------
% by Bedrich Sousedik, June 2016.

dim = 3;           % physical space dimension (only ndim = 2 supported)

% (mean) value of the coefficient
coeffs_mean = 1;

[nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_test01_3D_T_P1;
[nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_test02_3D_T_P1;
[nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_my_own_3D_T_P1;
%[nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_Matlab_3D_T_P1;

% adjust: facewise data for RT0 element
[elems2faces,faces2nodes] = get_faces(elems2nodes');
elems2faces = elems2faces'; 
faces2nodes = faces2nodes';
nfaces = max(max(elems2faces));
ngdof = nfaces + nelem;

% test:
nodes2coord(1,:) = nodes2coord(1,:) * 1/4;
nodes2coord(3,:) = nodes2coord(3,:) * 3/5;

% plot mesh
figure(1); show_mesh(elems2nodes',nodes2coord'); title('mesh');
% ToDo (perhaps): pdeplot(model) % alternatively use PDE toolbox

% coefficient of the equation
% coeffs = coeffs_mean*ones(nelem,1);
coeffs = 1*ones(dim,nelem);%(1:nelem)';
% coeffs = [1 0 0 ;0 2 0; 0 0 3]*coeffs;
coeffs = rand(dim,nelem);

% --- setup stiffness matrix: ------------------------------------------
% 1. standard approach: 
% tic
[Kloc,signs] = setup_stiffness_matrix_Darcy_3D_T_RT0(...
             coeffs,nelem,elems2nodes,elems2faces,faces2nodes,nodes2coord);
Kloc0 = mat2cell(Kloc,size(Kloc,1),size(Kloc,2),ones(size(Kloc,3),1));
% 
% assemble global matrices
elems2nodes2 = [elems2faces; nfaces+(1:nelem)];
% K = assemble(Kloc,elems2nodes,nodes2dofs,1:nelem,ngdof);
K = assemble(Kloc0,elems2nodes2,1:ngdof,1:nelem,ngdof);
% time(1) = toc;


% 2. vectorized approach
tic 
% affine transformations (same 2d and 3d, but only for triangles or tetrahedrons)
[B_K,~,B_K_det] = affine_transformations(nodes2coord',elems2nodes');
% signs = signs_faces(nodes2coord',elems2faces',faces2nodes',B_K);

% mass matrix assembly for RT0 element 
[A,Aloc] = mass_matrix_RT0(elems2faces',B_K,B_K_det,signs',coeffs);  

% pressure-velocity matrix 
% based on:
% K_RT0 = stiffness_matrix_RT0(elems2faces,B_K_det,signs); 
[B,Bloc] = pv_matrix_RT0(elems2faces',B_K_det,signs');

% global matrix
K2 = [A B'; B sparse(size(B,1),size(B,1))];

% time(2) = toc;


err = zeros(nelem,1);
for i = 1:size(Kloc,3)
   err(i) = norm(Kloc(1:4,1:4,i)-Aloc(:,:,i));
end
disp(sum(err))


norm(full(K-K2))


return % end of function

% OLD: --------

% setup degrees of freedom
nodes2dofs = 1:(nx+1)*(ny+1); 
ngdof = max(max(nodes2dofs)); 

% coefficient of the equation
k = a_mean*ones(nelem,1);

% get boundary conditions
bc = get_bc_laplace_2D(nodes2dofs,nodes2coord);
% test: bc(2,:)=[1 2 3]
if isempty(bc), error('Vector of boundary conditions empty.'), end

% setup element matrices
[Kloc,Mloc,ngdof] = setup_element_matrices_laplace_2D(...
             k,nelem,ngdof,elems2nodes,nodes2coord);

Kloc = mat2cell(Kloc,size(Kloc,1),size(Kloc,2),ones(size(Kloc,3),1));
Mloc = mat2cell(Mloc,size(Mloc,1),size(Mloc,2),ones(size(Mloc,3),1));

% assemble global matrices 
K = assemble(Kloc,elems2nodes,nodes2dofs,1:nelem,ngdof);
M = assemble(Mloc,elems2nodes,nodes2dofs,1:nelem,ngdof);
f = ones(length(K),1);

% apply the boundary conditions
if ~isempty(bc)
%    %[K,f] = apply_bc(K,f,bcdof,bcval);
[K,f] = apply_bc(K,f,bc(1,:),bc(2,:),1);%iblock);
[M,~] = apply_bc(M,f,bc(1,:),bc(2,:),1);%iblock);

     % this might work for vibrations
%      activeDof=setdiff([1:ngdof]',[bc(1,:)]);
%      K = K(activeDof,activeDof);
%      if iblock == 1
%         M = M(activeDof,activeDof);
%         f = f(activeDof);
%      end
else
   error('Missing boundary conditions to prescribe. ')
end



% % export
% %fem.name='model_laplace';
% fem.ndim = 2; 
% fem.ngdof = ngdof;
% fem.K = K; % blocks of the global stochastic matrix
% fem.f = f; % blocks of the global stochastic right-hand side
% fem.ien = ien;   % node-element connectivity (0 if unused)
% fem.node_dofs = id; % dofs on node (0 if unused)
% fem.xyz = xyz; % nodal coordinates

% if exist('f','var') % rhs
%    fem.f         = f;
% else
%    fprintf('\nRHS vector not available, generating random one.\n') 
%    fem.f = rand(ngdof,1);   % do not satisfy b.c. yet
% end


% solve the problem


% plot the results
% ToDo