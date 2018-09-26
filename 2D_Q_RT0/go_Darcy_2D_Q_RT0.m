function go_Darcy_2D_Q_RT0
%
% for reference, the old naming is:
% nodes2coord (= xyz)
% nodes2dofs (= node_dofs)
% elems2nodes (= ien)
% ----------------------------------------------------------------------
% by Bedrich Sousedik, June 2016.

%tic

dim = 2;           % physical space dimension (only ndim = 2 supported)

% (mean) value of the coefficient
coeffs_mean = 1;

%[nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_my_own_2D_Q_P1;
%[nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_ted_test_2D_Q_P1;
[nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_ted_test2_2D_Q_P1;
% adjust: facewise data for RT0 element
% ToDo: elems2faces = get_edges(elems2nodes')
% etc.
%[elems2faces,edges2nodes] = get_edges3(elems2nodes');
[elems2faces,edges2nodes] = get_edges2(elems2nodes');
elems2faces = elems2faces';
nfaces = max(max(elems2faces));
elengdof = nfaces + nelem;
figure(1); show_mesh2(elems2nodes',nodes2coord'); title('mesh');%...

coeffs = rand(dim,nelem)*10;
faces2nodes = edges2nodes;
tic
[Kloc,signs] = setup_stiffness_matrix_Darcy_2D_Q_RT0(...
            coeffs,nelem,elems2nodes,elems2faces,faces2nodes,nodes2coord);
elems2nodes2 = [elems2faces; nfaces+(1:nelem)];
Kloc0 = mat2cell(Kloc,size(Kloc,1),size(Kloc,2),ones(size(Kloc,3),1));
K = assemble(Kloc0,elems2nodes2,1:elengdof,1:nelem,elengdof);
time(1) = toc;
%[Kloc,signs] = setup_stiffness_matrix_Darcy_2D_Q_RT0_V(...
             %coeffs,nelem,elems2nodes,elems2faces,faces2nodes,nodes2coord);
tic
         
[jac,detj] = getDeterminants(elems2nodes,nodes2coord,dim,nelem);
signs = ones(nelem,4);
%signs = signs_edges_Q(elems2nodes');
%mass_matrix_RT0(elems2faces,B_K,B_K_det,signs,coeffs)
[A,Aloc] = mass_matrix_RT0(elems2faces',jac,detj,signs,coeffs);
time(2) = toc;
%
%...

[B_Ktest,b_Ktest,B_K_dettest] = affine_transformations(nodes2coord',elems2nodes');

[B,Bloc] = pv_matrix_RT0(elems2faces',detj,signs);


% global matrix
K2 = [A B'; B sparse(size(B,1),size(B,1))];
% time(2) = toc;

[MPV] = mpv_matrix_RT0(elems2faces',jac,detj,signs,coeffs);
MPV0 = mat2cell(MPV,size(MPV,1),size(MPV,2),ones(size(MPV,3),1));
MPVA = assemble(MPV0,elems2nodes2,1:elengdof,1:nelem,elengdof);
err = zeros(nelem,1);
for i = 1:size(Kloc,3)
   err(i) = norm(Kloc(:,:,i)-MPV(:,:,i));
end
disp(sum(err))
err2 = zeros(nelem,1);
for i = 1:size(Kloc,3)
   err2(i) = norm(Kloc(1:4,1:4,i)-Aloc(:,:,i));
end
disp(sum(err2))

% err2 = zeros(nelem,1);
% for i = 1:size(Kloc,4)
%    err2(i) = norm(Kloc(1:4,1:4,i)-Aloc1(:,:,i));
% end
% disp(sum(err2))
% [maxe,I]=max(err);

norm(full(K-K2))

[A1, B1] = mass_pv_matrix_RT0(elems2faces',jac,detj,signs,coeffs);

% %y-indexes
nbasis = 4;
Y = reshape(repmat(elems2nodes',nbasis,1),nbasis,nbasis,nelem);    
% % x-indexes (exchanges row and cols and preserves the third/elems dim.) 
X = permute(Y,[2 1 3]); 
% % mass matrix
MASS = sparse(X(:),Y(:),MASS(:)); 
% 
% % export
% STIFF_loc = STIFF;
% 
% %commented for now
% Y = elems';% Y = reshape(repmat(elems',nbasis,1),nbasis,nbasis,nelems);    % y-indexes
Y = Y(:);                                    % x-indexes
% % new y indeces % ToDo
X = ones(nbasis,nelems)*diag(1:nelems);
STIFF = sparse(X(:),Y(:),STIFF(:));                           % stiffness matrix

% err = zeros(nelem,1);
% for i = 1:size(Kloc,4)
%    err(i) = norm(Kloc(1:4,1:4,i)-Aloc(:,:,i));
% end
% disp(sum(err))
return %end of function



% geometry of the physical domain
% ndim = 2;           % physical space dimension (only ndim = 2 supported)
dx = 1;              % dimension (size) of the domain in x-direction
dy = 1;              % dimension (size) of the domain in y-direction
nx = 2; ny = nx;     % number of elements in one direction (on a side)

% (mean) value of the coefficient
a_mean = 1;

hx = dx/nx; hy = dy/ny; % mesh step 
nelem = nx*ny;      % number of finite elements
% volD = dx*dy;       % volume(=area) of the physical domain

[elems2nodes,nodes2coord] = create_mesh_Q_2D_P1(hx,hy,nx,ny);

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

return

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

return % end of function