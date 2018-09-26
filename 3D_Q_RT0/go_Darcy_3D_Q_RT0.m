function go_Darcy_3D_Q_RT0
%
%
%
% ----------------------------------------------------------------
% by Bedrich Sousedik, June 2016.

%tic

dim = 3;

% (mean) value of the coefficient
coeffs_mean = 1;

[nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_my_own_3D_Q_P1;

% adjust: facewise data for RT0 element
% ToDo: elems2faces = get_edges(elems2nodes)
% etc.
[elems2faces,faces2nodes] = get_faces2(elems2nodes');
elems2faces = elems2faces'; 
faces2nodes = faces2nodes';
nfaces = max(max(elems2faces));
ngdof = nfaces + nelem;
coeffs = 1*ones(dim,nelem);
coeffs = rand(dim,nelem);
%...
tic;
[Kloc,signs] = setup_stiffness_matrix_Darcy_3D_Q_RT0(...
             coeffs,nelem,elems2nodes,elems2faces,faces2nodes,nodes2coord);
%...
Kloc0 = mat2cell(Kloc,size(Kloc,1),size(Kloc,2),ones(size(Kloc,3),1));
% 
% assemble global matrices
elems2nodes2 = [elems2faces; nfaces+(1:nelem)];
% K = assemble(Kloc,elems2nodes,nodes2dofs,1:nelem,ngdof);
K = assemble(Kloc0,elems2nodes2,1:ngdof,1:nelem,ngdof);
time(1) = toc;

tic
% affine transformations (same 2d and 3d, but only for triangles or tetrahedrons)
[B_K,~,B_K_det] = affine_transformations(nodes2coord',elems2nodes');
% mass matrix assembly for RT0 element 
[A,Aloc] = mass_matrix_RT0(elems2faces',B_K,B_K_det,signs',coeffs);  
[B,Bloc] = pv_matrix_RT0(elems2faces',B_K_det,signs');

% global matrix
K2 = [A B'; B sparse(size(B,1),size(B,1))];
time(2) = toc;
err = zeros(nelem,1);
for i = 1:size(Kloc,3)
   err(i) = norm(Kloc(1:6,1:6,i)-Aloc(:,:,i));
end
disp(sum(err))
disp(time(1) - time(2))

%norm(full(K-K2))

return %end of function



% --- set parameters of the problem ------------------------------------
% 
% geometry of the physical domain
dim = 3;             % dimension of the physical space
dx = 8;              % dimension (size) of the domain in x-direction
dy = 2;              % dimension (size) of the domain in y-direction
dz = 1;
nx = 4; ny = 2; nz = 2;     % number of elements in one direction (on a side)

coeffs_mean = 0.01;         % coefficient 
% --- end of set parameters of the problem -------------------------------


% getmesh = 'create_mesh'; :
hx = dx/nx; hy = dy/ny; hz = dz/nz; % mesh step 
[nelem,elems2nodes,nodes2coord] = create_mesh_3D_Q_P1(hx,hy,hz,nx,ny,nz);
% abandoned here: getmesh = 'read_pmd_ns';

[nodes2dofs,ngdof] = create_nodes2dofs_3D_Q_P1(elems2nodes,nx,ny,nz);

% setup element matrices
coeffs = coeffs_mean*ones(nelem,1);   % coefficient of the problem
Kloc = setup_element_matrices_Laplace_3D_Q_P1(coeffs,nelem,elems2nodes,nodes2coord);

% assemble global matrices
fprintf('\nAssembling the global matrix ...')
% Kloc = mat2cell(Kloc,size(Kloc,1),size(Kloc,2),ones(size(Kloc,3),1));
K = assemble(Kloc,elems2nodes,nodes2dofs,1:nelem,ngdof);
f = ones(length(K),1); % generating arbitrary rhs
fprintf('\n ... done: ')

% OLD : ^^^^^^

% check if the two assemblies give the same
% norm(full(K-Kglo))

% right-hand side vector
%if ~exist('f','var')
   f = zeros(ngdof,1);    
%end  

fprintf('\nGetting boundary conditions ...'), tic
% get boundary conditions
if ~exist('bc','var')
   bc = get_bc(nodes2dofs,nodes2coord,dx,dy,dz);
end
fprintf('\n ... done: '), toc

% apply boundary conditions
fprintf('\nApplying boundary conditions ...'), tic
[A,f] = apply_bc(Kglo,f,bc);
fprintf('\n ... done: '), toc

return % end of function

if cavity   
   % fix one pressure dof (assume uniform cubes),
   % and pick dof on the node with coordinates (h,h,h) 
   xcoord = dx/nx;
   ycoord = dy/ny;
   zcoord = dz/nz;
   % DEBUG:   xcoord = 0; ycoord = 0; zcoord = 0;
   pdof = id(intersect(intersect(find(xyz(:,1)==xcoord),...
       find(xyz(:,2)==ycoord)),find(xyz(:,3)==zcoord)),dim+1); 
   pval = A(1,1);
   A(pdof,pdof) = pval;
end    

% (global) solution
fprintf('\nThere is %d degrees of freedom in the system, ',ngdof)
if ngdof <= 2.5e5
   fprintf(' solving the system ...'), tic
   u = A\f; 
   fprintf('\n ... done: \n'), toc
else
   fprintf(' not solving the system.')    
end


% settings for substructuring
sube = get_sube_uniform(dim,nsub,Hh);
% DEBUG:   sube = [1 2 1 2];

% name of the problem
if strcmp(getmesh,'create_mesh')
   nsub_side = nthroot(nsub,dim);
   nss = strcat(num2str(nsub_side),'x',num2str(nsub_side));
   for i=1:dim-2   %in 2D does nothing
      nss = strcat(nss,'x',num2str(nsub_side));
   end    
   name = strcat('stokes',nss,'_Hh',num2str(Hh));
else %read from pmd
   name = 'pmd'; 
end
% DEBUG:   name = 'stokes_test';
   
% export 
jmf.ndim = dim;
jmf.Kloc      = Kloc;
jmf.ien       = ien;
jmf.node_dofs = id;
jmf.xyz       = xyz;
if exist('A','var') % global stiffness matrix
   if ~issparse(A)  
      A = sparse(A);
   end
   jmf.A      = A; 
end
jmf.f         = f;
if exist('u','var') % solution
   jmf.u      = u; 
end
jmf.bc = [bcdof; bcval]; %bc;
jmf.pdof = pdof;
jmf.pval = pval; 
jmf.sube = sube;
jmf.name = name;

% save into a file
savefile = strcat('data_jmf/',name);
save(savefile,'jmf','-V7.3') 

%toc