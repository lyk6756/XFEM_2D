% ********************************************************************
% *****                     2D XFEM                              *****
% *****                 Nguyen Vinh Phu                          *****
% *****                    July 2006                             *****
% ********************************************************************

% Simple code for solving fracture mechanics using the eXtended Finite
% Element Method, XFEM
% Crack is represented by level sets

% Solving problem: Finite plate with edge crack
% under remote traction along Y direction

% ---------------------------------------
clear all
clc
state = 0;
tic;
% ---------------------------------------

% ++++++++++++++++++++++
%   GLOBAL VARIABLES
% ++++++++++++++++++++++
% Use of global variables help the info sharing
% between functions easier. Also makes the list of
% function parameters shorter

global node element  

% +++++++++++++++++++++++++++++
%             INPUT
% +++++++++++++++++++++++++++++

% ---------------------------------------
% Dimension of the domain
% (it is simply a rectangular region D x L)
L = 10 ;
D = 10 ;

% Material properties

% Matrix
E1  = 1e3 ;
nu1 = 0.3 ;

% Inclusion  (inclusion:°üº¬Îï£¬¼ÐÔÓ)
E2  = 1   ;
nu2 = 0.3 ;

stressState='PLANE_STRAIN';

% Loading
sigmato = 1  ;

% Circular inclusion
r  = 1.0 ;
xc = L/2 ;
yc = D/2 ;

% ---------------------------------------

% +++++++++++++++++++++++++++++
%            MESHING
% +++++++++++++++++++++++++++++

% ---------------------------------------
disp([num2str(toc),'   MESH GENERATION'])
% Number of nodes along two directions
nnx = 20 ;
nny = 20 ;

% Four corner points
pt1 = [0 0] ;
pt2 = [L 0] ;
pt3 = [L D] ;
pt4 = [0 D] ;

% Uniform meshing with Q4 elements
% Data structures for nodes and elements
% Node = [x1 y1
%         x2 y2
%         ...
%         xn yn]
% element = [1 3   5 7
%            4 20 35 78
%            ...
%           ]
elemType = 'Q4' ;
[node,element] = meshRectangularRegion(...
    pt1, pt2, pt3, pt4, nnx,nny,elemType);

% compute number of nodes, of elements
numnode = size(node,1);
numelem = size(element,1);

% define essential boundaries
uln = nnx*(nny-1)+1;       % upper left node number
urn = nnx*nny;             % upper right node number
lrn = nnx;                 % lower right node number
lln = 1;                   % lower left node number
cln = nnx*(nny-1)/2+1;     % node number at (0,0)

topEdge  = [ uln:1:(urn-1); (uln+1):1:urn ]';
botEdge  = [ lln:1:(lrn-1); (lln+1):1:lrn ]';

% GET NODES ON DIRICHLET BOUNDARY AND ESSENTIAL BOUNDARY
botNodes   = unique(botEdge);
topNodes   = unique(topEdge);

dispNodes = [botNodes];
tracNodes = topNodes;

% ---------------------------------------

% +++++++++++++++++++++++++++++
%           PROCESSING
% +++++++++++++++++++++++++++++
% Procedure:
%  1. Level set intialization
%  2. Set up enriched nodes, elements cut by crack
%  3. Initialize stiffness matrix, force vector
%  4. Loop on elements and compute stiffness matrix
%     and assemble into the global ones
%  5. Forming the force vector
%  6. Imposing essential boundary conditions
%  7. Solution of linear equations

% -------------------------------------------------------------
% ***************************************
% Level set initialization for all nodes
% ***************************************

% Level sets (normal and tangent) are stored in matrix ls
% whose the first column is normal LS and the second column contains
% the tangent LS. The i-th row is for node I
disp([num2str(toc),'   LEVEL SET INITIALIZATION'])
for i = 1 : numnode
    x = node(i,1);
    y = node(i,2);
    d = sqrt((x-xc)^2+(y-yc)^2);
    ls(i) = d - r ;
end
% Choose enriched nodes...

% for one element, if max(phi)*min(phi) < 0

% Data structures for elements cut by discontinuity
% Array split_elem contains the number of elements which are

enrich_node = zeros(numnode,1);

count = 0;

for iel = 1 : numelem
    sctr = element(iel,:);
    phi  = ls(sctr);
    if ( max(phi)*min(phi) < 0 )
        count = count + 1 ; % ah, one split element
        split_elem(count) = iel;
        enrich_node(sctr)  = 1;
    end
end
split_nodes = find(enrich_node == 1);

% Plot mesh and enriched nodes to check
figure
hold on
plot_mesh(node,element,elemType,'b-');
theta = -pi:0.1:pi;
x = xc + r*cos(theta);
y = yc + r*sin(theta);
cir = plot(x,y,'k-');
set(cir,'LineWidth',2)
n1 = plot(node(split_nodes,1),node(split_nodes,2),'r*');
set(n1,'MarkerSize',12);
plot(node(dispNodes,1),node(dispNodes,2),'ks');
axis off

pause % watch the figure before continuing...

% -------------------------------------------------------------

% -------------------------------------------
% Initialize stiffness matrix, force vector

% Each split node is enriched by ONE function, H(x)
% Each tip node is enriched by FOUR functions, B_i(x)
% then, the total dofs is :
% total dof = numnode*nsd + numsplitnode*1*nsd + numstipnode*4*nsd
% here, two dimension, nsd = 2

total_unknown = numnode*2 + size(split_nodes,1)*1*2 ;

K = sparse(total_unknown,total_unknown);
f = zeros(total_unknown,1);
% -------------------------------------------

% ***********************************
%    Stiffness matrix computation
% ***********************************


% Due to the presence of additional dofs, the assembly is a little
% bit difficult than in FEM. We use fictitious nodes to handle these
% additional dofs. At a H(x) enriched node, we add one fantom node. These fictitious nodes
% are numbered from the total number of true nodes, ie, from numnode+1 ...

pos = zeros(numnode,1);
nsnode = 0 ;
for i = 1 : numnode
    if (enrich_node(i) == 1)
        pos(i) = (numnode + nsnode*1) + 1 ;
        nsnode = nsnode + 1 ;
    end
end

q = [];
disp([num2str(toc),'   STIFFNESS MATRIX COMPUTATION'])
% -----------------
% Loop on elements
% -----------------
for iel = 1 : numelem
    sctr = element(iel,:); % element connectivity
    nn   = length(sctr);   % number of nodes per element
    % -----------------------------------------------
    % Choose Gauss quadrature rules for elements
    if (ismember(iel,split_elem))     % split element
        order = 3 ;
        phi   = ls(sctr);
        [W,Q] = discontQ4quad(order,phi);
    else
        order = 2 ;
        [W,Q] = quadrature(order,'GAUSS',2);
    end
    % -----------------------------------------------

    % -----------------------------------------------
    % Transform these Gauss points to global coords
    % for plotting only
    for igp = 1 : size(W,1)
        gpnt = Q(igp,:);
        [N,dNdxi]=lagrange_basis('Q4',gpnt);
        Gpnt = N' * node(sctr,:); % global GP
        q = [q;Gpnt];
    end
    % -----------------------------------------------

    % Stiffness matrix Ke = B^T C B
    % B = [ Bfem | Bxfem ]

    % The nodal parameters are stored in the following manner:
    % u = [u1 v1 u2 v2 ... un vn | a1 b1 ... an bn]';
    % It means that the additional unknowns are placed at the end of vector
    % u. Hence, the true nodal displacement is simply :
    % u_x = u(1:2:2*numnode) and u_y = u(2:2:2*numnode)

    % There are two kinds of element:
    % 1. Non-enriched one (usual element)          : calculer comme d'habitude
    % 2. Enriched one (fully and partially as well): special traitement

    % Determine the position in the global matrix K
    sctrB = assembly(iel,enrich_node,pos);

    % ---------------------
    % Loop on Gauss points
    % ---------------------
    for kk = 1 : size(W,1)
        pt = Q(kk,:);                             % quadrature point
        % B matrix

        [N,dNdxi] = lagrange_basis(elemType,pt);  % element shape functions
        J0 = node(sctr,:)'*dNdxi;                 % element Jacobian matrix
        invJ0 = inv(J0);
        dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY
        Gpt = N' * node(sctr,:);                  % GP in global coord, used

        % Bfem is always computed
        Bfem = zeros(3,2*nn);
        Bfem(1,1:2:2*nn)  = dNdx(:,1)' ;
        Bfem(2,2:2:2*nn)  = dNdx(:,2)' ;
        Bfem(3,1:2:2*nn)  = dNdx(:,2)' ;
        Bfem(3,2:2:2*nn)  = dNdx(:,1)' ;

        % Switch between non-enriched and enriched elements
        if ( any(enrich_node(sctr)) == 0 ) % Non-enriched elements
            B = Bfem;
        else                               % Enriched elements
            Bxfem = [] ;
            % loop on nodes, check node is enriched ...
            for in = 1 : nn
                if ( enrich_node(sctr(in)) == 1)     % H(x) enriched node
                    % compute the enrichment function and derivatives
                    % at gauss point
                    x = Gpt(1);
                    y = Gpt(2);
                    d = sqrt((x-xc)^2+(y-yc)^2);
                    F    = abs(d - r);
                    dFdx = sign(F)*(x-xc)/d;
                    dFdy = sign(F)*(y-yc)/d;

                    % compute the enrichment function at node               
                    x = node(sctr(in),1);
                    y = node(sctr(in),2);
                    d = sqrt((x-xc)^2+(y-yc)^2);
                    FI = abs(d - r);
                           
                    % Bxfem at node "in"
                    aa = dNdx(in,1)*(F-FI) + N(in)*dFdx ;
                    bb = dNdx(in,2)*(F-FI) + N(in)*dFdy ;                    
                    BI_enr = [aa 0; 0 bb; bb aa];
                    
                    % Add to the total Bxfem
                    Bxfem = [Bxfem BI_enr];
                    clear BI_enr ;
                end
            end          % end of loop on nodes
            % B matrix
            B = [ Bfem Bxfem ];
            clear Bfem; clear Bxfem;
        end              % end of switch between enriched and non-enriched elements

        % Compliance matrix
        x = pt(1);
        y = pt(2);
        d = sqrt((x-xc)^2+(y-yc)^2);
        if d < r
            E  = E2 ;
            nu = nu2;
        else
            E  = E1 ;
            nu = nu1;
        end

        if ( strcmp(stressState,'PLANE_STRESS') )
            C=E/(1-nu^2)*[ 1   nu 0;
                nu  1  0 ;
                0   0  0.5*(1-nu) ];
        else
            C=E/(1+nu)/(1-2*nu)*[ 1-nu  nu  0;
                nu    1-nu 0;
                0     0  0.5-nu ];
        end

        % Stiffness matrix
        K(sctrB,sctrB) = K(sctrB,sctrB) + B'*C*B*W(kk)*det(J0);
    end                  % end of looping on GPs
end                      % end of looping on elements

% -------------------------------------
% Plot GPs for checking
figure
hold on
plot_mesh(node,element,elemType,'b-');
plot(q(:,1),q(:,2),'r*');
theta = -pi:0.1:pi;
x = xc + r*cos(theta);
y = yc + r*sin(theta);
cir = plot(x,y,'k-');
set(cir,'LineWidth',2)
axis off

% -------------------------------------

% ************************
%    NODAL FORCE VECTOR
% ************************
disp([num2str(toc),'   NODAL FORCE VECTOR COMPUTATION'])

% The top edge is applied a traction along Y direction
[W,Q]=quadrature(1,'GAUSS',1);
for e = 1:size(topEdge,1)
    sctr = topEdge(e,:);
    sctry = sctr.*2 ;

    for q=1:size(W,1)
        pt = Q(q,:);
        wt = W(q);
        N  = lagrange_basis('L2',pt);
        J0 = abs( node(sctr(2))-node(sctr(1)) )/2;
        f(sctry)=f(sctry)+N*sigmato*det(J0)*wt;
    end   % of quadrature loop
end       % of element loop
% -------------------------------------


% **********************************
%    ESSENTIAL BOUNDARY CONDITION
% **********************************
disp([num2str(toc),'   IMPOSING ESSENTIAL BOUNDARY CONDITION'])

bcwt = mean(diag(K)); % a measure of the average size of an element in K
% used to keep the conditioning of the K matrix

udofs = zeros(length(dispNodes),1);
vdofs = zeros(length(dispNodes),1);

vdofs = dispNodes.*2;
udofs = [1]; % for lower left corner node

f(udofs) = 0 ;
f(vdofs) = 0 ;

K(udofs,:) = 0;   % zero out the rows and columns of the K matrix
K(vdofs,:) = 0;
K(:,udofs) = 0;
K(:,vdofs) = 0;
K(udofs,udofs) = bcwt*speye(length(udofs)); % put ones*bcwt on the diagonal
K(vdofs,vdofs) = bcwt*speye(length(vdofs));


% **********************************
%       SOLUTION OF EQUATIONS
% **********************************
disp([num2str(toc),'   SOLUTION'])
u   = K\f;
u_x = u(1:2:2*numnode) ;
u_y = u(2:2:2*numnode) ;

% **********************************
%          POST PROCESSING
% **********************************
disp([num2str(toc),'   POST PROCESSING'])

disp([num2str(toc),'      Deformed configuration'])
% --------------------------------------------------
% Plot numerical deformed configuration

figure
hold on
fac = 550;
plot_mesh(node+fac*[u_x u_y],element,'Q4','k')
axis off
title(' Scaled deformed configuration ')

figure
hold on
fac = 550;
plot_field(node+fac*[u_x u_y],element,elemType,u_y)
colorbar
axis off
title(' Y displacement on deformed configuration ')
% --------------------------------------------------

% ---------------------------------------------
% Compute stress at nodes and plot
disp([num2str(toc),'      Stress computation'])

stressPoints=[-1 -1;1 -1;1 1;-1 1];

for e = 1 : numelem
    sctr = element(e,:);
    nn   = length(sctr);
    U     = element_disp(e,pos,enrich_node,u);
    for q=1:nn
        pt = stressPoints(q,:);
        % B matrix
        [N,dNdxi] = lagrange_basis(elemType,pt);  % element shape functions
        J0 = node(sctr,:)'*dNdxi;                 % element Jacobian matrix
        invJ0 = inv(J0);
        dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY
        Gpt = N' * node(sctr,:);                  % GP in global coord, used

        % Bfem is always computed
        Bfem = zeros(3,2*nn);
        Bfem(1,1:2:2*nn)  = dNdx(:,1)' ;
        Bfem(2,2:2:2*nn)  = dNdx(:,2)' ;
        Bfem(3,1:2:2*nn)  = dNdx(:,2)' ;
        Bfem(3,2:2:2*nn)  = dNdx(:,1)' ;

        % Switch between non-enriched and enriched elements
        if ( any(enrich_node(sctr)) == 0 ) % Non-enriched elements
            B = Bfem;
        else                               % Enriched elements
            Bxfem = [] ;
            % loop on nodes, check node is enriched ...
            for in = 1 : nn
                if ( enrich_node(sctr(in)) == 1)     % H(x) enriched node
                    % compute the enrichment function and derivatives
                    % at gauss point
                    x = Gpt(1);
                    y = Gpt(2);
                    d = sqrt((x-xc)^2+(y-yc)^2);
                    F    = abs(d - r);
                    dFdx = sign(F)*(x-xc)/d;
                    dFdy = sign(F)*(y-yc)/d;

                    % compute the enrichment function at node               
                    x = node(sctr(in),1);
                    y = node(sctr(in),2);
                    d = sqrt((x-xc)^2+(y-yc)^2);
                    FI = abs(d - r);
                           
                    % Bxfem at node "in"
                    aa = dNdx(in,1)*(F-FI) + N(in)*dFdx ;
                    bb = dNdx(in,2)*(F-FI) + N(in)*dFdy ;                    
                    BI_enr = [aa 0; 0 bb; bb aa];
                    
                    % Add to the total Bxfem
                    Bxfem = [Bxfem BI_enr];
                    clear BI_enr ;
                end
            end          % end of loop on nodes
            % B matrix
            B = [ Bfem Bxfem ];
            clear Bfem; clear Bxfem;
        end              % end of switch between enriched and non-enriched elements
        
        % Compliance matrix
        x = pt(1);
        y = pt(2);
        d = sqrt((x-xc)^2+(y-yc)^2);
        if d < r
            E  = E2 ;
            nu = nu2;
        else
            E  = E1 ;
            nu = nu1;
        end

        if ( strcmp(stressState,'PLANE_STRESS') )
            C=E/(1-nu^2)*[ 1   nu 0;
                nu  1  0 ;
                0   0  0.5*(1-nu) ];
        else
            C=E/(1+nu)/(1-2*nu)*[ 1-nu  nu  0;
                nu    1-nu 0;
                0     0  0.5-nu ];
        end

        strain = B*U;
        stress(e,q,:) = C*strain;
    end
end

stressComp=1;
figure
clf
plot_field(node+fac*[u_x u_y],element,elemType,stress(:,:,stressComp));
colorbar
title('Stress plot, \sigma_{xx}')

% ---------------------------------------------



        
          
        

% -------------------------------------------------------------------------

