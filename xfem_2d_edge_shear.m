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

% +++++++++++++++++++++++++++++
%             INPUT
% +++++++++++++++++++++++++++++

% ---------------------------------------
% Dimension of the domain
% (it is simply a rectangular region D x L)
L = 16 ;
D = 7 ;

% Material properties
E  = 3e7 ;
nu = 0.25 ;
stressState='PLANE_STRAIN';

% Crack properties
% Crack is represented as straight segments
% with points stored in matrix xCr
% Array xTip contains coordinates of the crack tip
% Array seg is the segment containg the tip
% It serves in computing alpha, the inclination angle
% of the crack w.r.t the X axis
a = 3.5 ;               % modeled crack length
xCr   = [0 0; a 0];
xTip  = [a 0];
seg   = xCr(2,:) - xCr(1,:);   % tip segment
alpha = atan2(seg(2),seg(1));  % inclination angle
QT    =[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];

% Loading
sigmato = 1  ;

% ---------------------------------------

% +++++++++++++++++++++++++++++
%            MESHING
% +++++++++++++++++++++++++++++

% ---------------------------------------
disp([num2str(toc),'   MESH GENERATION'])
% Number of nodes along two directions
nnx = 26 ;
nny = 46 ;

% Four corner points
pt1 = [ 0 -L/2];
pt2 = [ D -L/2];
pt3 = [ D L/2];
pt4 = [ 0 L/2];

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
x0  = xCr(1,1); y0 = xCr(1,2);
x1  = xCr(2,1); y1 = xCr(2,2);
t   = 1/norm(seg)*seg; %这里是将向量单位化
for i = 1 : numnode
    x = node(i,1);
    y = node(i,2);
    l   = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)) ;
    phi = (y0-y1)*x + (x1-x0)*y + (x0*y1-x1*y0);
    ls(i,1) = phi/l;            % normal LS 这里是单元节点到裂纹直线的距离
    ls(i,2) = ([x y]-xTip)*t';  % tangent LS 这里是节点和裂尖形成向量在裂纹直线上的投影
end

% Choose enriched nodes...

% for one element, if max(phi)*min(phi) < 0
% and max(psi) < 0, then it is a split element
% If max(phi)*min(phi) < 0 and max(psi)*min(psi) < 0, it is
% tip element

% Data structures for elements cut by crack
% Array split_elem contains the number of elements which are completely
% cut by crack. Similarly, array tip_elem stores number of tip element

enrich_node = zeros(numnode,1);

count1 = 0;
count2 = 0;
for iel = 1 : numelem
    sctr = element(iel,:);
    phi  = ls(sctr,1);
    psi  = ls(sctr,2);
    if ( max(phi)*min(phi) < 0 )
        if max(psi) < 0
            count1 = count1 + 1 ; % ah, one split element
            split_elem(count1) = iel;
            enrich_node(sctr)   = 1;
        elseif max(psi)*min(psi) < 0
            count2 = count2 + 1 ; % ah, one tip element
            tip_elem(count2) = iel;
            enrich_node(sctr)   = 2;
        end
    end
end
split_nodes = find(enrich_node == 1);
tip_nodes   = find(enrich_node == 2);

% Plot mesh and enriched nodes to check
figure
hold on
plot_mesh(node,element,elemType,'b-');
cr = plot(xCr(:,1),xCr(:,2),'r-');
set(cr,'LineWidth',3);
n1 = plot(node(split_nodes,1),node(split_nodes,2),'r*');
n2 = plot(node(tip_nodes,1),node(tip_nodes,2),'rs');
set(n1,'MarkerSize',15);
set(n2,'MarkerSize',15);
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

total_unknown = numnode*2 + size(split_nodes,1)*1*2 + ...
    size(tip_nodes,1)*4*2;

K = sparse(total_unknown,total_unknown);
f = zeros(total_unknown,1);
% -------------------------------------------

% ***********************************
%    Stiffness matrix computation
% ***********************************

% Compliance matrix C
if ( strcmp(stressState,'PLANE_STRESS') )
    C = E/(1-nu^2)*[ 1   nu 0;
        nu  1  0 ;
        0   0  0.5*(1-nu) ];
else
    C = E/(1+nu)/(1-2*nu)*[ 1-nu  nu  0;
        nu    1-nu 0;
        0     0  0.5-nu ];
end

% Due to the presence of additional dofs, the assembly is a little
% bit difficult than in FEM. We use fictitious nodes to handle these
% additional dofs. At a H(x) enriched node, we add one fantom node and
% at tip enriched node, four fantom nodes are added. These fictitious nodes
% are numbered from the total number of true nodes, ie, from numnode+1 ...

pos = zeros(numnode,1);
nsnode = 0 ;
ntnode = 0 ;
for i = 1 : numnode
    if (enrich_node(i) == 1)
        pos(i) = (numnode + nsnode*1 + ntnode*4) + 1 ;
        nsnode = nsnode + 1 ;
    elseif (enrich_node(i) == 2)
        pos(i) = (numnode + nsnode*1 + ntnode*4) + 1 ;
        ntnode = ntnode + 1 ;
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
    ke   = 0 ;             % elementary stiffness matrix
    % -----------------------------------------------
    % Choose Gauss quadrature rules for elements
    if (ismember(iel,split_elem))     % split element
        order = 2 ;
        phi   = ls(sctr,1);
        [W,Q] = discontQ4quad(order,phi);
    elseif (ismember(iel,tip_elem))   % tip element
        order = 7;
        phi   = ls(sctr,1);
        nodes = node(sctr,:);
        [W,Q] = disTipQ4quad(order,phi,nodes,xTip);
    elseif ( any(intersect(tip_nodes,sctr)) ~= 0)% having tip enriched nodes
        order = 4 ;
        [W,Q] = quadrature(order,'GAUSS',2);
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

    for k = 1 : nn
        sctrBfem(2*k-1) = 2*sctr(k)-1 ;
        sctrBfem(2*k)   = 2*sctr(k)   ;
    end

    if ( any(enrich_node(sctr)) == 0 ) % Non-enriched elements
        sctrB = sctrBfem ;
    else
        sn = size(find(enrich_node(sctr) == 1),1);
        tn = size(find(enrich_node(sctr) == 2),1);
        sctrBxfem = zeros(1,2*(sn*1+tn*4));
        cnt = 0 ;
        for k = 1 : nn
            if ( enrich_node(sctr(k)) == 1)
                cnt = cnt + 1 ;
                sctrBxfem(2*cnt - 1) = 2 * pos(sctr(k)) - 1;
                sctrBxfem(2*cnt    ) = 2 * pos(sctr(k))    ;
            elseif ( enrich_node(sctr(k)) == 2)
                cnt = cnt + 1 ;
                sctrBxfem(2*cnt - 1) = 2 * pos(sctr(k)) - 1;
                sctrBxfem(2*cnt    ) = 2 * pos(sctr(k))    ;

                cnt = cnt + 1 ;
                sctrBxfem(2*cnt - 1) = 2 * (pos(sctr(k))+1) - 1;
                sctrBxfem(2*cnt    ) = 2 * (pos(sctr(k))+1)    ;

                cnt = cnt + 1 ;
                sctrBxfem(2*cnt - 1) = 2 * (pos(sctr(k))+2) - 1;
                sctrBxfem(2*cnt    ) = 2 * (pos(sctr(k))+2)    ;

                cnt = cnt + 1 ;
                sctrBxfem(2*cnt - 1) = 2 * (pos(sctr(k))+3) - 1;
                sctrBxfem(2*cnt    ) = 2 * (pos(sctr(k))+3)    ;
            end
        end
        sctrB = [ sctrBfem sctrBxfem ];
    end
    % ---------------------
    % Loop on Gauss points
    % ---------------------
    for kk = 1 : size(W,1)
        pt = Q(kk,:);                             % quadrature point
        [N,dNdxi] = lagrange_basis(elemType,pt);  % element shape functions
        J0 = node(sctr,:)'*dNdxi;                 % element Jacobian matrix
        invJ0 = inv(J0);
        dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY
        Gpt = N' * node(sctr,:);                  % GP in global coord, used
        % to compute H(x),B(x)

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
                    % Enrichment function, H(x) at global Gauss point
                    dist = signed_distance(xCr,Gpt);
                    Hgp  = heaviside(dist);
                    % Enrichment function, H(x) at node "in"
                    dist = signed_distance(xCr,node(sctr(in),:));
                    Hi   = heaviside(dist);
                    % Bxfem at node "in"
                    BI_enr = [dNdx(in,1)*(Hgp - Hi) 0 ;
                        0 dNdx(in,2)*(Hgp - Hi) ;
                        dNdx(in,2)*(Hgp - Hi) dNdx(in,1)*(Hgp - Hi)];
                    % Add to the total Bxfem
                    Bxfem = [Bxfem BI_enr];
                    clear BI_enr ;
                elseif ( enrich_node(sctr(in)) == 2) % B(x) enriched node
                    % compute branch functions at Gauss point
                    xp    = QT*(Gpt-xTip)';           % local coordinates
                    r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
                    theta = atan2(xp(2),xp(1));
                    if ( theta > pi | theta < -pi)
                        disp (['something wrong with angle ',num2str(thet)]);
                    end
                    [Br,dBdx,dBdy] = branch(r,theta,alpha);

                    % compute branch functions at node "in"
                    xp    = QT*(node(sctr(in),:)-xTip)';
                    r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
                    theta = atan2(xp(2),xp(1));
                    if ( theta > pi | theta < -pi)
                        disp (['something wrong with angle ',num2str(thet)]);
                    end
                    [BrI] = branch_node(r,theta);

                    % composants of Benr matrix

                    aa = dNdx(in,1)*(Br(1)-BrI(1)) + N(in)*dBdx(1) ;
                    bb = dNdx(in,2)*(Br(1)-BrI(1)) + N(in)*dBdy(1) ;
                    B1_enr = [aa 0 ; 0 bb ; bb aa];

                    aa = dNdx(in,1)*(Br(2)-BrI(2)) + N(in)*dBdx(2) ;
                    bb = dNdx(in,2)*(Br(2)-BrI(2)) + N(in)*dBdy(2) ;
                    B2_enr = [aa 0 ; 0 bb ; bb aa];

                    aa = dNdx(in,1)*(Br(3)-BrI(3)) + N(in)*dBdx(3) ;
                    bb = dNdx(in,2)*(Br(3)-BrI(3)) + N(in)*dBdy(3) ;
                    B3_enr = [aa 0 ; 0 bb ; bb aa];

                    aa = dNdx(in,1)*(Br(4)-BrI(4)) + N(in)*dBdx(4) ;
                    bb = dNdx(in,2)*(Br(4)-BrI(4)) + N(in)*dBdy(4) ;
                    B4_enr = [aa 0 ; 0 bb ; bb aa];

                    BI_enr = [B1_enr B2_enr B3_enr B4_enr];
                    clear B1_enr; clear B2_enr; clear B3_enr; clear B4_enr;
                    Bxfem = [Bxfem BI_enr];
                    clear BI_enr ;
                end
            end          % end of loop on nodes
            % B matrix
            B = [ Bfem Bxfem ];
            clear Bfem; clear Bxfem;
        end              % end of switch between enriched and non-enriched elements
        % Stiffness matrix
        %ke = ke + B'*C*B*W(kk)*det(J0);
        K(sctrB,sctrB) = K(sctrB,sctrB) + B'*C*B*W(kk)*det(J0);
    end                  % end of looping on GPs
end                      % end of looping on elements

% -------------------------------------
% Plot GPs for checking
figure
hold on
plot_mesh(node,element,elemType,'b-');
plot(q(:,1),q(:,2),'r*');
cr = plot(xCr(:,1),xCr(:,2),'r-');
set(cr,'LineWidth',3);
axis off
clear q
% -------------------------------------

% ************************
%    NODAL FORCE VECTOR
% ************************
disp([num2str(toc),'   NODAL FORCE VECTOR COMPUTATION'])

% The top edge is applied a traction along Y direction
[W,Q]=quadrature(1,'GAUSS',1);
for e = 1:size(topEdge,1)
    sctr = topEdge(e,:);
    sctrx = sctr.*2-1 ;

    for q=1:size(W,1)
        pt = Q(q,:);
        wt = W(q);
        N  = lagrange_basis('L2',pt);
        J0 = abs( node(sctr(2))-node(sctr(1)) )/2;
        f(sctrx)=f(sctrx)+N*sigmato*det(J0)*wt;
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
u = K\f;
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
fac = 50;
plot_mesh(node+fac*[u_x u_y],element,elemType,'b-');
title(' Numerical deformed configuration ')
% --------------------------------------------------

disp([num2str(toc),'      Stress intensity factors computation'])
% -------------------------------------------------------------------------
% Compute the Stress Intensity Factors
% Using the Interaction integral method

% Steps :
% 1- detection of the elements on which we integrate
% 2- loop over these elements
% 3- loop over Gauss points
% 4- computation of stress, strain... in local coordinates !!!   ATTENTION
% 5- computation of the auxilliary fields: AuxStress and AuxEps and AuxGradDisp
% 6- computation of I1 and I2

% -------------------------------------
% calculation of the area of the tip element
x = node(element(tip_elem,:),:);
% Area = sum of areas of each sub-triangle
x0 = x(1,1);
y0 = x(1,2);

x1 = x(2,1);
y1 = x(2,2);

x2 = x(3,1);
y2 = x(3,2);

x3 = x(4,1);
y3 = x(4,2);

A1 = 0.5 * ((x0-x2)*(y1-y2) - (x1-x2)*(y0-y2)) ;
A2 = 0.5 * ((x0-x3)*(y2-y3) - (x2-x3)*(y0-y3)) ;
area = A1 + A2 ;

% J radius = fac * sqrt(area);
fac    = 3.0;
radius = fac * sqrt(area);
center = 1/4*sum(x)  ;

r=[];
% Distance from the center of tip element
for i = 1 : numnode
    sctr = node(i,:);
    rho  = sqrt((sctr(1)-center(1))^2+(sctr(2)-center(2))^2);
    r    = [r,rho];
end
test = r-radius;
test = test(element)';
test = max(test).*min(test);
Jdomain = element(find(test<=0),:);

test1 = r-radius;
test1 = test1(Jdomain)';
test1 = (test1<=0);
qnode = test1';

% plot
figure
hold on
% plot the circle
theta = -pi:0.1:pi;
xo = center(1) + radius*cos(theta) ;
yo = center(2) + radius*sin(theta) ;
plot(xo,yo,'k-');
plot_mesh(node,element,'Q4','b-')
plot_mesh(node,Jdomain,'Q4','r-')
plot(center(1),center(2),'b*')
cr = plot(xCr(:,1),xCr(:,2),'k-');
set(cr,'LineWidth',2);
% -------------------------------------

%---------------------------------------------
% Compute interaction integral

I1 = 0;
I2 = 0;
I  = [zeros(2,1)];
Ideb=[];
% ---------------------------
% Starting LOOP over ELEMENTS
%----------------------------
for e = 1 : size(Jdomain,1)
    sctr = Jdomain(e,:);
    nn   = length(sctr);
    % Choose Gauss quadrature rule
    if (ismember(iel,split_elem))     % split element
        order = 7 ; % 13 GPs for each sub-triangle
        phi   = ls(sctr,1);
        [W,Q] = discontQ4quad(order,phi);
    else
        order = 6 ; % 5x5 GPs
        [W,Q] = quadrature(order,'GAUSS',2);
    end
    % -----------------------------
    % start loop over Gauss points
    % -----------------------------
    for q = 1:size(W,1)
        pt = Q(q,:);
        wt = W(q);
        [N,dNdxi] = lagrange_basis(elemType,pt);
        J0    = node(sctr,:)'*dNdxi;
        invJ0 = inv(J0);
        dNdX  = dNdxi*invJ0;
        Gpt = N' * node(sctr,:);     % GP in global coord
        
        % +++++++++++++++++++++++++ 
        % Gradient of displacement
        % +++++++++++++++++++++++++ 
        % need to compute u,x u,y v,x v,y, stored in matrix H
    
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
                    % Enrichment function, H(x) at global Gauss point
                    dist = signed_distance(xCr,Gpt);
                    Hgp  = heaviside(dist);
                    % Enrichment function, H(x) at node "in"
                    dist = signed_distance(xCr,node(sctr(in),:));
                    Hi   = heaviside(dist);
                    % Bxfem at node "in"
                    BI_enr = [dNdx(in,1)*(Hgp - Hi) 0 ;
                        0 dNdx(in,2)*(Hgp - Hi) ;
                        dNdx(in,2)*(Hgp - Hi) dNdx(in,1)*(Hgp - Hi)];
                    % Add to the total Bxfem
                    Bxfem = [Bxfem BI_enr];
                    clear BI_enr ;
                elseif ( enrich_node(sctr(in)) == 2) % B(x) enriched node
                    % compute branch functions at Gauss point
                    xp    = QT*(Gpt-xTip)';           % local coordinates
                    r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
                    theta = atan2(xp(2),xp(1));
                    if ( theta > pi | theta < -pi)
                        disp (['something wrong with angle ',num2str(thet)]);
                    end
                    [Br,dBdx,dBdy] = branch(r,theta,alpha);

                    % compute branch functions at node "in"
                    xp    = QT*(node(sctr(in),:)-xTip)';
                    r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
                    theta = atan2(xp(2),xp(1));
                    if ( theta > pi | theta < -pi)
                        disp (['something wrong with angle ',num2str(thet)]);
                    end
                    [BrI] = branch_node(r,theta);

                    % composants of Benr matrix

                    aa = dNdx(in,1)*(Br(1)-BrI(1)) + N(in)*dBdx(1) ;
                    bb = dNdx(in,2)*(Br(1)-BrI(1)) + N(in)*dBdy(1) ;
                    B1_enr = [aa 0 ; 0 bb ; bb aa];

                    aa = dNdx(in,1)*(Br(2)-BrI(2)) + N(in)*dBdx(2) ;
                    bb = dNdx(in,2)*(Br(2)-BrI(2)) + N(in)*dBdy(2) ;
                    B2_enr = [aa 0 ; 0 bb ; bb aa];

                    aa = dNdx(in,1)*(Br(3)-BrI(3)) + N(in)*dBdx(3) ;
                    bb = dNdx(in,2)*(Br(3)-BrI(3)) + N(in)*dBdy(3) ;
                    B3_enr = [aa 0 ; 0 bb ; bb aa];

                    aa = dNdx(in,1)*(Br(4)-BrI(4)) + N(in)*dBdx(4) ;
                    bb = dNdx(in,2)*(Br(4)-BrI(4)) + N(in)*dBdy(4) ;
                    B4_enr = [aa 0 ; 0 bb ; bb aa];

                    BI_enr = [B1_enr B2_enr B3_enr B4_enr];
                    clear B1_enr; clear B2_enr; clear B3_enr; clear B4_enr;
                    Bxfem = [Bxfem BI_enr];
                    clear BI_enr ;
                end
            end          % end of loop on nodes
            % B matrix
            B = [ Bfem Bxfem ];
            clear Bfem; clear Bxfem;
        end
        leB = size(B,2);
        % nodal displacement of current element
        % taken from the total nodal parameters u
        idx = 0 ;
        U = zeros(2*nn,1);
        for in = 1 : nn
            idx = idx + 1;
            nodeI = sctr(in) ;
            U(2*idx-1) = u(2*nodeI-1);
            U(2*idx)   = u(2*nodeI  );
        end
        if ( any(enrich_node(sctr)) == 0 ) % Non-enriched elements
            U = U ;
        else                               % having enriched DOFs
            A = [];
            for in = 1 : nn
                nodeI = sctr(in) ;
                if (enrich_node(nodeI) == 1)     % H(x) enriched node
                    AA = [u(2*pos(nodeI)-1);u(2*pos(nodeI))];    
                    A  = [A;AA];
                elseif (enrich_node(nodeI) == 2) % B(x) enriched node
                    AA = [u(2*pos(nodeI)-1);
                          u(2*pos(nodeI));
                          u(2*(pos(nodeI)+1)-1);
                          u(2*(pos(nodeI)+1));
                          u(2*(pos(nodeI)+2)-1);                          
                          u(2*(pos(nodeI)+2));  
                          u(2*(pos(nodeI)+3)-1);                          
                          u(2*(pos(nodeI)+3));                                                    
                          ];  
                    A  = [A;AA];
                end                
            end
            U = [U;A];
        end
      
        H(1,1) = B(1,1:2:leB)*U(1:2:leB);    % u,x
        H(1,2) = B(2,2:2:leB)*U(1:2:leB);    % u,y
        H(2,1) = B(1,1:2:leB)*U(2:2:leB);    % v,x
        H(2,2) = B(2,2:2:leB)*U(2:2:leB);    % v,y
        
        % ++++++++++++++
        % Gradient of q
        % ++++++++++++++        
        q     = qnode(e,:);
        gradq = q*dNdX;

        % ++++++++++++++
        % Stress at GPs
        % ++++++++++++++ 
        
        epsilon = B*U ;
        sigma   = C*epsilon;
        
        % +++++++++++++++++++++++++++++++++++
        % Transformation to local coordinate
        % +++++++++++++++++++++++++++++++++++ 
        
        voit2ind    = [1 3;3 2];
        gradqloc    = QT*gradq';
        graddisploc = QT*H*QT';
        stressloc   = QT*sigma(voit2ind)*QT';

        % ++++++++++++++++++
        %  Auxiliary fields
        % ++++++++++++++++++ 
        
        xp    = QT*(Gpt-xTip)';           % local coordinates
        r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
        theta = atan2(xp(2),xp(1));

        K1 = 1.0 ;
        K2 = K1  ;

        mu = E/(2.+ nu + nu);
        kappa = 3-4*nu;    %Kolosov coeff, Plain strain
        
        SQR  = sqrt(r);
        CT   = cos(theta);
        ST   = sin(theta);
        CT2  = cos(theta/2);
        ST2  = sin(theta/2);
        C3T2 = cos(3*theta/2);
        S3T2 = sin(3*theta/2);

        drdx = CT;
        drdy = ST;
        dtdx = -ST/r;
        dtdy = CT/r;

        FACStress1 = sqrt(1/(2*pi));
        FACStress2 = FACStress1;

        FACDisp1 = sqrt(1/(2*pi))/(2*mu);
        FACDisp2 = FACDisp1;

        AuxStress   = zeros(2,2);
        AuxGradDisp = zeros(2,2);
        AuxEps      = zeros(2,2);

        for mode = 1:2
            if (mode == 1)
           
                AuxStress(1,1) = K1*FACStress1/SQR*CT2*(1-ST2*S3T2);
                AuxStress(2,2) = K1*FACStress1/SQR*CT2*(1+ST2*S3T2);
                AuxStress(1,2) = K1*FACStress1/SQR*ST2*CT2*C3T2;
                AuxStress(2,1) = AuxStress(1,2);

                u1    = K1*FACDisp1*SQR*CT2*(kappa - CT);
                du1dr = K1*FACDisp1*0.5/SQR*CT2*(kappa - CT);
                du1dt = K1*FACDisp1*SQR*(-0.5*ST2*(kappa - CT) + CT2*ST);

                u2    = K1*FACDisp1*SQR*ST2*(kappa - CT);
                du2dr = K1*FACDisp1*0.5/SQR*ST2*(kappa - CT);
                du2dt = K1*FACDisp1*SQR*(0.5*CT2*(kappa - CT) + ST2*ST);

                AuxGradDisp(1,1) = du1dr*drdx + du1dt*dtdx;
                AuxGradDisp(1,2) = du1dr*drdy + du1dt*dtdy;
                AuxGradDisp(2,1) = du2dr*drdx + du2dt*dtdx;
                AuxGradDisp(2,2) = du2dr*drdy + du2dt*dtdy;

                AuxEps(1,1) = AuxGradDisp(1,1);
                AuxEps(2,1) = 0.5*(AuxGradDisp(2,1) + AuxGradDisp(1,2));
                AuxEps(1,2) = AuxEps(2,1);
                AuxEps(2,2) = AuxGradDisp(2,2);                
            elseif (mode == 2)
                AuxStress(1,1) = -K2*FACStress2/SQR*ST2*(2-CT2*C3T2);
                AuxStress(2,2) = K2*FACStress2/SQR*ST2*CT2*C3T2;
                AuxStress(1,2) = K2*FACStress2/SQR*CT2*(1-ST2*S3T2);
                AuxStress(2,1) = AuxStress(1,2);

                u1    = K2*FACDisp2*SQR*ST2*(kappa + 2 + CT);
                du1dr = K2*FACDisp2*0.5/SQR*ST2*(kappa + 2 + CT);
                du1dt = K2*FACDisp2*SQR*(0.5*CT2*(kappa + 2 + CT) - ST2*ST);

                u2    = -K2*FACDisp2*SQR*CT2*(kappa - 2 + CT);
                du2dr = -K2*FACDisp2*0.5*(1/SQR)*CT2*(kappa - 2 + CT);
                du2dt = -K2*FACDisp2*SQR*(-0.5*ST2*(kappa - 2 + CT) - CT2*ST);

                AuxGradDisp(1,1) = du1dr*drdx + du1dt*dtdx;
                AuxGradDisp(1,2) = du1dr*drdy + du1dt*dtdy;
                AuxGradDisp(2,1) = du2dr*drdx + du2dt*dtdx;
                AuxGradDisp(2,2) = du2dr*drdy + du2dt*dtdy;

                AuxEps(1,1) = AuxGradDisp(1,1);
                AuxEps(2,1) = 0.5*(AuxGradDisp(2,1) + AuxGradDisp(1,2));
                AuxEps(1,2) = AuxEps(2,1);
                AuxEps(2,2) = AuxGradDisp(2,2);
            end
            
            % +++++++++++++++
            %   J integral
            % +++++++++++++++
            I1= (stressloc(1,1) * AuxGradDisp(1,1) + stressloc(2,1) * AuxGradDisp(2,1) ) * gradqloc(1) + ...
                (stressloc(1,2) * AuxGradDisp(1,1) + stressloc(2,2) * AuxGradDisp(2,1) ) * gradqloc(2);

            I2= (AuxStress(1,1) * graddisploc(1,1) + AuxStress(2,1) * graddisploc(2,1) ) * gradqloc(1) + ...
                (AuxStress(2,1) * graddisploc(1,1) + AuxStress(2,2) * graddisploc(2,1) ) * gradqloc(2);

            StrainEnergy = 0;
            for i=1:2 %size(AuxEpsm1,1)
                for j=1:2  %size(AuxEpsm1,2)
                    StrainEnergy= StrainEnergy+  stressloc(i,j)*AuxEps(i,j);
                end
            end
            
            % Interaction integral I
            I(mode,1) = I(mode,1) + (I1 + I2 - StrainEnergy*gradqloc(1))*det(J0)*wt;
        end   %loop on mode

    end       % of quadrature loop
end           % end of element loop

% Compute SIFs from I integral
Knum = I.*E/(2*(1-nu^2)); % plain strain 

% Compute the exact SIFs

Kexact = [34 4.55]

% -------------------------------------------------------------------------

