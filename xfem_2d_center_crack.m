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
L = 2 ;
D = 6 ;

% Material properties
E  = 1e3 ;
nu = 0.3 ;
stressState='PLANE_STRAIN';

% Crack properties
% Crack is represented as straight segments
% with points stored in matrix xCr
% Array xTip contains coordinates of the crack tip
% Array seg is the segment containg the tip
% It serves in computing alpha, the inclination angle
% of the crack w.r.t the X axis
a = 0.25 ;                     % crack length
xCr   = [-a 0; a 0];
xTip  = [-a 0; a 0];
seg1   = xCr(2,:) - xCr(1,:);   % tip segment
alpha1 = atan2(seg1(2),seg1(1));  % inclination angle
QT1    = [cos(alpha1) sin(alpha1); -sin(alpha1) cos(alpha1)];
seg2   = xCr(1,:) - xCr(2,:);   % tip segment
alpha2 = atan2(seg2(2),seg2(1));  % inclination angle
QT2    = [cos(alpha2) sin(alpha2); -sin(alpha2) cos(alpha2)];

% Loading
sigmato = 1  ;

% ---------------------------------------

% +++++++++++++++++++++++++++++
%            MESHING
% +++++++++++++++++++++++++++++

% ---------------------------------------
disp([num2str(toc),'   MESH GENERATION'])
% Number of nodes along two directions
nnx = 24 ;
nny = 50 ;

% Four corner points
pt1 = [ -L/2 -D/2];
pt2 = [  L/2 -D/2];
pt3 = [  L/2  D/2];
pt4 = [ -L/2  D/2];

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
t   = 1/norm(seg1)*seg1;
for i = 1 : numnode
    x = node(i,1);
    y = node(i,2);
    l   = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)) ;
    phi = (y0-y1)*x + (x1-x0)*y + (x0*y1-x1*y0);
    ls(i,1) = phi/l;            % normal LS
    ls(i,2) = ([x y]-xTip(1))*t';  % tangent LS
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
    sctrB = assembly(iel,enrich_node,pos);
    
    % ---------------------
    % Loop on Gauss points
    % ---------------------
    for kk = 1 : size(W,1)
        pt = Q(kk,:);                             % quadrature point
        % B matrix
        [B,J0] = xfemBmatrix(pt,elemType,iel,enrich_node,xCr,xTip,alpha);
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
fac = 50;
plot_mesh(node+fac*[u_x u_y],element,elemType,'b-');
title(' Numerical deformed configuration ')
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
        [B,J0] = xfemBmatrix(pt,elemType,e,enrich_node,xCr,xTip,alpha);
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

% Determine J domain and weight function
[Jdomain,qnode,radius] = Jdomain(tip_elem,xTip);

% plot
figure
hold on
% plot the circle
theta = -pi:0.1:pi;
xo = xTip(1) + radius*cos(theta) ;
yo = xTip(2) + radius*sin(theta) ;
plot(xo,yo,'k-');
plot_mesh(node,element,'Q4','b-')
plot_mesh(node,element(Jdomain,:),'Q4','r-')
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
for iel = 1 : size(Jdomain,2)
    e = Jdomain(iel) ; % current element
    sctr = element(e,:);
    nn   = length(sctr);
    % Choose Gauss quadrature rule
    if (ismember(e,split_elem))     % split element
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
        dNdx  = dNdxi*invJ0;
        Gpt = N' * node(sctr,:);     % GP in global coord
        
        % +++++++++++++++++++++++++ 
        % Gradient of displacement
        % +++++++++++++++++++++++++ 
        % need to compute u,x u,y v,x v,y, stored in matrix H
        [B,J0] = xfemBmatrix(pt,elemType,e,enrich_node,xCr,xTip,alpha);
        leB = size(B,2);
        
        % nodal displacement of current element
        % taken from the total nodal parameters u
        U = element_disp(e,pos,enrich_node,u);
        
        % compute derivatives of u w.r.t xy
        H(1,1) = B(1,1:2:leB)*U(1:2:leB);    % u,x
        H(1,2) = B(2,2:2:leB)*U(1:2:leB);    % u,y
        H(2,1) = B(1,1:2:leB)*U(2:2:leB);    % v,x
        H(2,2) = B(2,2:2:leB)*U(2:2:leB);    % v,y
        
        % ++++++++++++++
        % Gradient of q
        % ++++++++++++++        
        q     = qnode(iel,:);
        gradq = q*dNdx;

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

F = 1.12 - 0.23*(a/L) + 10.55*(a/L)^2 - 21.72*(a/L)^3 + 30.39*(a/L)^4; 
Kexact = F * sigmato*sqrt(pi*a);

% -------------------------------------------------------------------------

