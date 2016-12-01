function [W,Q]=discontQ4quad(order,phi)


corner = [1 2 3 4 1];
node   = [-1 -1; 1 -1; 1 1; -1 1];

% loop on element edges
for i = 1 : 4
    n1 = corner(i);
    n2 = corner(i+1);
    if ( phi(n1)*phi(n2) < 0 )
        r    = phi(n1)/(phi(n1)-phi(n2));
        pnt  = (1-r)*node(n1,:)+r*node(n2,:);
        node = [node;pnt]; %将裂纹与单元边的交点也列出来
    end
end

% get decompused triangles
tri = delaunay(node(:,1),node(:,2));%一中三角化的方法，将已知坐标的一些点连接成一系列的三角形，保证三角形的外接圆内没有其他坐标点
tri = tricheck(node,tri);

% loop over subtriangles to get quadrature points and weights
pt = 1;
for e = 1:size(tri,1)   %还是一到四
    [w,q]=quadrature(order,'TRIANGULAR',2); % 这里order==2
    % transform quadrature points into the parent element
    coord = node(tri(e,:),:);
    a = det([coord,[1;1;1]])/2;
    if ( a<0 )  % need to swap connectivity
        coord = [coord(2,:);coord(1,:);coord(3,:)];   %为什么又要进行行的调整，估计是两者都要满足，雅克比行列式和面积都不能为负值
        a = det([coord,[1;1;1]])/2;   %好象是三角形的面积
    end

    if ( a~=0 )
        for n=1:length(w)
            N=lagrange_basis('T3',q(n,:));
            Q(pt,:) = N'*coord;
            W(pt,1) = 2*w(n)*a;      %2a的值等于雅可比行列式的值
            pt = pt+1;
        end
    end

end








