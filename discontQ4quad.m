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
        node = [node;pnt]; %�������뵥Ԫ�ߵĽ���Ҳ�г���
    end
end

% get decompused triangles
tri = delaunay(node(:,1),node(:,2));%һ�����ǻ��ķ���������֪�����һЩ�����ӳ�һϵ�е������Σ���֤�����ε����Բ��û�����������
tri = tricheck(node,tri);

% loop over subtriangles to get quadrature points and weights
pt = 1;
for e = 1:size(tri,1)   %����һ����
    [w,q]=quadrature(order,'TRIANGULAR',2); % ����order==2
    % transform quadrature points into the parent element
    coord = node(tri(e,:),:);
    a = det([coord,[1;1;1]])/2;
    if ( a<0 )  % need to swap connectivity
        coord = [coord(2,:);coord(1,:);coord(3,:)];   %Ϊʲô��Ҫ�����еĵ��������������߶�Ҫ���㣬�ſ˱�����ʽ�����������Ϊ��ֵ
        a = det([coord,[1;1;1]])/2;   %�����������ε����
    end

    if ( a~=0 )
        for n=1:length(w)
            N=lagrange_basis('T3',q(n,:));
            Q(pt,:) = N'*coord;
            W(pt,1) = 2*w(n)*a;      %2a��ֵ�����ſɱ�����ʽ��ֵ
            pt = pt+1;
        end
    end

end








