function [Jdomain,qnode,radius] = Jdomain(tip_elem,xTip)

global node element

numnode = size(node,1);
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
y3 = x(4,2);  %是裂尖单元四个节点的整体坐标（逆时针）

A1 = 0.5 * ((x0-x2)*(y1-y2) - (x1-x2)*(y0-y2)) ;%%？？？？？？？？？
A2 = 0.5 * ((x0-x3)*(y2-y3) - (x2-x3)*(y0-y3)) ;
area = A1 + A2 ;%裂尖单元的面积

% J radius = fac * sqrt(area);
fac    = 2.5;
radius = fac * sqrt(area);
center = xTip; %以裂纹尖端为圆心，以一个常数乘以面积开方为半径

r=[];
% Distance from the center of tip element
for i = 1 : numnode  %所有的节点
    sctr = node(i,:); %所有节点的整体坐标
    rho  = sqrt((sctr(1)-center(1))^2+(sctr(2)-center(2))^2);%任一节点和裂尖的距离
    r    = [r,rho];
end
test = r-radius;%这里是一个列向量，每个元素是距离和半径的差值
test = test(element)';%这里是一个253*4的矩阵，每一行元素是一个单元各个节点相对应的值（转置）
test = max(test).*min(test);
Jdomain = find(test<=0);%J积分区域就是这样一个单元的集合，即刚才画的圆所穿过单元（有22个）的集合，而不包括完全在圆内的单元
test1 = r-radius;
test1 = test1(element(Jdomain,:))';%这里是一个22*4的矩阵，即每个单元节点对应的值（转置）
test1 = (test1<=0);
qnode = test1'; %是一个逻辑矩阵，22*4%%%%%%%%%%22个单元的节点在圆内的置1，否则置0

