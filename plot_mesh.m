function plot_mesh(X,connect,elem_type,se)
                 %(node,element,elemType,'b-'线的类型，子程序中用变量se表示)
% function plot_mesh(X,connect,elem_type,linespec)
% 
% plots a nodal mesh and an associated connectivity.  X is
% teh nodal coordinates, connect is the connectivity, and
% elem_type is either 'L2', 'L3', 'T3', 'T6', 'Q4', or 'Q9' 
% depending on the element topology.

%(hold on:使当前坐标系和图形保留；hold off:使当前坐标系和图形不保留；hold:在以上两个命令中切换)  
if ( nargin < 4 )
   se='w-';
end

holdState=ishold; %如果坐标系和图形保留ishold值取1，如果不保留ishold值取0（When HOLD is ON, the current plot and all axis properties
                  %are held so that subsequent graphing commands add to the existing graph）
hold on

% fill X if needed
if (size(X,2) < 3)  %二维或一维
   for c=size(X,2)+1:3
      X(:,c)=[zeros(size(X,1),1)];
   end
end

for e=1:size(connect,1)  %总的单元个数
  
   if ( strcmp(elem_type,'Q9') )       % 9-node quad element
      ord=[1,5,2,6,3,7,4,8,1];
   elseif ( strcmp(elem_type,'Q8') )  % 8-node quad element
      ord=[1,5,2,6,3,7,4,8,1];
   elseif ( strcmp(elem_type,'T3') )  % 3-node triangle element
      ord=[1,2,3,1];
   elseif ( strcmp(elem_type,'T6') )  % 6-node triangle element
      ord=[1,4,2,5,3,6,1];
   elseif ( strcmp(elem_type,'Q4') )  % 4-node quadrilateral element
      ord=[1,2,3,4,1];     %逆时针
   elseif ( strcmp(elem_type,'L2') )  % 2-node line element
      ord=[1,2];   
   elseif ( strcmp(elem_type,'L3') )  % 3-node line element
      ord=[1,3,2];   
   elseif ( strcmp(elem_type,'H4') )  % 4-node tet element
      ord=[1,2,4,1,3,4,2,3];   
   elseif ( strcmp(elem_type,'B8') )  % 8-node brick element
      ord=[1,5,6,2,3,7,8,4,1,2,3,4,8,5,6,7];   
   end
   
   for n=1:size(ord,2)
      xpt(n)=X(connect(e,ord(n)),1);
      ypt(n)=X(connect(e,ord(n)),2);      
      zpt(n)=X(connect(e,ord(n)),3);
   end
   plot3(xpt,ypt,zpt,se)
end

rotate3d on
axis equal  %sets the aspect ratio so that equal tick mark(核对符号) increments on the x-,y- and z-axis are equal in size.？
      
if ( ~holdState )
  hold off
end
