function plot_field(X,connect,elem_type,field)
                  %(变形后的节点坐标，element,elemType,平面内各节点的应力)
% function plot_field(X,connect,elem_type,field)
 
if ( nargin == 4 )
  nodesoff=0;  %变量名
end
  
if ( size(field) == size(connect) )
  elementalField=1;
else
  elementalField=0;  
end

% fill X if needed  补0
if (size(X,2) < 3)
   for c=size(X,2)+1:3
      X(:,c)=[zeros(size(X,1),1)];
   end
end

holdState=ishold;
hold on

% plot elements
if ( strcmp(elem_type,'Q9') )      % Q9 element
  ord=[1,5,2,6,3,7,4,8,1];
elseif ( strcmp(elem_type,'T3') )  % T3 element
  ord=[1,2,3,1];
elseif ( strcmp(elem_type,'T4') )  % T4 element
  ord=[1,2,3,1];
elseif ( strcmp(elem_type,'T6') )  % T6 element
  ord=[1,4,2,5,3,6,1];
elseif ( strcmp(elem_type,'Q4') )  % Q4 element
  ord=[1,2,3,4,1];
elseif ( strcmp(elem_type,'Q8') )  % T3 element
   ord=[1,5,2,6,3,7,4,8,1];
elseif ( strcmp(elem_type,'L2') )  % L2 element
  ord=[1,2];   
elseif ( strcmp(elem_type,'L3') )  % L3 element
  ord=[1,3,2];   
end

for e=1:size(connect,1)  %单元循环
  
   xpt=X(connect(e,ord),1);
   ypt=X(connect(e,ord),2);      
   zpt=X(connect(e,ord),3);%单元四个节点的坐标，次序是[1，2，3，4，1]
   
   if ( elementalField )
     fpt=field(e,ord); %field是应力场
   else
     fpt=field(connect(e,ord));
   end
   
   fill3(xpt,ypt,zpt,fpt) %FILL(X,Y,C) fills the 2-D polygon defined by vectors X and Y with the color specified by C.
end

shading interp  %Color shading mode
axis equal
      
if ( ~holdState )
  hold off
end
