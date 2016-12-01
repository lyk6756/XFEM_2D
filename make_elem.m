function element=make_elem(node_pattern,num_u,num_v,inc_u,inc_v)

% function element=make_elem(node_pattern,num_u,num_v,inc_u,inc_v)
%
% creates a connectivity list

if ( nargin < 5 )
   disp(['Not enough parameters specified for make_elem function'])
end

inc=[zeros(1,size(node_pattern,2))];  %设置增量矩阵
e=1;%给单元编号赋初值，不能放到循环内部，否则每次循环都会赋为1
element=zeros(num_u*num_v,size(node_pattern,2));%定义单元矩阵大小

for row=1:num_v%对列循环
   for col=1:num_u%对行循环
      element(e,:)=node_pattern+inc;
      inc=inc+inc_u;
      e=e+1;
   end
   inc=row*inc_v;
end  %对于长方形区域，单元编号增量会有间断，设置上面的嵌套循环比较方便，熟悉这种编程技巧
