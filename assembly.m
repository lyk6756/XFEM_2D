function sctrB = assembly(e,enrich_node,pos)

global node element

sctr = element(e,:); %四个节点的序号（逆时针）
nn   = length(sctr);

for k = 1 : nn
    sctrBfem(2*k-1) = 2*sctr(k)-1 ;
    sctrBfem(2*k)   = 2*sctr(k)   ;
end

if ( any(enrich_node(sctr)) == 0 ) % Non-enriched elements四个数全是零的时候
    sctrB = sctrBfem ;
else
    sn = size(find(enrich_node(sctr) == 1),1); %一个单元中有H(x)加强的节点个数
    tn = size(find(enrich_node(sctr) == 2),1); %一个单元中有裂尖加强的节点个数
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

