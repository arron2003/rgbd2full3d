function order = occlusion2depthOrdering(bndinfo, elab, glabels)
% order = occlusion2depthOrdering(bndinfo, elab, glabels)
% hacky algorithm for finding best ordering of vertical objects
% where front objects should have low numbers


order = (1:bndinfo.nseg);

currcost = Inf;

mc = 0;
improved = true;
while improved || currcost <= 0 
    improved = false;
    
    %disp('iter')
    
    for k = 1:bndinfo.nseg %1:bndinfo.nseg
        
        lastorder = order;
        
        
        for pos = 1:bndinfo.nseg      
             
            ind = find(lastorder==k);
            neworder = lastorder;
            neworder(ind) = 0;
            
            neworder = [neworder(1:pos-1) k neworder(pos:end)];
            neworder(neworder==0) = [];
            
            cost = inconsistencies(neworder, bndinfo, elab, glabels, mc);            
            
            if cost < currcost
                currcost = cost;
                improved = true;
                order = neworder;
                lastorder = order;
            end
            
        end
        
    end   
    
    if ~improved
        
        currcost = inconsistencies(order, bndinfo, elab, glabels, 0);
        
        if currcost > 0 && mc==0
            best = order;           
            break;
        elseif currcost > 0
            break;
        end
        best = order;
        if mc < numel(order)
            break;
        end        
        
        mc = mc - 3;

        disp(num2str([mc currcost]))
    end
end

order = best;

%keyboard;
%disp(num2str(currcost))
%inconsistencies(order, bndinfo, elab, glabels, 0)

%% cost is number of places away from valid ordering            
function cost = inconsistencies(order, bndinfo, elab, glabels, mc)

cost = 0;

spLR = bndinfo.edges.spLR;

for k = 1:numel(elab)
    if all(glabels(spLR(k, :))==2)
        pos1 = find(order==spLR(k, 1));
        pos2 = find(order==spLR(k, 2));
        if elab(k)==1 
            cost = cost + max(mc, pos1-pos2); %  (order(spLR(k, 1))-order(spLR(k, 2))>0); % 
        elseif elab(k)==2
            cost = cost + max(mc, pos2-pos1); % (order(spLR(k, 2))-order(spLR(k, 1))>0); %max(0, order(spLR(k, 2))-order(spLR(k, 1)));
        end
    end
end