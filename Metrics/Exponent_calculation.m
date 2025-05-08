function [B ,log_exponent_all] =Exponent_calculation(filepath_ascii)
[edge_network,vert_network,point_network,edge, point,vertex]=ultimate_amira_read(filepath_ascii);

%% Check and fix edges
s=edge_network.EdgeConnectivity_EDGE{1}(:,1);
     t=edge_network.EdgeConnectivity_EDGE{1}(:,2);
     G=digraph(s+1,t+1);
     G.Nodes.Isterminal = G.outdegree == 0;
     figure
     plot(G, 'NodeCdata', G.Nodes.Isterminal+1,'Layout','layered');
     rootIDs= find(G.outdegree==0)-1;

 corrected=0;

     if size(rootIDs,1)>1 % check if there is more than one root node.
         corrected=1;
         prompt='input the correct Root_node ID'

         answer = (inputdlg(prompt))
         answer=str2num(cell2mat(answer))
         if ismember(answer,rootIDs)
             rootIDs=answer;
         end
         [new_edge_nodes,bad_edge_indices] =Find_bad_edges(cell2mat(edge_network.EdgeConnectivity_EDGE),rootIDs,0);
         %clear network
        network=table(new_edge_nodes);   %% put the new nodes into the table network to be used for the rest of the Strahler stuff
        network.Properties.VariableNames= [{'edge_nodes'} ];
     else
         original_edges=cell2mat(edge_network.EdgeConnectivity_EDGE);
         network=table([original_edges]);
         network.Properties.VariableNames= [{'edge_nodes'} ];
     end
s2=network.edge_nodes(:,1);
t2=network.edge_nodes(:,2);
G_reverse=digraph(t2+1,s2+1);
G_reverse.Nodes.Isterminal = G_reverse.outdegree == 0;
figure
plot(G_reverse, 'NodeCdata', G_reverse.Nodes.Isterminal+1,'Layout','layered');


edges = unique(network.edge_nodes);
counts = histc(network.edge_nodes(:), edges); %this is the co-ordination number
idx = find(counts == 1);  % find the branching point nodes

exponent_data=NaN([length(edges),2]);

for i= 1:length(edges)
    downstream=nearest(G_reverse,edges(i)+1,Inf);
    %find the segment where that is 
    ind2=find(network.edge_nodes(:,1)==edges(i));
    rad=edge_network.MeanRadius_EDGE{1}(ind2);
    %sum the number of terminal ends in the downstream each node check if each node's co-ordination number
    ds_tips=sum(counts(downstream)==1);

    exponent_data(i,:)=([ds_tips,rad]);
end

exponent_data(exponent_data == 0) = NaN;

log_exponent=log(exponent_data);
log_exponent_all=vertcat(log_exponent_cc1,log_exponent_cc4,log_exponent_cc9);
%% concaternate all the log exponent datasets now from the different connected components.

[B,BINTR,BINTJM]=gmregress(log_exponent_all(:,1),log_exponent_all(:,2),0.05); 

end
