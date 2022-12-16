%% Load data and see the iterative walks until 6
load("email-Eu-core.mat")
% load("wiki-topcats.mat")
%% Save the data to the A_1
A_1 = Problem.A;
%% Remove the self-loops
diag_A = diag(diag(A_1));
spy(diag_A);
A_1 = A_1 - diag_A;
%% Look at the removed graph
diag_A = diag(diag(A_1));
spy(diag_A);
%% Remove the isolated nodes
A_1_graph = digraph(A_1);
indeg_A_1 = indegree(A_1_graph);
outdeg_A_1 = outdegree(A_1_graph);
mixdeg_A_1 = indeg_A_1 + outdeg_A_1;
zero_pos = find(~mixdeg_A_1);
%% Remove the empty nodes
A_1_graph = rmnode(A_1_graph, zero_pos);
A_1 = adjacency(A_1_graph);
%% PLot the graph
A_1_graph = digraph(A_1);
plot(A_1_graph)
%% Continue the walks until N steps/hops
A = A_1;
A_org = A_1;
A_sum = A_org;
N = 4;
while N > 1
    A = A*A_org;
    A_sum = A_sum + A;
    N = N - 1;
end
zero_elements_num = length(A_sum)*length(A_sum) - nnz(A_sum);
%% Here is the idea, I will see what nodes are still not connected on step five,
% excluding the 0 indegree and 0 outdegree nodes, it can also be step 4,5,6
A_sum_Connectivity = (A_sum ~= 0);
A_sum_Connectivity_graph = digraph(A_sum_Connectivity);
A_sum_Connectivity_OutDegree = outdegree(A_sum_Connectivity_graph);
A_sum_Connectivity_InDegree = indegree(A_sum_Connectivity_graph);
zero_pos_InDegree = find(~A_sum_Connectivity_InDegree);
zero_pos_OutDegree = find(~A_sum_Connectivity_OutDegree);
zero_pos = find(~A_sum_Connectivity);
%% The point is to assume that nodes are not connected in the step 4,5,6 
% are very unlikely to be in the same community, excluding those have no
% out/in at all nodes
% So first step is to find unconnected nodes list
% We iterate through all the nodes(zero ones), if it belongs to zero_InDegree or
% zero_Outdegree, we set it to be 0, otherwise we set it to be 1 
% the one in Prob_A means node i is not connected to j, (but not in special
% case where indegree = 0 or outdegree  = 0)
Prob_A = zeros(length(A_1),length(A_1));
for i = 1:length(zero_pos)
% for i = 1:900
    x = ceil(zero_pos(i)/length(A_1));
%     disp("x = "+x)
    y = zero_pos(i) - (x-1)*length(A_1);
%     disp("y = "+y)
%     break
    if ismember(y, zero_pos_InDegree) | ismember(x, zero_pos_OutDegree)
%         disp("x = "+x+" y = "+y+" continued");
        continue;
    else
%         disp("x = "+x+" y = "+y+" unconnected");
        Prob_A(x,y) = 1;
    end
end
%% Now we have the unconnected node pair in five steps, let's see how many
% zero edges between nodes;
% Afterwards, we will want to record the unique node index that has zero
% value in i,j as well as j,i, which represents pair of nodes that are not
% connected with each other in either way (in and out)
total_zeros = sum(Prob_A, 'all');
Prob_A_sym = floor((Prob_A + Prob_A.')/2);
L = tril(Prob_A_sym);
L_total_num = sum(L, 'all');
Disconnected_node_id = find(L);
%% Record the position of those nodes, and return unique node indexes
unique_disconnected_node = [0];
index = 1;
for i = 1:length(Disconnected_node_id)
    x = ceil(Disconnected_node_id(i)/length(A_1));
    y = Disconnected_node_id(i) - (x-1)*length(A_1);
    unique_disconnected_node(index) = x;
    unique_disconnected_node(index+1) = y;
    index = index + 2;
end
% The total number of disconnected nodes
size(unique_disconnected_node);
% The total number of unique disconnected nodes
unique_disconnected_node = unique(unique_disconnected_node);
size(unique_disconnected_node);
%% With the unique disconnected nodes, we assign them into different 
% community labels, i.e. from 1 to M with M equals to the distinct nodes
Community_label = 1;
Community_Prob_A = zeros(length(A_1), length(unique_disconnected_node));
for y = 1:length(A_1)
    if ismember(y, unique_disconnected_node)
        Community_Prob_A(y,Community_label) = 1.0;
%         disp("y = "+y+" Assign community label "+Community_label);
        Community_label = Community_label + 1;
    else
%         disp("Assign default label "+1.0/length(Class_Prob))
        Community_Prob_A(y,:) = Community_Prob_A(y,:) + 1.0/length(unique_disconnected_node);
    end   
end
%% Now the community labels are assigned, we want to calculate from node
% N to 1, we start from adding probabilities to its neibour but
% proportional to its outdegrees. 
EPOCHS = 100;
for epochs = 1:EPOCHS
    Community_Prob_A_temp = zeros(length(A_1), length(unique_disconnected_node));
    for node = length(A_1):-1:1
%         disp("We enter the node = "+node);
        successor_nodes = successors(A_1_graph, node);
        influencer_prob = Community_Prob_A(node, :);
        for i = 1:length(successor_nodes)
            successor_node = successor_nodes(i);
%             disp("successor_node = "+ successor_node+" is changed");
            if ismember(successor_node, unique_disconnected_node)
                continue;
            else
                % We increase the probability of the node accordingly and
                % store in the temp column
                Community_Prob_A_temp(successor_node, :) = Community_Prob_A_temp(successor_node, :) + 1.0*influencer_prob/length(successor_nodes);
%                 Community_Prob_A(successor_node, :) = Community_Prob_A(successor_node, :) + 1.0*influencer_prob/length(successor_nodes);
%                 Community_Prob_A(successor_node, :) = Community_Prob_A(successor_node, :) / sum(Community_Prob_A(successor_node, :));
%                 disp("The prob of successor = "+successor_node+" is increased to "+Community_Prob_A(successor_node, :));
            end
        end
    end
    % Update all the nodes after iteration
    for node = length(A_1):-1:1
        Community_Prob_A(node, :) = Community_Prob_A(node, :) + Community_Prob_A_temp(node, :);
        Community_Prob_A(node, :) = Community_Prob_A(node, :) / sum(Community_Prob_A(node, :) );
    end
end
[Community_Prob, Community_label] = max(Community_Prob_A, [], 2);
%% We observe the distribution of Community_label and compare it with the real distribution
histfit(Community_label);
[prob, label] = max(Community_Prob_A, [], 2);
counts = hist(label, [1:length(unique_disconnected_node)]);
%% Calculate the modularity for predicted and ground-truth labels
% First we need to generate the null model (constraining the degree)
InDegree = indegree(A_1_graph);
OutDegree = outdegree(A_1_graph);
Null_A = zeros(length(A_1), length(A_1));
Perm_A = randperm(length(A_1));
for node = 1:length(A_1)
    temp_node = Perm_A(node);
    % We connect this node to others until it reaches max outdegree
    while sum(Null_A(temp_node, :)) < OutDegree(temp_node)
        
    end
    successors_list = succcessors(A_1_graph, temp_node);
    successors_index_random_list = randperm(length(successors_list));
    successors_index_random_list_index = 1;
    % Calculate the outdegree to see if it equals to the maximum
    while sum(Null_A(temp_node, :)) < OutDegree(temp_node)
        successor_index = successors_index_random_list(successors_index_random_list_index);
        % Calculate the indegree to see if it equals to maximum
        if sum(Null_A(:, successors_list(successor_index))) < Indegree(successors_list(successor_index))
            Null_A(temp_node, sucessors_list(successor_index)) = 1;
            successor_index = successor_index + 1
        else
            successor_index = successor_index + 1
        end
    end
    
    for i = 1:length(successors_list)
        
    end
end


