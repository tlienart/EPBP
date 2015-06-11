#
# additional support functions
#
get_message(from,to) = messages[get_edge_idx(from,to),:]
#
# --------------------------------------------------------
#
# PBP_NODE_UPDATE(NODE):
#   Update of a node following the PBP method where
#   the proposals are the last estimated beliefs
#   sampled using Metropolis Hastings
#
function pbp_node_update(node)
    #
    neighbors = node_info(node)
    K         = length(neighbors)
    #
    # STEP 1(A): sample from proposal
    #
    node_p     = particles[node,:]  # (size 1,N)
    old_belief = b_evals[node,:]    # (size 1,N)
    #
    for iter=1:MHIter
        cand_p     = sampleMHP(node_p) # (size 1,N)
        cur_belief = pbp_eval_belief(node,cand_p)
        #
        # acceptance ratio
        #
        alpha  = cur_belief./old_belief
        accept = rand(1,N).<alpha
        #
        node_p[1,accept[:]]     = cand_p[accept[:]]
        old_belief[1,accept[:]] = cur_belief[accept]
        #
    end
    old_belief           /= sum(old_belief)
    pbp_particles[node,:] = node_p
    pbp_beliefs[node,:]   = old_belief
    #
    # STEP 1(B): evaluate incoming messages at new points
    #
    for k = 1:K
        neighb   = neighbors[k]
        neighb_p = particles[neighb,:]
        mess     = zeros(1,N)
        #
        for j=1:N # incoming message
            tmp     = eval_edge_pot(node,neighb,node_p[j],neighb_p)
            tmp   .*= eval_node_pot(neighb,neighb_p)
            tmp   ./= get_message(node,neighb)
            mess[j] = sum(tmp)
        end
        mess /= sum(mess)
        # store
        messages[get_edge_idx(neighb,node),:] = mess
    end
    #
    # STEP 2: evaluate outgoing messages
    #
    for k = 1:K
        neighb   = neighbors[k]
        neighb_p = particles[neighb,:]
        mess     = zeros(1,N)
        #
        for j=1:N
            tmp     = eval_edge_pot(neighb,node,neighb_p[j],node_p)
            tmp   .*= eval_node_pot(node,node_p)
            tmp   ./= get_message(neighb,node)
            mess[j] = sum(tmp)
        end
        mess /= sum(mess)
        # store
        messages[get_edge_idx,:] = mess
    end
end