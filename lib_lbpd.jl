#
# LBPD_NODE_UPDATE(NODE):
#   Update of a node following LBP on discretized grid.
#
function lbpd_node_update(node)
    neighbors = get_neighbors(node)
    K         = length(neighbors)
    #
    # STEP 1: EVALUATE LOCAL BELIEF + INCOMING MESSAGES
    #
    cur_belief = init_beliefs[node,:]
    #
    for k=1:K
        neighb = neighbors[k]
        # > compute the INCOMING (neighb>node) message
        inc_mess = lbpd_eval_message(neighb,node)
        # > multiply to cur_belief
        cur_belief .*= inc_mess
        # > store message
        eidx             = get_edge_idx(neighb,node)
        messages[eidx,:] = inc_mess
    end
    # > store (and normalize to avoid under/over flow)
    beliefs[node,:] = cur_belief/sum(cur_belief)
    #
    # STEP 2: COMPUTE OUTGOING MESSAGES
    #
    for k=1:K
        neighb   = neighbors[k]
        out_mess = lbpd_eval_message(node,neighb)
        # > store message
        eidx             = get_edge_idx(node,neighb)
        messages[eidx,:] = out_mess
    end
end
#
# LBPD_EVAL_MESSAGE(FROM,TO):
#   Evaluate message on determinstic grid
#
function lbpd_eval_message(from,to)
    mess   = zeros(1,Ngrid)
    eidx_r = get_edge_idx(to,from)
    if HOMOG_EDGE_POT
        # m_uv(x_v^j) = sum_i psi_uv(x_u^i,x_v^j)B_u(x_u^i)/m_vu(x_u^i)
        for j=1:Ngrid
            tmp     = edge_pot_grid[:,j]' # psiuv(x_u^:,x_v^j)
            tmp   .*= beliefs[from,:]     # B_u(x_u^:)
            tmp   ./= messages[eidx_r,:]  # /m_vu(x_u^:)
            mess[j] = sum(tmp)
        end
    else 
        # m_uv(x_v^j) = sum_i psi_uv(x_u^i,x_v^j)B_u(x_u^i)/m_vu(x_u^i)
        for j=1:Ngrid
            tmp     = eval_edge_pot(from,to,grid,grid[j]) # psiuv(x_u^:,x_v^j)
            tmp   .*= beliefs[from,:]                     # B_u(x_u^:)
            tmp   ./= messages[eidx_r,:]                  # /m_vu(x_u^:)
            mess[j] = sum(tmp)
        end
    end
    # > return (after normalizing to avoid overflow)
    return mess/sum(mess)
end