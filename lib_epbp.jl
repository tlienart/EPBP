#
# additional support functions
#
get_edge_weights(from,to) = e_weights[get_edge_idx(from,to),:]
get_edge_eta(from,to)     = eta_moments[get_edge_idx(from,to),:]
#
# --------------------------------------------------------
#
# EPBP_NODE_UPDATE(NODE):
#   Update of a node following the BP/EP method
#
function epbp_node_update(node,fastmode=false)
    neighbors = get_neighbors(node)
    K         = length(neighbors)
    #
    # STEP 1: SAMPLING
    #   >> sample points at node from current proposal
    #
    # > sample (Normal)
    node_p = q_moments[node,1] + q_moments[node,2] * randn(1,N) # size (1,N)
    # > store
    particles[node,:] = node_p # size (1,N)
    #
    # STEP 2: EVALUATE INCOMING MESSAGES + LOCAL BELIEF (Normal)
    #
    inmess_eval,belief_eval = epbp_eval_belief(node,node_p,fastmode)
    # > compute importance weights
    proposal_eval  = pdf(Normal(q_moments[node,1],q_moments[node,2]),node_p)
    belief_weights = belief_eval./proposal_eval
    # > normalize to avoid under/over - flow
    belief_weights /= sum(belief_weights)
    # > store
    b_weights[node,:] = belief_weights
    b_evals[node,:]   = belief_eval
    #
    # STEP 3: EVALUATE OUTGOING MESSAGES
    #
    for k=1:K
        neighb   = neighbors[k]
        # > M_uv = B_u/m_vu = for outmess m_uv
        tmp_e_w  = belief_weights./inmess_eval[k,:]  
        # > normalize to avoid under/over - flow
        tmp_e_w /= sum(tmp_e_w)
        # > store
        e_weights[get_edge_idx(node,neighb),:] = tmp_e_w
    end
    #
    # STEP 4a: EP PROJECTION - PtA (node)
    #
    # 
    # STEP 4b: EP PROJECTION - PtB (node)
    #
    M            = length(integ_pts)
    outmess_eval = zeros(K,M)
    #
    for k=1:K
        neighb   = neighbors[k]
        cavity   = q_moments[neighb,:]
        eta_out  = get_edge_eta(node,neighb)
        prev_mom = cavity
        if bool(prod(eta_out)) # initially they're set to 0
            cavity = normal_div(cavity,eta_out) # this should always be ok
        end
        outmess_eval_tp = epbp_eval_message(node,neighb,integ_pts)
        #
        # here: hack, using MLE (should not)
        #
        eta_out_new         = params(fit_mle(Normal,integ_pts,outmess_eval_tp))
        eidx                = get_edge_idx(node,neighb)
        q_moments[neighb,:] = normal_prod(cavity,eta_out_new)
        # > store new eta
        eta_moments[eidx,:] = [m for m in eta_out_new]
        #
    end
end
#
# --------------------------------------------------------
#
# EPBP_EVAL_BELIEF(NODE,EVAL_POINTS):
#   Evaluate the current estimator of the beliefs at
#   a given node and for given points.
#   For that, all the incoming messages are evaluated
#   and the product is taken.
#
#   Complexity: O(K*N*M) where
#       M = length(eval_points) (typically N)
#       N = number of samples per node
#       K = number of neighbors
#
function epbp_eval_belief(node,eval_points,fastmode=false)
    neighbors   = get_neighbors(node)
    K,M         = length(neighbors),length(eval_points)
    inmess_eval = zeros(K,M)
    #
    for k=1:K
        inmess_eval[k,:] = epbp_eval_message(neighbors[k],node,eval_points,fastmode) # size (1,M)
    end
    return inmess_eval, eval_node_pot(node,eval_points).*prod(inmess_eval,1)
end
#
# --------------------------------------------------------
#
# EPBP_EVAL_MESSAGE(FROM,TO,EVAL_POINTS):
#   Evaluate the current message estimator FROM=>TO
#   at given points.
#
#   Complexity: O(N*M) where
#       M = length(eval_points) (typically N)
#       N = number of samples per node
#
function epbp_eval_message(from,to,eval_points,fastmode=false)
    M,from_p = length(eval_points), particles[from,:]
    #
    curedge_w    = get_edge_weights(from,to) # size (1,N)
    message_eval = zeros(1,M)                # size (1,M)
    #
    if fastmode
        idx_comp = rand(Multinomial(C,curedge_w[:]))
        comp_idx = int(zeros(C,1))
        i = 1
        for k=1:length(idx_comp)
            for l=1:idx_comp[k]
                comp_idx[i]  = k
                i           += 1
            end
        end
        for comp=1:C 
            message_eval += eval_edge_pot(from,to,from_p[comp_idx[comp]],eval_points)
        end
    else
        for i=1:M
            message_eval[i] += sum(curedge_w.*eval_edge_pot(from,to,eval_points[i],from_p))
        end
    end
    return message_eval # size (1,M)
end