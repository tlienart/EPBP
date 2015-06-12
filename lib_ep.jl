#
# EP_NODE_UPDATE(NODE):
#   Update of a node following LBP on discretized grid.
#
function ep_node_update(node)
	neighbors = get_neighbors(node)
	K  		  = length(neighbors)
	#
	for k=1:K
		neighb 	= neighbors[k]
		# cavity distr
		node_cav   = q_moments[node,:]
		neighb_cav = q_moments[neighb,:]
		eta_out    = get_edge_eta(node,neighb)
		eta_in 	   = get_edge_eta(neighb,node)
		#
		if bool(prod(eta_out))
			node_cav = normal_div(node_cav,eta_out)
		end
		if bool(prod(eta_in))
			neighb_cav = normal_div(neighb_cav,eta_in)
		end
		#
		node_cav_distr   = Normal(node_cav[1],node_cav[2])
		neighb_cav_distr = Normal(neighb_cav[1],neighb_cav[2])
		#
		eval_grid = zeros(Ninteg,Ninteg)
		for i=1:Ninteg
			xi = integ_pts[i]
			for j=1:Ninteg
				xj = integ_pts[j]
				eval_grid[i,j] = pdf(node_cav_distr,xi) *
									pdf(neighb_cav_distr,xj) *
										eval_edge_pot(node,neighb,xi,xj)
			end
		end
		#
		node_marg  	 = sum(eval_grid,2)
		node_marg 	/= sum(node_marg)
		neighb_marg  = sum(eval_grid,1)
		neighb_marg /= sum(neighb_marg)
		#
		mom_node_marg   = params(fit_mle(Normal,integ_pts,node_marg))
		mom_neighb_marg = params(fit_mle(Normal,integ_pts,neighb_marg))
		#
		if mom_node_marg[2] < node_cav[2]
			new_eta_in 				= normal_div(mom_node_marg,node_cav)
			q_moments[node,:] 		= normal_prod(node_cav,new_eta_in)
			edge_idx 				= get_edge_idx(neighb,node)
			eta_moments[edge_idx,:] = [new_eta_in[1] new_eta_in[2]]
		end
		if mom_neighb_marg[2] < neighb_cav[2]
			new_eta_out 			= normal_div(mom_neighb_marg,neighb_cav)
			q_moments[neighb,:]		= normal_prod(neighb_cav,new_eta_out)
			edge_idx 				= get_edge_idx(node,neighb)
			eta_moments[edge_idx,:] = [new_eta_out[1] new_eta_out[2]]
		end
	end
end
