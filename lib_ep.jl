#
# 	Code released under the MIT license, see corresponding LICENSE file
#	(c) 2015, Thibaut Lienart
#
# --------------------------------------------------------------------------------------------------
#
# EP_NODE_UPDATE(NODE):
#   Update of a node following pure EP
#
function ep_node_update(node)
	# STEP 1: EP NODE POT PROJECTION
	ep_node_proj(node)
	# STEP 2: EP EDGE POT PROJECTION
	for neighb in get_neighbors(node)
		ep_edge_proj(node,neighb)
	end
end
#
function ep_node_proj(node)
	node_cavity = q_moments[node,:]
	eta_node 	= get_node_eta(node)
	if eta_node[2]>node_cavity[2]
		node_cavity = normal_div(node_cavity,eta_node)
	end
	node_eval = eval_node_pot(node,integ_pts)	#
	try
		tilted_eval  = node_eval .* pdf(Normal(node_cavity[1],node_cavity[2]),integ_pts)
		tilted_node	 = params(fit_mle(Normal,integ_pts,tilted_eval))
		new_eta_node = normal_div(tilted_node,node_cavity)
		#
		tmp_mom = normal_prod(new_eta_node,node_cavity)
		if tmp_mom[2]>0.01
			q_moments[node,:] 		 = tmp_mom
			eta_node_moments[node,:] = [m for m in new_eta_node]
		end
	end
end

function ep_edge_proj(from,to)
	from_cavity = q_moments[from,:]
	to_cavity   = q_moments[to,:]
	eta_out     = get_edge_eta(from,to)
	eta_in 	    = get_edge_eta(to,from)
	#
	if eta_in[2]>from_cavity[2]
		from_cavity = normal_div(from_cavity,eta_in)
	end
	if eta_out[2]>to_cavity[2]
		to_cavity = normal_div(to_cavity,eta_out)
	end
	#
	from_cavity_d = Normal(from_cavity[1],from_cavity[2])
	to_cavity_d   = Normal(to_cavity[1],to_cavity[2])
	#
	eval_grid = zeros(Ninteg,Ninteg)
	for i=1:Ninteg
		xi = integ_pts[i]
		for j=1:Ninteg
			xj = integ_pts[j]
			eval_grid[i,j] = pdf(from_cavity_d,xi) *
								pdf(to_cavity_d,xj) *
									eval_edge_pot(from,to,xi,xj)
		end
	end
	#
	from_marg  = sum(eval_grid,2) # integrate out "to" (columns)
	from_marg /= sum(from_marg)
	to_marg    = sum(eval_grid,1) # integrate out "from" (rows)
	to_marg   /= sum(to_marg)
	#
	try
		from_marg_mom = params(fit_mle(Normal,integ_pts,from_marg))
		new_eta_in 	  = normal_div(from_marg_mom,from_cavity)
		tmp_mom 	  = normal_prod(new_eta_in,from_marg_mom)
		#
		if tmp_mom[2] > sigma_thresh
			q_moments[from,:] 		= tmp_mom
			edge_idx 		  		= get_edge_idx(to,from)
			eta_moments[edge_idx,:] = [m for m in new_eta_in]
		end
	end
	try
		to_marg_mom = params(fit_mle(Normal,integ_pts,to_marg))
		new_eta_out	= normal_div(to_marg_mom,to_cavity)
		tmp_mom 	= normal_prod(new_eta_out,to_marg_mom)
		#
		if tmp_mom[2] > sigma_thresh
			q_moments[to,:] 		= tmp_mom
			edge_idx 				= get_edge_idx(from,to)
			eta_moments[edge_idx,:] = [m for m in new_eta_out]
		end
	end
end
