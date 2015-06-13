function doLBPD()
    _start_lbpd = time()
    println("LBPD sim ($expname::$Ngrid)")
    # > pre-allocation of storage space
    global messages = ones(2*nedges,Ngrid)
    global beliefs  = zeros(nnodes,Ngrid)
    # > initial beliefs: just node pot (mess=1)
    for node=1:nnodes
        neighbors   	= get_neighbors(node)
        cur_beliefs 	= eval_node_pot(node,grid)
        beliefs[node,:] = cur_beliefs/sum(cur_beliefs)
    end
    # > pre-computation in case of homogeneous edge pot
    global	edge_pot_grid = []
    if HOMOG_EDGE_POT
        edge_pot_grid = zeros(Ngrid,Ngrid)
        for gi=1:Ngrid
            for gj=1:Ngrid
                edge_pot_grid[gi,gj]=eval_edge_pot(0,0,grid[gi],grid[gj])
            end
        end
    end
    # > keep init belief, useful for later comp
    init_beliefs = copy(beliefs)
    #
    for loop = 1:nloops
        print(">loop: ",loop)
        _start_loop = time()
        for i=1:length(scheduling)
            lbpd_node_update(scheduling[i])
        end
        println(" [completed in ",get_time(_start_loop),"s]")
    end
    println("LBPD completed in ",get_time(_start_lbpd),"s.")
    #
    writecsv("$expname/$expname\_lbpd_origbel_np$Ngrid.dat",init_beliefs)
    writecsv("$expname/$expname\_lbpd_grid_np$Ngrid.dat",	grid)
    writecsv("$expname/$expname\_lbpd_beliefs_np$Ngrid.dat",beliefs)
end

function doEPBP()
    _start_epbp = time()
    println("EPBP sim ($expname::$N) [run::$R]")
    # > pre-allocation of storage space
    global particles   		= zeros(nnodes,N)
    global b_weights   		= zeros(nnodes,N)
    global b_evals   		= zeros(nnodes,N)
    global q_moments   		= zeros(nnodes,2)
    global e_weights   		= zeros(2*nedges,N)
    global eta_moments 		= zeros(2*nedges,2)
    global eta_node_moments = zeros(nnodes,2)
    #
    # > initial proposals & particles [!USER!]
    for node = 1:nnodes
        mu_node   		  = obs_values[node]
        q_moments[node,:] = [ mu_node s_init ]
        particles[node,:] = mu_node + s_init*randn(N,1)
    end
    #
    # > initial edge weights
    for edge = 1:2*nedges
        from     = edge_list[edge,1]
        weights  = eval_node_pot(from,particles[from,:])
        weights /= sum(weights)
        #
        e_weights[edge,:] = weights
    end
    #
    orig_q_moments   = copy(q_moments)
    orig_eta_moments = copy(eta_moments)
    #
    for loop = 1:nloops
        print(">loop: ",loop)
        _start_loop = time()
        for i=1:length(scheduling)
            epbp_node_update(scheduling[i])
        end
        println(" [completed in ",get_time(_start_loop),"s]")
    end
    println("EPBP completed in ",get_time(_start_epbp),"s.")
    print("...eval est. beliefs on mesh...")
    _start_epbp_estbel = time()
    epbp_est_beliefs   = zeros(nnodes,Ngrid)
    for node = 1:nnodes
        t1,t2  = epbp_eval_belief(node,grid)
        epbp_est_beliefs[node,:] = t2
    end
    println(" [done in ",get_time(_start_epbp_estbel),"s]")
    #
    writecsv("$expname/$expname\_epbp_est_beliefs_np$N\_r$R.dat",epbp_est_beliefs)
    writecsv("$expname/$expname\_epbp_particles_np$N\_r$R.dat",	 particles)
    writecsv("$expname/$expname\_epbp_weights_np$N\_r$R.dat",	 b_weights)
    writecsv("$expname/$expname\_epbp_evals_np$N\_r$R.dat",		 b_evals)
    writecsv("$expname/$expname\_epbp_qmom_np$N\_r$R.dat",		 q_moments)
end

function doFEPBP()
    _start_fepbp = time()
    println("FEPBP sim ($expname::$N/$C) [run::$R]")
    # > pre-allocation of storage space
    global particles   = zeros(nnodes,N)
    global b_weights   = zeros(nnodes,N)
    global b_evals     = zeros(nnodes,N)
    global q_moments   = zeros(nnodes,2)
    global e_weights   = zeros(2*nedges,N)
    global eta_moments = zeros(2*nedges,2)
    #
    # > initial proposals & particles [!USER!]
    for node = 1:nnodes
        mu_node   		  = obs_values[node]
        q_moments[node,:] = [ mu_node s_init ]
        particles[node,:] = mu_node + s_init*randn(N,1)
    end
    #
    # > initial edge weights
    for edge = 1:2*nedges
        from     = edge_list[edge,1]
        weights  = eval_node_pot(from,particles[from,:])
        weights /= sum(weights)
        #
        e_weights[edge,:] = weights
    end
    #
    orig_q_moments   = copy(q_moments)
    orig_eta_moments = copy(eta_moments)
    #
    for loop = 1:nloops
        print(">loop: ",loop)
        _start_loop = time()
        for i=1:length(scheduling)
            epbp_node_update(scheduling[i],true) # with fastmode
        end
        println(" [completed in ",get_time(_start_loop),"s]")
    end
    println("FEPBP completed in ",get_time(_start_fepbp),"s.")
    print("...eval est. beliefs on mesh...")
    _start_fepbp_estbel = time()
    fepbp_est_beliefs   = zeros(nnodes,Ngrid)
    for node = 1:nnodes
        t1,t2  = epbp_eval_belief(node,grid) # complete eval (not fast)
        fepbp_est_beliefs[node,:] = t2
    end
    println(" [done in ",get_time(_start_fepbp_estbel),"s]")
    #
    writecsv("$expname/$expname\_fepbp_est_beliefs_np$N\_nc$C\_r$R.dat",fepbp_est_beliefs)
    writecsv("$expname/$expname\_fepbp_particles_np$N\_nc$C\_r$R.dat",  particles)
    writecsv("$expname/$expname\_fepbp_weights_np$N\_nc$C\_r$R.dat",    b_weights)
    writecsv("$expname/$expname\_fepbp_evals_np$N\_nc$C\_r$R.dat",      b_evals)
    writecsv("$expname/$expname\_fepbp_qmom_np$N\_nc$C\_r$R.dat",       q_moments)
end

function doPBP()
    _start_pbp = time()
    println("PBP sim ($expname::$N) [run::$R]")
    #
    # > pre-allocation of storage space
    global particles   = zeros(nnodes,N)
    global b_evals     = zeros(nnodes,N)
    global messages    = ones(2*nedges,N)
    #
    for node=1:nnodes
        node_p 			  = obs_values[node]+s_init*randn(1,N)
        bel_p  			  = eval_node_pot(node,node_p)
        particles[node,:] = node_p
        b_evals[node,:]   = bel_p/sum(bel_p)
    end
    #
    for loop=1:nloops
        print(">loop: ",loop)
        _start_loop = time()
        for i=1:length(scheduling)
            pbp_node_update(scheduling[i]) # with fastmode
        end
        println(" [completed in ",get_time(_start_loop),"s]")
    end
    #
    println("PBP completed in ",get_time(_start_pbp),"s.")
    print("...eval est. beliefs on mesh...")
    _start_pbp_estbel = time()
    pbp_est_beliefs   = zeros(nnodes,Ngrid)
    for node = 1:nnodes
        pbp_est_beliefs[node,:] = pbp_eval_belief(node,grid)
    end
    println(" [done in ",get_time(_start_pbp_estbel),"s]")
    #
    writecsv("$expname/$expname\_pbp_est_beliefs_np$N\_r$R.dat",pbp_est_beliefs)
    writecsv("$expname/$expname\_pbp_particles_np$N\_r$R.dat",  particles)
    writecsv("$expname/$expname\_pbp_evals_np$N\_r$R.dat",      b_evals)
end

function doEP()
    _start_ep = time()
    println("EP sim ...")
    #
    global eta_moments 		= zeros(2*nedges,2)
    global eta_node_moments	= zeros(nnodes,2)
    global q_moments   		= zeros(nnodes,2)
    #
    for node=1:nnodes
        q_moments[node,:] = [obs_values[node] s_init]
    end
    for loop=1:nEPloops
        print(">loop ",loop)
        _start_loop = time()
        old_moms 	= copy(q_moments)
        for i=1:length(scheduling)
            ep_node_update(scheduling[i])
        end
        println(" [completed in ",get_time(_start_loop),"s; raw diff ",
                    round(norm(q_moments-old_moms),3),"]")
    end
    #
    println("EP completed in ",get_time(_start_ep),"s.")
end
