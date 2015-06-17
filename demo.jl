#
# 	Code released under the MIT license, see corresponding LICENSE file
#	(c) 2015, Thibaut Lienart
#
# --------------------------------------------------------------------------------------------------
#
include("lib_support.jl")
include("lib_epbp.jl")
include("lib_lbpd.jl")
include("lib_pbp.jl")
include("lib_ep.jl")
include("lib_doSim.jl")
#
demoNames = ["demoGrid","demoTree","demoImg","demoChain"]
expname   = demoNames[2] # choice of demo (reason for this syntax is to show possibilities)
# make directory to store stuff
if ~isdir(expname)
    mkdir(expname)
end
#
# ==================================================================================================
if expname == "demoGrid"
    #
    # SIMULATIONS TO BE RUN
    #
    RELOAD = false  # re-generate everything
    LBPD   = false  # LBP on deterministic grid
    EPBP   = false  # EPBP
    FEPBP  = false  # Fast-EPBP
    PBP    = true  # PBP with MH sampling
    EP 	   = false   # straight EP
    #
    # VERBOSITY
    NODE_PROGRESS = false
    #
    # SIMULATION PARAMETERS [!USER!]
    #
    Nlist	  = [50] # (list) number of particles per node
    Clist 	  = [7]    # (list) number of components for FEPBP, need to be of same dim as NLIST
    Ninteg    = 30		  # number of integration points for EP proj
    Ngrid     = 200		  # number of points in the discretization
    nloops    = 10 		  # number of loops through scheduling
    nEPloops  = 25 		  # number of EP iterations
    nruns     = 1  		  # number of time we run the whole thing
    #
    # EP PROJECTION MODE, default is KL ignoring update if moments not valid
    EP_PROJ_MLE  = false    # use MLE projection instead of KL-EP
    #
    # Additional parameters for PBP
    #
    MHIter 	   = 20 		  	# number of MH iterations
    MHProposal = Normal(0,1) 	# form of the MH proposal
    PARACHAINS = true
    #
    # DECLARE GM
    #
    # > Declare a 5x5 grid
    m,n = 5,5
    nnodes,nedges,edge_list = gm_grid(m,n)
    # > declare scheduling
    scheduling = gm_grid_scheduling(m,n)
    # > declare edge and node potential
    HOMOG_EDGE_POT = true # edge pot is symmetric & same everywhere (eg: image)
    # > side functions to define edge/node potential
    node_potential = MixtureModel([Normal(-2,1),Gumbel(2,1.3)],[0.6,0.4])
    edge_potential = Laplace(0,2)
    # > definition of edge/node potential functions
    eval_edge_pot(from,to,xfrom,xto) = pdf(edge_potential,xfrom-xto)
    eval_node_pot(node,xnode)        = pdf(node_potential,obs_values[node]-xnode)
    # > (PBP) sampling from MH?
    sampleMHP(old) = old+rand(MHProposal,N)'
    LMHCHAIN = 1000
    ENDCHUNK = 300
    sampleMHP2(old) = old+rand(Normal(0,0.1),1)
    # > initial values on the graph
    orig_values = zeros(nnodes,1) + 2
    #
    est_range    = (-5,15) 	# > estimated 1D-range for integration
    sigma_thresh = 0.01
    if RELOAD
    	# > generate observations
    	obs_values = orig_values + rand(node_potential,nnodes)
    	writecsv("$expname/$expname\_orig_values.dat",orig_values)
    	writecsv("$expname/$expname\_obs_values.dat",obs_values)
    	#
    	obs_var = sqrt(var(obs_values))
    	s_init  = 4*obs_var
    end
# ==================================================================================================
elseif expname == "demoTree"
    #
    # SIMULATIONS TO BE RUN
    #
    RELOAD = true  # re-generate everything
    LBPD   = true  # LBP on deterministic grid
    EPBP   = true  # EPBP
    FEPBP  = false  # Fast-EPBP
    PBP    = true  # PBP with MH sampling
    EP 	   = false   # straight EP
    #
    # VERBOSITY
    NODE_PROGRESS = false
    #
    # SIMULATION PARAMETERS [!USER!]
    #
    Nlist	  = [100] # (list) number of particles per node
    Clist 	  = [7]    # (list) number of components for FEPBP, need to be of same dim as NLIST
    Ninteg    = 30		  # number of integration points for EP proj
    Ngrid     = 200		  # number of points in the discretization
    nloops    = 10 		  # number of loops through scheduling
    nEPloops  = 25 		  # number of EP iterations
    nruns     = 1  		  # number of time we run the whole thing
    #
    # EP PROJECTION MODE, default is KL ignoring update if moments not valid
    EP_PROJ_MLE  = false    # use MLE projection instead of KL-EP
    #
    # Additional parameters for PBP
    #
    MHIter 	   = 20 		  	# number of MH iterations
    MHProposal = Normal(0,1) 	# form of the MH proposal
    PARACHAINS = true
    #
    # DECLARE GM
    #
    nnodes,nedges = 8,7
    edge_list = [   1 2; 1 3;                           # top 1-{2,3}
                    2 1; 2 4; 2 5; 2 6; 3 1; 3 7; 3 8;  # middle 2-{1,4,5,6} and 3-{1,7,8}
                    4 2; 5 2; 6 2; 7 3; 8 3; ]          # leaves 4-2,5-2,6-2,7-3,8-3
    sched1     = hcat(1,2,4,5,6,3,7,8)
    sched2     = hcat(4,5,6,2,7,8,3,1)
    scheduling = vcat(sched1[:],sched2[:])
    # > declare edge and node potential
    HOMOG_EDGE_POT = true # edge pot is symmetric & same everywhere (eg: image)
    # > side functions to define edge/node potential
    node_potential = MixtureModel([Normal(-2,1),Gumbel(1,1.5)],[0.7,0.3])
    edge_potential = Cauchy(0,1)
    # > definition of edge/node potential functions
    eval_edge_pot(from,to,xfrom,xto) = pdf(edge_potential,xfrom-xto)
    eval_node_pot(node,xnode)        = pdf(node_potential,obs_values[node]-xnode)
    # > (PBP) sampling from MH?
    sampleMHP(old) = old+rand(MHProposal,N)'
    LMHCHAIN = 1000
    ENDCHUNK = 300
    sampleMHP2(old) = old+rand(Normal(0,0.1),1)
    # > initial values on the graph
    orig_values = zeros(nnodes,1) + 2
    #
    est_range    = (-7,7) 	# > estimated 1D-range for integration
    sigma_thresh = 0.01
    if RELOAD
    	# > generate observations
    	obs_values = orig_values + rand(node_potential,nnodes)
    	writecsv("$expname/$expname\_orig_values.dat",orig_values)
    	writecsv("$expname/$expname\_obs_values.dat",obs_values)
    	#
    	obs_var = sqrt(var(obs_values))
    	s_init  = 4*obs_var
    end
# ==================================================================================================
elseif expname == "demoImg"
    #
    # SIMULATIONS TO BE RUN
    #
    RELOAD = false  # re-generate everything
    LBPD   = false  # LBP on deterministic grid
    EPBP   = false  # EPBP
    FEPBP  = true   # Fast-EPBP
    PBP    = false  # PBP with MH sampling
    EP 	   = false  # straight EP
    #
    # VERBOSITY
    NODE_PROGRESS = false
    #
    # SIMULATION PARAMETERS [!USER!]
    #
    Nlist	  = [70]	    # (list) number of particles per node
    Clist 	  = [5]		    # (list) number of components for FEPBP, need to be of same dim as NLIST
    Ninteg    = 50			# number of integration points for EP proj
    Ngrid     = 200			# number of points in the discretization
    nloops    = 2 			# number of loops through scheduling
    nEPloops  = 3 			# number of EP iterations
    nruns     = 1  			# number of time we run the whole thing
    #
    # EP PROJECTION MODE, default is KL ignoring update if moments not valid.
    EP_PROJ_MLE  = false    # use MLE projection instead of KL-EP
    #
    # Additional parameters for PBP
    #
    MHIter 	   = 20 		  	# number of MH iterations
    MHProposal = Normal(0,.1) 	# form of the MH proposal
    PARACHAINS = true
    #
    # DECLARE GM
    #
    obs        = readdlm("ex_squareNoisy.dat")
    obs_values = obs[:]
    obs_var    = sqrt(var(obs_values))
    s_init     = 4*obs_var
    #
    est_range = (-1,1.5)
    sigma_thresh = 0.01
    #
    m,n = 50,50
    nnodes,nedges,edge_list = gm_grid(m,n)
    # > declare scheduling
    scheduling = gm_grid_scheduling(m,n)
    #
    HOMOG_EDGE_POT = true
    node_potential = Normal(0,0.1)
    #
    function eval_edge_pot(from,to,xfrom,xto)
        d = abs(xfrom-xto)
        if length(d)==1
            r = d
            if d>0.25
                v = exp(-0.25/0.035)
                r = d*0+v
            else
                r = exp(-d/0.035)
            end
        else
            r=copy(d)
            for i =1:length(d)
                if d[i]>0.25
                    v = exp(-0.25/0.035)
                    r[i] = d[i]*0+v
                else
                    r[i] = exp(-d[i]/0.035)
                end
            end
        end
        return r
    end
    #
    function eval_node_pot(node,xnode)
        return pdf(node_potential,obs_values[node]-xnode)
    end
# ==================================================================================================
elseif expname == "demoChain"
    #
    # SIMULATIONS TO BE RUN
    #
    RELOAD = true  # re-generate everything
    LBPD   = true  # LBP on deterministic grid
    EPBP   = true  # EPBP
    FEPBP  = false  # Fast-EPBP
    PBP    = false  # PBP with MH sampling
    EP 	   = true   # straight EP
    #
    # VERBOSITY
    NODE_PROGRESS = false
    #
    # SIMULATION PARAMETERS [!USER!]
    #
    Nlist	  = [100,200] # (list) number of particles per node
    Clist 	  = [7,10]    # (list) number of components for FEPBP, need to be of same dim as NLIST
    Ninteg    = 30		  # number of integration points for EP proj
    Ngrid     = 200		  # number of points in the discretization
    nloops    = 10 		  # number of loops through scheduling
    nEPloops  = 25 		  # number of EP iterations
    nruns     = 1  		  # number of time we run the whole thing
    #
    # EP PROJECTION MODE, default is KL ignoring update if moments not valid.
    EP_PROJ_MLE  = false    # use MLE projection instead of KL-EP
    #
    # Additional parameters for PBP
    #
    MHIter 	   = 20 		  	# number of MH iterations
    MHProposal = Normal(0,.1) 	# form of the MH proposal
    PARACHAINS = true
    #
    # DECLARE GM
    #
    T = 5
    nnodes,nedges,edge_list = gm_chain(T)
    # > declare scheduling
    scheduling = gm_chain_scheduling(T,true) # forward only
    # > declare edge and node potential
    HOMOG_EDGE_POT = false # if edge pot is symmetric & the same everywhere (eg: image)
    #
    node_noise = Normal(0,2)
    node_mult  = 0.7
    edge_noise = Normal(0,1)
    edge_mult  = 0.5
    #
    eval_edge_pot(from,to,xfrom,xto) = pdf(edge_noise,xto-edge_mult*xfrom)
    eval_node_pot(node,xnode)        = pdf(node_noise,obs_values[node]-node_mult*xnode)
    #
    # > (PBP) sampling from MH?
    sampleMHP(old) = old+rand(MHProposal,N)'
    LMHCHAIN = 1000
    ENDCHUNK = 300
    sampleMHP2(old) = old+rand(Normal(0,0.1),1)
    #
    est_range    = (-10,10) 	# > estimated 1D-range for integration
    sigma_thresh = 0.01
    # > generate observations
    if RELOAD
    	orig_values = zeros(nnodes,1)
    	obs_values  = zeros(nnodes,1)
    	#
    	orig_values[1] = rand(Normal(0,1))
    	obs_values[1]  = rand(Normal(node_mult*orig_values[1],1))
    	#
    	for node=2:T
    		orig_values[node] = rand(Normal(edge_mult*orig_values[node-1],2))
    		obs_values[node]  = rand(Normal(node_mult*orig_values[node],1))
    	end
    	writecsv("$expname/$expname\_orig_values.dat",orig_values)
    	writecsv("$expname/$expname\_obs_values.dat",obs_values)
    	#
    	# > to start
    	obs_var = sqrt(var(obs_values))
    	s_init  = 4*obs_var
    end
end

#
# ======== RUN SIMULATIONS =========================================================================
#
integ_pts = linspace(est_range[1],est_range[2],Ninteg)' # ! leave the transpose
grid      = linspace(est_range[1],est_range[2],Ngrid)'
# > generate observations
#
# --------------------------------------------------------------------------------------------------
# (cf. lib_doSim.jl)
#
if LBPD
	doLBPD()
end

for N_index in 1:length(Nlist)
	global N = Nlist[N_index]
	global C = Clist[N_index] # for FEPBP
    global R
	for R = 1:nruns
		if EPBP
			doEPBP()
		end
		if FEPBP
			doFEPBP()
		end
		if PBP
			doPBP()
		end
	end
end

if EP
    doEP()
end
