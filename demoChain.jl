EP_PROJ_MLE  = false
DEB_MIX		 = false
#
# Accomodates: 1D, Normal approx (2 moments)
# > should be flexible wrt these params in future
#
# things that should be determined by the user in this code
# > everything related to declaration of graphical model
# > how to initialize the proposals
# >> check [!USER!]
#
include("lib_support.jl")
include("lib_epbp.jl")
include("lib_lbpd.jl")
include("lib_pbp.jl")
include("lib_ep.jl")
include("lib_doSim.jl")
#
RELOAD = false  # re-generate everything
LBPD   = false  # LBP on deterministic grid
EPBP   = true   # EPBP
FEPBP  = false  # Fast-EPBP
PBP    = false  # PBP with MH sampling
EP 	   = false  # straight EP
#
expname = "demoChain"
#
# SIMULATION PARAMETERS [!USER!]
#
Nlist	 = [500]		    # number of particles per node
Clist 	 = [10]				# number of components for FEPBP
Ninteg   = 30				# number of integration points for EP proj
Ngrid    = 200				# number of points in the discretization
nloops   = 25 				# number of loops through scheduling
nEPloops = 30 				# number of EP iterations
nruns    = 1  				# number of time we run the whole thing
#
MHIter 	   = 20 		  	# number of MH iterations
MHProposal = Normal(0,.1) 	# form of the MH proposal
#
est_range = (-10,10) 		# > estimated range for integration
#
# DECLARE GRAPHICAL MODEL [!USER!]
#
# > declare underlying structure (cf. LIB_SUPPORT for eg: Grid)
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
# > sampling from MH?
sampleMHP(old) = old+rand(MHProposal,N)'
#
# make directory to store stuff
if ~isdir(expname)
    mkdir(expname)
end
#
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
#
# ==================================================================================================
# ======== RUN SIMULATIONS =========================================================================
# ==================================================================================================
#
integ_pts = linspace(est_range[1],est_range[2],Ninteg)' # ! leave the transpose
grid      = linspace(est_range[1],est_range[2],Ngrid)'
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
