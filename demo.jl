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
expname = "demoGrid" # results will be stored in $expname/
#
# SIMULATIONS TO BE RUN
#
RELOAD = true   # re-generate everything
LBPD   = true   # LBP on deterministic grid
EPBP   = true   # EPBP
FEPBP  = false  # Fast-EPBP
PBP    = false  # PBP with MH sampling
EP 	   = false  # straight EP
#
# EP PROJECTION MODE, default is KL ignoring update if moments not valid.
#
EP_PROJ_MLE  = false    # use MLE projection instead of KL-EP
#
# SIMULATION PARAMETERS [!USER!]
#
Nlist	  = [50,100]	    # number of particles per node
Clist 	  = [7,10]		# number of components for FEPBP
Ninteg    = 30			# number of integration points for EP proj
Ngrid     = 200			# number of points in the discretization
nloops    = 10 			# number of loops through scheduling
nEPloops  = 30 			# number of EP iterations
nruns     = 1  			# number of time we run the whole thing
est_range = (-5,15) 	# > estimated 1D-range for integration
#
# Additional parameters for PBP
#
MHIter 	   = 20 		  	# number of MH iterations
MHProposal = Normal(0,.1) 	# form of the MH proposal
#
# DECLARE GRAPHICAL MODEL [!USER!]
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
# > initial values on the graph
orig_values = zeros(nnodes,1) + 2
#
# ======== RUN SIMULATIONS =========================================================================
#
# make directory to store stuff
if ~isdir(expname)
    mkdir(expname)
end
#
integ_pts = linspace(est_range[1],est_range[2],Ninteg)' # ! leave the transpose
grid      = linspace(est_range[1],est_range[2],Ngrid)'
# > generate observations
if RELOAD
	# > generate observations
	obs_values = orig_values + rand(node_potential,nnodes)
	writecsv("$expname/$expname\_orig_values.dat",orig_values)
	writecsv("$expname/$expname\_obs_values.dat",obs_values)
	#
	obs_var = sqrt(var(obs_values))
	s_init  = 4*obs_var
end
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
