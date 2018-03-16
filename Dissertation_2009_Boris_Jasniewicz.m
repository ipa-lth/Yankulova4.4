clear

# loading YALMIP
init_yalmip

# loading Fusions Reaktor
init_fusion_reactor

# select solver
sol = 'sdpt3'
#sol = 'sedumi'
#sol = 'mosek'

# Init values (may come from optimization loop)
pmin = 1/20
mu1 = 1.5
zeta0 = 2.5
zeta1 = 5

##############
# Löse A.15  #
##############
jasniewicz_step1

###############################
# Optimisation problem (4.40) #
###############################
jasniewicz_step2

disp("end")
