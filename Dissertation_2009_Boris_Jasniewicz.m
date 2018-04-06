clear

# loading YALMIP
init_yalmip

# loading Fusions Reaktor
init_fusion_reactor
#init_hyd_actor

# select solver
sol = 'sdpt3'
#sol = 'sedumi'
#sol = 'mosek'

# Init values (may come from optimization loop)

#pmin = 1.0/20.0
pmin = 0.02
mu1 = 1.5
zeta0 = 2.5
zeta1 = 5.0


##############
# Löse A.15  #
##############
jasniewicz_step1

###############################
# Optimisation problem (4.40) #
###############################
#jasniewicz_step2

instruct.n = n;
instruct.p_min = pmin;
instruct.M_l = -1;
instruct.M_u = 0;
instruct.N_l = -1;
instruct.N_u = 0;
instruct.bet = 7.5; # ??
instruct.A = A;
instruct.B = b;
instruct.x0 = X0; # ??
instruct.k = {k0_1, k0_0};
instruct.ks = {k1_1, k1_0};

[R_sol, sol, res] = solveLMIs(instruct);

if sol.problem ~= 0 % numerical problems?
    f = realmax;%-1/res;%
    disp(sol.info)
    return
end

disp(sol.info)
