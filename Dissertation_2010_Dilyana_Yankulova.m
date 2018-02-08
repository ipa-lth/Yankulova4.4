# loading YALMIP
init_yalmip

# loading Hydraulisher Aktor
init_hyd_actor


# Parameter
p_min = 0.1;

# Define variables
Q = sdpvar(n, n);  # Q = R^(-1)
z = sdpvar(n, 1); 

# Constants
N   = diag(-n:-1); # N = diag(-n,..., -2,-1)
M   =  diag(0: n-1); # M = diag(0,1 ..., n-1); 

# Constraints

constraints_451 = [Q >= 0]; # 4.59

constraints_451 = [constraints_451,
                   Q*(A' + a*b') + (A + b*a')*Q - z*b' - b*z' <= 0]; #4.60

constraints_451 = [constraints_451, 
                   Q*N + N*Q <= 0]; # 4.61

for i = 1: X0
  constraints_451 = [constraints_451,
                    [[1,          X0(:,i)'];
                     [X0(:,i),    Q       ]] > 0]; # 4.62
end

constraints_451 = [constraints_451,
                  [[u_max^2-a'*Q*a + 2*a'*z, z']; # no blank at minus or it's like two cols...
                   [z                      , Q ]] >=0]; # 4.63

                   
###########################################
# Constraint variant (4.64)               #
###########################################

# page 50: m <= 2*n -1
m = n;

for i = 0 : m
  constraints_451 = [constraints_451,
                     func_constraint464(i, Q, a, M, N, n, z) >= 0];
end

###########################################
# Objective                               #
###########################################

%objective = -logdet(Q);
objective = -geomean(Q);


###########################################
# Solve/Optimize                          #
###########################################
sol = optimize(constraints_451, objective)

disp("Status:")
disp(sol.info)
if sol.problem==0
  display('Feasible');
else
  display('Infeasible');
end

Q_451     = value(Q)
R1_451    = inv(Q_451)
a_hat_451 = R1_451 * value(z)

disp("z = ")
disp(value(z))
