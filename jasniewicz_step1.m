######################
# Löse A.15 mit mu=1 #
######################

# Define Variables
Q = sdpvar(n, n);  # Q = R^(-1)
#Q = sdpvar(n, n, 'full');  # Q = R^(-1)

z0 = sdpvar(n, 1);
z1 = sdpvar(n, 1);


# Bisection parameter
#g = cvxpy.Parameter(sign='positive')
sdpvar t;

mu = 1;
#mu = cvxpy.Parameter(sign='positive')
          
# Define Constraints
constraints_A15 = [Q >= 0]; # Q is Semidefinite

# (A.10)
for i = 1:size(X0,2)
  constraints_A15 = [constraints_A15,
                    [Q,       [X0(:,i)];
                    [X0(:,i)',       1]] >= 0];
end

# (A.11)
constraints_A15 = [constraints_A15,
                  [Q,   z0;
                   z0',  1 ] >= 0];

# (A.12
constraints_A15 = [constraints_A15,
                  [Q,   z1;
                   z1',  mu**2 ] >= 0];

# (A.13)
constraints_A15 = [constraints_A15,
                  Q*A+A'*Q-b*z0'-z0*b'<=0];

# (A.14)
constraints_A15 = [constraints_A15,
                  Q*A + A'*Q - b*z1' - z1*b' + 2*t*Q<=0];


ops = sdpsettings('solver','bisection','bisection.solver', sol);
Objective = -t; #Maximize
diagnostics = optimize(constraints_A15, Objective, ops);
#test_it(value(Q), value(z0), value(z1), mu, value(t), A, b, X0)
#check(constraints_A15)

disp("l*(mu=1)=")
l1m0 = inv(value(Q))*value(z1)

disp("lambda^{hat}(p=1)=")
roots_m0_p1 = eig(A-b*l1m0')'



#######################
# Löse A.15 mit mu>=1 #
#######################

# Define Variables
Q = sdpvar(n, n);  # Q = R^(-1)
#Q = sdpvar(n, n, 'full');  # Q = R^(-1)

z0 = sdpvar(n, 1);
z1 = sdpvar(n, 1);


# Bisection parameter
#g = cvxpy.Parameter(sign='positive')
sdpvar t;

mu = mu1;
#mu = cvxpy.Parameter(sign='positive')
          
# Define Constraints
constraints_A15 = [Q >= 0]; # Q is Semidefinite

# (A.10)
for i = 1:size(X0,2)
  constraints_A15 = [constraints_A15,
                    [Q,       [X0(:,i)];
                    [X0(:,i)',       1]] >= 0];
end

# (A.11)
constraints_A15 = [constraints_A15,
                  [Q,   z0;
                   z0',  1 ] >= 0];

# (A.12
constraints_A15 = [constraints_A15,
                  [Q,   z1;
                   z1',  mu**2 ] >= 0];

# (A.13)
constraints_A15 = [constraints_A15,
                  Q*A+A'*Q-b*z0'-z0*b'<=0];

# (A.14)
constraints_A15 = [constraints_A15,
                  Q*A + A'*Q - b*z1' - z1*b' + 2*t*Q<=0];


ops = sdpsettings('solver','bisection','bisection.solver', sol);
Objective = -t; #Maximize
diagnostics = optimize(constraints_A15, Objective, ops);
#test_it(value(Q), value(z0), value(z1), mu, value(t), A, b, X0)
#check(constraints_A15)

disp("l*(mu=1)=")
l1m1 = inv(value(Q))*value(z1)

disp("lambda^{hat}(p=1)=")
roots_m1_p1 = eig(A-b*l1m1')'

disp("lambda^{hat}(p=pmin)=")
roots_m0_pmin = zeta0 * real(roots_m0_p1) + 1i*imag(roots_m0_p1)
#roots_m0_pmin = zeta0 * roots_m0_p1  # Shifting on both real and imaginary axis in zeta

disp("lambda^{hat}_*(p=pmin)=")
roots_m1_pmin = zeta1 * real(roots_m1_p1) + 1i*imag(roots_m1_p1)
#roots_m1_pmin = zeta0 * roots_m1_p1  # Shifting on both real and imaginary axis in zeta

[k0_0, k0_1] = k_explizit_Ab2(roots_m0_p1, roots_m0_pmin, pmin, A, b)
figure;
plot_moving_poles(A, b, c, d, k0_0, k0_1, pmin)

[k1_0, k1_1] = k_explizit_Ab2(roots_m1_p1, roots_m1_pmin, pmin, A, b)
figure;
plot_moving_poles(A, b, c, d, k1_0, k1_1, pmin)
