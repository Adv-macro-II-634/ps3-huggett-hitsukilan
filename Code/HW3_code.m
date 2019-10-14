%%%%%HW3 Huggett's model %%%%%
clear;
clc;

% PARAMETERS
beta = .9932; % discount factor 
sigma = 1.5; % coefficient of risk aversion
b = 0.5; % replacement ratio (unemployment benefits)
y_s = [1, b]; % endowment in employment states
PI = [.97 .03; .5 .5]; % transition matrix


% ASSET VECTOR
a_lo = -2; %lower bound of grid points
a_hi = 5;% upper bound of grid points
num_a = 1000;

a = linspace(a_lo, a_hi, num_a); % asset (row) vector

% INITIAL GUESS FOR q
q_min = 0.98;
q_max = 1.2;
q_guess = (q_min + q_max) / 2;

% ITERATE OVER ASSET PRICES
aggsav = 1 ;
while abs(aggsav) >= 0.01 
    
    % CURRENT RETURN (UTILITY) FUNCTION
    cons = bsxfun(@minus, a', q_guess * a);
    cons = bsxfun(@plus, cons, permute(y_s, [1 3 2]));
    ret = (cons .^ (1-sigma)) ./ (1 - sigma); % current period utility
    
    % INITIAL VALUE FUNCTION GUESS
    v_guess = zeros(2, num_a);
    
    % VALUE FUNCTION ITERATION
    v_tol = 1;
    
    while v_tol >.0001
        % CONSTRUCT RETURN + EXPECTED CONTINUATION VALUE
        val_mat = ret + beta * repmat(permute((PI*v_guess),[3 2 1]), [num_a 1 1]);
        % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
        [vfn, pol_index] = max(val_mat, [], 2);
        
        v_tol = max(abs(permute(vfn, [3 1 2]) - v_guess));
        v_tol = max(v_tol(:));
       
        v_guess = permute(vfn, [3 1 2]);
  
    end
    
    % KEEP DECSISION RULE
    
    pol_indx = permute(pol_index,[3 1 2]);
    pol_fn = a(pol_indx);
    
    % SET UP INITITAL DISTRIBUTION
    
    mu = zerox(2,num_a);
    mu(:) = 1/(2 * num_a);
    
    dis = 1;
    while dis > 0.00001
        
    
        % ITERATE OVER DISTRIBUTIONS
        MuNew = zeros(size(Mu));
        [emp_ind, a_ind, mass] = find(Mu > 0); % find non-zero indices
       
        for ii = 1:length(emp_ind)
            
        apr_ind = pol_indx(emp_ind(ii), a_ind(ii)); % which a prime does the policy fn prescribe?
        
        MuNew(:, apr_ind) = MuNew(:, apr_ind) + ... % which mass of households goes to which exogenous state?
            (PI(emp_ind(ii), :) * mass)';
        end
        
     dis = max(max(abs(mu-MuNew)));
     mu = MuNew;
    end
    % market clearing 
    
    aggsav = mu(1,:) * pol_fn(1,:)' + mu(2,:) * pol_fn(2,:)';
    
    if aggsav > 0
        q_min = q_guess;
    else
        q_max = q_guess;
    end
        
end

% Question 2 figure

figure(1)
subplot(1,2,1)
plot(a,v_guess(1,:),'DisplayName','employed', 'LineWidth',2);
hold on;
plot(a,v_guess(2,:),'DisplayName','unemployed','LineWidth',2);
hold off;
title('Value Function','FontSize',18);
lgd = legend;
lgd.Location ='Southeast';
lgd.FontSize = 14;

subplot(1,2,2)
plot(a,pol_fn(1,:),'DisplayName','employed', 'LineWidth',2);
plot(a,pol_fn(2,:),'DisplayName','unemployed', 'LineWidth',2);
title('Policy Function','FontSize',18);

% Question 3 Gini coefficient and lorenz curve


