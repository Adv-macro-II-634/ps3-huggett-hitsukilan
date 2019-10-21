%%%%%HW3 Huggett's model %%%%%
clear;
clc;

% PARAMETERS
beta = .994; % discount factor 
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


aggsav = 1 ;

while (abs(aggsav) >= 0.01)
    
    q_guess = (q_min + q_max) / 2;
    
    % CURRENT RETURN (UTILITY) FUNCTION
    
    cons = bsxfun(@minus, a', q_guess * a);
    
    cons = bsxfun(@plus, cons, permute(y_s, [1 3 2]));
    
    ret = (cons .^ (1-sigma)) ./ (1 - sigma); % current period utility
    ret(cons<0) = -Inf;
    
    % INITIAL VALUE FUNCTION GUESS
    
    v_guess = zeros(2, num_a);
    
    % VALUE FUNCTION ITERATION
    v_tol = 1;
    
    while v_tol >.000001
        
        % CONSTRUCT RETURN + EXPECTED CONTINUATION VALUE
        value_mat=ret+beta*repmat(permute((PI*v_guess),[3 2 1]), [num_a 1 1]);

        % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
       [vfn, pol_indx] = max(value_mat, [], 2); %max for each row
        
       v_tol = max(abs(permute(vfn, [3 1 2])-v_guess));
       v_tol = max(v_tol(:));
       
       v_guess = permute(vfn, [3 1 2]);
    end
    
    % KEEP DECSISION RULE
    pol_indx = permute(pol_indx, [3 1 2]);
    pol_fn = a(pol_indx);
    
    % SET UP INITITAL DISTRIBUTION
    mu = zeros(2,num_a);
    mu(:) = 1/(2*num_a);
    
    dis=1;
  while dis>0.0000001 
      % ITERATE OVER DISTRIBUTIONS
      MuNew = zeros(size(mu));
     [emp_ind, a_ind, mass] = find(mu); % find non-zero indices
    
    for ii = 1:length(emp_ind)
        apr_ind = pol_indx(emp_ind(ii),a_ind(ii)); % which a prime does the policy fn prescribe?
        
        MuNew(:, apr_ind) = MuNew(:, apr_ind) + ... % which mass of households goes to which exogenous state?
            (PI(emp_ind(ii), :)*mass(ii))';
        
    end
    dis = max(max(abs(mu-MuNew)));
    mu = MuNew;
  end
   %Market clears
   aggsav = mu(1,:)*pol_fn(1,:)'+mu(2,:)*pol_fn(2,:)';
   if aggsav>0
       q_min = q_guess;
   else
       q_max = q_guess;
   end
end

% Question 3 Gini coefficient and lorenz curve

popu = reshape(mu', [2 * num_a, 1]);
wealth = reshape([a+y_s(1);a+y_s(2)]',[2 * num_a, 1]);
earning = reshape([repmat(y_s(1),[1, num_a]); repmat(y_s(2),[1, num_a])]', [2 * num_a,1]);

nonnegwealth = wealth;
nonnegwealth(nonnegwealth<0) = 0;

figure(1)
subplot(1,2,1)
gini_wealth = gini(popu, nonnegwealth, true);
title({'Lorenz Curve for Wealth','Gini coefficient = ',num2str(gini_wealth)});

subplot(1,2,2)
gini_earnings = gini(popu, earning, true);
title({'Lorenze Curve for Earning','Gini coefficient = ',num2str(gini_earnings)});

% Extra Credit

PI_LR = PI^1000;
C_bar = PI_LR(1,:) * y_s';
W_FB = (C_bar^(1-sigma))/(1-sigma)/(1-beta);
lambda = (v_guess.^(-1).* W_FB).^(1/(1-sigma))-1;
Hh = sum(sum((lambda >0) .* mu));
WFG = sum(sum(lambda .* mu));

figure(2)

plot(a,lambda(1,:),'DisplayName','Employed', 'LineWidth',2);
hold on;
plot(a,lambda(2,:),'DisplayName','Unemployed','LineWidth',2);
hold off;
title('Consumption Equivalent','FontSize',18);
lgd = legend;
lgd.Location ='Northeast';
lgd.FontSize = 14;






