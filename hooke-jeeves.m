% Name - Utkarsh Yashvardhan
% BITS ID - 2019B4A70704P

% Problem statement - Implement the Hooke and Jeeves method in python/matlab to minimize the function
% over R2 with starting point (1, 0).

% The given program computes the minimum point of a multi-dimensional function using the method of Hooke-Jeeves using discrete steps.
% 
% Given below is a sample run to determine the optimum value for the function given below:
% 
% y = (X(1) - 2)^4 + (X(1) - 2*X(2))^2
% 
% The search for the optimum value of the 2 variables has the initial estimate of [1 0], 
% initial step size of [0.1 0.1] with a minimum step size of [1e-9 1e-9]. 
% The search has upper bound of 1000 over the number of iterations and a function tolerance limit of 1e-9:
% 
% >> [X, optimal_value, no_of_iterations] = hooke_jeeves_method(2, [0 3], [.1 .1], [1e-9 1e-9], 1e-9, 1000, 'function1')
% 
% X =
% 
%     1.9714    0.9857
% 
% 
% optimal_value =
% 
%    6.6847e-07
% 
% 
% no_of_iterations =
% 
%    335


% -------------------------------------------------------------------------
% Just run the given code on MATLAB (with the file having .m extension) and 
% find the optimal value of the following function (function1) along with 
% optimal points (and the number of iterations it took to find optimal objective function value).
% -------------------------------------------------------------------------
[X, optimal_value, no_of_iterations] = hooke_jeeves_method(2, [1 0], [1e-3 1e-3], [1e-9 1e-9], 1e-9, 1000, 'function1')
% [X, optimal_value, no_of_iterations] = hooke_jeeves_method(2, [1 0], [.1 .1], [1e-6 1e-6], 1e-7, 500, 'function2')



% True optimal value of this function is 0 and it occurs at X(1) = 2 and X(2) =
% 1
function y = function1(X, ~)
    y = (X(1) - 2)^4 + (X(1) - 2*X(2))^2;
end

% True optimal value of this function is 10 and it occurs at X(1) = 2 and X(2) =
% -5
% function y = function2(X, ~)
%     y = 10 + (X(1) - 2)^2 + (X(2) + 5)^2;
% end

function[X, optimal_value, no_of_iterations] = hooke_jeeves_method(no_of_variables, X, step_size, min_step_size, tolerance, iterations_upper_bound, objective_function)

    % 'hooke_jeeves_method()' performs multi-dimensional optimization using the
    % method of Hooke-Jeeves with discrete steps.
    %
    % Input:
    %
    % no_of_variables - number of variables
    % X - array of initial estimates
    % step_size - search step sizes for each dimension
    % min_step_size - minimum step sizes for each dimension
    % tolerance - tolerance limit for absolute difference between consecutive function values
    % iterations_upper_bound - upper bound on the number of iterations
    % algorithm takes to compute optimal value
    % objective_function - function to be optimized
    %
    % Output:
    %
    % X - values of optimized variables
    % optimal_value - optimal function value
    % no_of_iterations - total number of iterations the algorithm underwent
    %
    
    new_X = X;
    optimal_value = feval(objective_function, new_X, no_of_variables);
    last_best_function_value = 100 * optimal_value + 100;
    
    iterate_again = true;
    no_of_iterations = 0;
    
    while iterate_again
      
      no_of_iterations = no_of_iterations + 1;

      % stop if the 'no_of_iterations' exceed the
      % 'iterations_upper_bound'
      if no_of_iterations > iterations_upper_bound
        break;
      end

      % exploratory search phase of Hooke Jeeves algorithm
      X = new_X;
      movement = zeros(no_of_variables);
      for i = 1 : no_of_variables
        better_function_value = true;
        while better_function_value
          new_X_temp = new_X(i);
          new_X(i) = new_X_temp + step_size(i);
          new_value = feval(objective_function, new_X, no_of_variables);
          if new_value < optimal_value
            optimal_value = new_value;
            movement(i) = 1;
          else
            new_X(i) = new_X_temp - step_size(i);
            new_value = feval(objective_function, new_X, no_of_variables);
            if new_value < optimal_value
              optimal_value = new_value;
              movement(i) = 1;
            else
              new_X(i) = new_X_temp;
              better_function_value = false;
            end
          end
        end
      end
    
      any_movement = sum(movement);

      % pattern search phase of Hooke Jeeves algorithm
      if any_movement > 0
        delta_X = new_X - X;
        lambda = 0.5;
        lambda = linsearch(X, no_of_variables, delta_X, lambda, objective_function);
        new_X = X + lambda * delta_X;
      end
    
      optimal_value = feval(objective_function, new_X, no_of_variables);
    
      % reduce the step size for the dimensions where no move took place
      for i = 1 : no_of_variables
        if movement(i) == 0
          step_size(i) = step_size(i) / 2;
        end
      end
    
      % stop if the difference between 2 consecutive optimal function
      % values is less than 'tolerance'
      if abs(optimal_value - last_best_function_value) < tolerance
        break
      end
    
      last_best_function_value = optimal_value;

      % stop if the 'step_size' for any variable becomes less than the
      % 'min_step_size' corresponding to that variable
      halt = true;
      for i = 1 : no_of_variables
        if step_size(i) >= min_step_size(i)
          halt = false;
        end
      end
      iterate_again = ~halt;

    end

end

% This function returns the value of 'objective_function' computed at 'X +
% lambda * delta_X'
function y = function_evaluation(no_of_variables, X, delta_X, lambda, objective_function)

    X = X + lambda * delta_X;
    y = feval(objective_function, X, no_of_variables);

end

% Using Newton's Line Search method to find the optimal
% 'objective_function' value with initial point 'X' and in the direction of 'direction'.
function lambda = linsearch(X, no_of_variables, direction, lambda, objective_function)

    iterations_upper_bound = 1000;
    tolerance = 1e-9;
    
    iter = 0;
    iterate_again = true;
    while iterate_again
        iter = iter + 1;
        if iter > iterations_upper_bound
          lambda = 0;
          break
        end
        
        h = 0.01 * (1 + abs(lambda));
        f_0 = function_evaluation(no_of_variables, X, direction, lambda, objective_function);
        f_1 = function_evaluation(no_of_variables, X, direction, lambda+h, objective_function);
        f_2 = function_evaluation(no_of_variables, X, direction, lambda-h, objective_function);
        deriv1 = (f_1 - f_2) / 2 / h;
        deriv2 = (f_1 - 2 * f_0 + f_2) / h ^ 2;
        difference = deriv1 / deriv2;
        lambda = lambda - difference;

        if abs(difference) < tolerance
          iterate_again = false;
        end
    end

end
