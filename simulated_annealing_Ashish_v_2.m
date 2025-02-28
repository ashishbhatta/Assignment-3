clear;
clc;

%minimum and maximum size of PV1 and PV2
bounds = ([[0, 22710]; [0, 22710]]);

%Initial Temperature
t_init = 0.1;

%Final Temperature
t_final = 0.01;

%Cooling rates
cooling_rates = [0.99];

% Initialize a structure to store the results
results = struct();

% Loop over the cooling rates
for i = 1:length(cooling_rates)
    cooling_rate = cooling_rates(i);
    
    % Call the simulated annealing function with the current cooling rate
    [iterations, current_solutions, final_solution, final_quality, temperatures] = ... 
    simulated_annealing(bounds, t_init, t_final, cooling_rate);


    filename = strcat('Results_for_cooling_rate_', num2str(cooling_rate), '.xlsx');
    
    % Concatenating data into one matrix
    data = [temperatures', current_solutions'];
    
    % Writing to a single sheet
    writematrix(data, filename, 'Sheet', 'Combined Data');

     
    % Display the result
    disp(['Result for cooling rate ', num2str(cooling_rate), ' is:']);
    disp(final_solution);
    
    % Store the results in the structure
    % field_name = sprintf('cooling_rate_%.2f', cooling_rate);
    % field_name = strrep(field_name, '.', '_');  % Replace '.' with '_'
    % results.(field_name) = struct('solution', final_solution, 'quality', final_quality, 'temperatures', temperatures);
end

function [iterations, current_solutions, current_params, current_solution, temperatures] = ... 
    simulated_annealing(bounds, t_init, t_final, cooling_rate)
    
    % Initialize the current solution and temperature
    num_params = size(bounds,1);
    current_params = zeros(1, num_params);  % Initialize the array for parameters

    for i = 1:num_params
        current_params(1,i) = rand .* (bounds(i, 2) - bounds(i, 1)) + bounds(i, 1);
    end
    
    % rand gives a number between 0 and 1 with normal distribution

    current_solution = power_factory(current_params);
    t_current = t_init;
    temperatures = t_current;
    current_solutions = current_solution;
    iterations = 0;

    % Print initial state
    fprintf('temp = %f\t', t_current);
    fprintf('PV_size = %.6f,%.6f\t ', current_params(1,1),current_params(1,2));
    fprintf('current Loss %.6f\n', current_solution);

    
     % Iterate until the temperature is below the final temperature
    while (t_current > t_final) 
        % Randomly perturb the current solution
        for i = 1:num_params
        perturbed_params(1,i) = rand * (bounds(i, 2) - bounds(i, 1)) + bounds(i, 1);
        end
        
        perturbed_solution = power_factory(perturbed_params);

        % Calculate the change in solution quality
        delta = perturbed_solution - current_solution;

        % If the perturbed solution is better, accept it as the new current solution
        if (delta < 0)
            current_params = perturbed_params;
            current_solution = perturbed_solution;
        % If the perturbed solution is worse, accept it with a certain probability
        else
            probability = exp(-delta / t_current);
            rand_val = rand;
            if rand_val < probability
                current_params = perturbed_params;
                current_solution = perturbed_solution;
            end
        end
        % Decrease the temperature according to the cooling rate
        t_current = t_current*cooling_rate;
        temperatures = [temperatures, t_current];
        iterations = [iterations,i];
        current_solutions = [current_solutions,current_solution];
        
        fprintf('temp = %f\t', t_current);
        fprintf('PV_size = %.6f,%.6f\t ', current_params(1,1),current_params(1,2));
        fprintf('current Loss %.6f\n', current_solution);
    end
end


function current_solution = power_factory(variables)
    %write the variables to a power file
    file1=fopen('Power.txt','w'); %w means it will overwrite all previous content of file
    fprintf(file1, '%f ', variables(1:end-1));  % Write all but the last with space
    fprintf(file1, '%f\n', variables(end));  % Write the last variable without a space
    fclose(file1);


    % a =0 means signals powerfactory to perform loadflow
    % a = 1 means power factory has perfomed loadflow and result has been
    % exported
    %a = 2 means stop the program in powerfactory

    %change the flag of  file to signal power factory to start Load Flow
    a=0; 
    file_2=fopen('Couple.txt','w');
    fprintf(file_2,'%d',a);
    fclose(file_2);


    % Wait until PowerFactory completes the Loadflow
    while a==0
    pause(0.1);
    a=load('Couple.txt');
    end

    % Read loss value from Power Factory
    file_3 = fopen('Loss.txt','r');
    current_solution=fscanf(file_3,'%f'); %Read loss value which act as o/p of function
    fclose(file_3);
end
