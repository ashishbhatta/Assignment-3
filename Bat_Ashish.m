clear;
clc;
 
n = 25;   % Population size
t_max = 50;
A = 100;   
r0 = 100;   
alpha = 0.95;
gamma = 0.1;

Freq_min = 0;
Freq_max = 100;

d = 2;   % Number of Variables in Objective Function

Freq = zeros(n,1);
v = zeros(n,d);
Lb = zeros(1,d);
Ub = (222710 - 200000)* ones(1,d); %in kW

Sol = zeros(n,d);
Fitness = zeros(n,1);  % Ensure Fitness is a scalar for each bat

% Initialize the population/solutions
for i = 1:n
    Sol(i,:) = Lb + (Ub - Lb) .* rand(1,d);
    Fitness(i) = power_factory(Sol(i,:));  % Ensure function returns a scalar
end

% Find the best solution of the initial population
[fmin,I]=min(Fitness);
best=Sol(I,:);

t = 0;
% Main Body for BatAlgorithm
while (t<t_max)
   % Varying loundness (A) and pulse emission rate (r)
   r=r0*(1-exp(-gamma*t));
   A=alpha*A;
   % Loop over all bats/solutions
    for i=1:n,
       Freq(i)=Freq_min+(Freq_max-Freq_min)*rand;
       v(i,:)=v(i,:)+(Sol(i,:)-best)*Freq(i);
       S(i,:)=Sol(i,:)+v(i,:);
   % Check a switching condition
   if rand<r,
       S(i,:)=best+0.1*randn(1,d)*A;
   end

   % Check if the new solution is within the simple bounds
   S(i,:)=simplebounds(S(i,:),Lb,Ub);
   % Evaluate new solutions
   Fnew=power_factory(S(i,:));
   % If the solution improves or not too loudness
    if ((Fnew<=Fitness(i)) & (rand>A)),
       Sol(i,:)=S(i,:);
       Fitness(i)=Fnew;
    end
   % Update the current best solution
    if Fnew<=fmin,
       best=S(i,:);
       fmin=Fnew;
    end
   end % end of for i
  min_f(t+1,1) = t+1;
  min_f(t+1,2) = fmin;
  disp(['Iteration=',num2str(t), '  Best =',num2str(best),' fmin=',num2str(fmin)]);
  t=t+1;  % Update iteration counter
end

% Output the best solution
disp(['Best =',num2str(best),' fmin=',num2str(fmin)]);

plot(min_f(:,1),min_f(:,2));
xlabel('Iteration');
ylabel('Values');
title(['Convergence Plot ', num2str(k)]);
grid on;
saveas(gcf, ['convergence_plot_' num2str(k) '.svg']);







% Application of simple bounds/constraints
function s=simplebounds(s,Lb,Ub)
  % Apply the lower bound
  ns_tmp=s;
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb(I);
  
  % Apply the upper bounds 
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub(J);
  % Update this new move 
  s=ns_tmp;
end

function fitness = power_factory(variables)
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
    fitness=fscanf(file_3,'%f'); %Read loss value which act as o/p of function
    fclose(file_3);
end

% function fitness = power_factory(variables)
% fitness=sum((variables-2).^2);      % Optimal solution fmin=0 at (2,2,...,2)
% end

