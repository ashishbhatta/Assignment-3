
% Initialization
clear;

n=25;           % Population size, number of Bats
t_max=50;       % Maximum number of iterations
A=1;            % Initial loudness (constant or decreasing)
r0=1;           % The initial pulse rate (constant or decreasing)
alpha=0.97;     % Parameter alpha
gamma=0.1;      % Parameter gamma

% Frequency range
Freq_min=0;     % Frequency minimum
Freq_max=2;     % Frequency maximum
t=0;            % Initialize iteration counter

% Dimensions of the search variables
d=2;            %Number of Variables in Objective Function

% Initialization of all the arrays
Freq=zeros(n,1);   % Frequency-tuning array
v=zeros(n,d);      % Equivalnet velocities or increments
Lb=zeros(1,d);   % Lower bounds
Ub=222710*ones(1,d);    % Upper bounds 222710 kW maximum size of PV

%Obtaining location of each bat
for i=1:n
    Sol(i,:) =Lb+(Ub-Lb).*rand(1,d);

    %Identifying bats location
    Fitness(i) = power_factory(Sol(i,:))
    %power_factory(sol(i,:))
end

% Find the initial best solution for Initial Population
[fmin,I]=min(Fitness);
best=Sol(I,:);
fprintf("Initial best solution is %f \n",best)

fmin_history = zeros(1, t_max);  

% Start the iterations -- the Bat Algorithm (BA) -- main loop
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

       % Evaluate new solutions/Finding new 
       Fnew = power_factory(S(i,:))
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
  fprintf("Iteration:%d best:%f fmin: %f",t,best,fmin)

  fmin_history(t+1) = fmin;
  t=t+1;  % Update iteration counter

  % % Display the results every 100 iterations
  % if ~mod(t,100),
  %    disp(strcat('Iteration = ',num2str(t)));    
  %    best, fmin 
  % end
end  % End of the main loop

% Output the best solution
disp(['Best =',num2str(best),' fmin=',num2str(fmin)]);

%Shutdown Powerfactory
a=2;
t=fopen('Couple.txt','w');
fprintf(t,'%d',a);
fclose(t);

%plot
figure;
plot(1:t_max, fmin_history, 'b-', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Fitness Value (Loss)');
title('Convergence of Bat Algorithm');
grid on;








function Fitness = power_factory(variables)
    %write the variables to a power file
    fprintf('%f,',variables)
    file1=fopen('Power.txt','w'); %w means it will overwrite all previous content of file
    fprintf("File_1 Opened \n")
    fprintf(file1, '%f ', variables(1:end-1));  % Write all but the last with space
    fprintf(file1, '%f\n', variables(end));  % Write the last variable without a space
    fclose(file1);
    fprintf("File_1 Closed \n")

    % a =0 means signals powerfactory to perform loadflow
    % a = 1 means power factory has perfomed loadflow and result has been
    % exported
    %a = 2 means stop the program in powerfactory

    %change the flag of  file to signal power factory to start Load Flow
    a=0; 
    file_2=fopen('Couple.txt','w');
    fprintf("File_2 Opened \n")
    fprintf(file_2,'%d',a);
    fclose(file_2);
    fprintf("File_2 Closed \n")

    % Wait until PowerFactory perform Loadflow
    fprintf("waiting the value of a to change to non zero\n")
    while a==0
    pause(0.1);
    a=load('Couple.txt');
    end


    fprintf("Vaue of a changes to %f\n",a)
    % Read Objective Function Values 
    file_3 = fopen('Loss.txt','r');
    Fitness=fscanf(file_3,'%f'); %Read loss value which act as o/p of function
    fclose(file_3);
end

function s = simplebounds(s, Lb, Ub)
    % Ensure that solution s remains within the specified lower and upper bounds
    s(s < Lb) = Lb(s < Lb);  % If values are below Lb, set them to Lb
    s(s > Ub) = Ub(s > Ub);  % If values are above Ub, set them to Ub
end




% Initialize the population/solutions for all bats
%for i=1:n,
% Sol(i,:)=Lb+(Ub-Lb).*rand(1,d);
%  Fitness(i)=Fun(Sol(i,:));
%end

% The cost function or objective function is given by the PowerFactory
% function z=Fun(x)
% z=sum((x-2).^2);      % Optimal solution fmin=0 at (2,2,...,2)



