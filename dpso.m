
clc;
clear;
close all;



nPop = 80;
nVar = 10;
VarSize = [1 nVar];         % Matrix Size of Decision Variables
VarMin = 1;     % Lower Bound of Decision Variables
VarMax = 5;     % Upper Bound of Decision Variables
MaxVelocity = 0.2*(VarMax-VarMin);
MinVelocity = -MaxVelocity;
MaxIt = 100;
wmax = 0.9;
wmin = 0.4;
[w1,w2]=deal(wmax);
c1=1.4962;
c2=1.4962;
c3=1.4962;
c4=1.4962;

CostFunction = @(x) Sphere(x);



    empty_particle.Position = [];
    empty_particle.Velocity = [];
    empty_particle.Cost = [];
    empty_particle.Best.Position = [];
    empty_particle.Best.Cost = [];

    % Create Population Array
    particle = repmat(empty_particle, nPop, 1);
	
    % Initialize Global Best
	
	% Globalbest of all particles(three sub_swarms)
    GlobalBest.Cost = inf;			
	% Globalbest of each subswarm
	emptyGlobalBest.Cost=inf;
    emptyGlobalBest.Position=[];
    GlobalBest1 = repmat(emptyGlobalBest,3,1);               % Initialize GlobalBest of each sub_swarm
    BestCosts=zeros(MaxIt,1);
    ShowIterInfo = true;
 
 %% Initial all particles
 for i=1:nPop

        % Generate Random Solution
        particle(i).Position = unifrnd(VarMin, VarMax, VarSize);

        % Initialize Velocity
        particle(i).Velocity = zeros(VarSize);

        % Evaluation
        particle(i).Cost = CostFunction(particle(i).Position);

        % Update the Personal Best
        particle(i).Best.Position = particle(i).Position;
        particle(i).Best.Cost = particle(i).Cost;

        % Update Global Best
        if particle(i).Best.Cost < GlobalBest.Cost
            GlobalBest = particle(i).Best;
        end

 end
%  cat(2,particle.Cost)
%% rank the particles by fitness value
	for i=1:nPop-1
        for j=i+1:nPop
            if particle(i).Cost>particle(j).Cost
                t=particle(j);
                particle(j)=particle(i);
                particle(i)=t;
            end
        end
	end
	
%% main loop
t=0;	%parameter to calculate avg2
for it=1:MaxIt
	% calculate favg favg1 favg2 and GlobalBest1
    for i = 1:nPop
		fmin = particle(1).Cost;
		fmax = particle(nPop).Cost;
		
		
        fsum =0;
        for j = 1:nPop
            fsum = fsum + particle(j).Cost;
        end
        favg = fsum/nPop;
		
        fsum1=0;
        fsum2=0;
		for j=1:nPop
			if particle(i).Cost < favg
			fsum1 = fsum1+particle(j).Cost;
			t=t+1;
            else
			fsum2 = fsum2+particle(j).Cost;
			end
		end
		
		favg1=fsum1/t;
		favg2=fsum2/(nPop-t);
	
		for j = 1:nPop
			if particle(i).Cost <=favg1
				if particle(j).Cost < GlobalBest1(1).Cost
					GlobalBest1(1) = particle(i).Best;
				end

			elseif favg1<particle(i).Cost <=favg2
				if particle(j).Cost < GlobalBest1(2).Cost
					GlobalBest1(2) = particle(i).Best;
				end
			elseif particle(i).Cost > favg2
				if particle(j).Cost < GlobalBest1(3).Cost
					GlobalBest1(3) = particle(i).Best;
				end
			end
		end
	
	%evaluate particles
	
		% superior sub-swarm
        if particle(i).Cost <=favg1
			
			w1 = w1-(w1-wmin)*abs((particle(i).Cost-favg1)/(fmin-favg1)); 		% several mistakes here
			% Update Velocity
            particle(i).Velocity = w1*particle(i).Velocity ...
                + c1*rand(VarSize).*(particle(i).Best.Position - particle(i).Position) ...
                + c2*rand(VarSize).*(GlobalBest1(1).Position - particle(i).Position);

            % Apply Velocity Limits
            particle(i).Velocity = max(particle(i).Velocity, MinVelocity);
            particle(i).Velocity = min(particle(i).Velocity, MaxVelocity);
            
            % Update Position
            particle(i).Position = particle(i).Position + particle(i).Velocity;
            
            % Apply Lower and Upper Bound Limits
            particle(i).Position = max(particle(i).Position, VarMin);
            particle(i).Position = min(particle(i).Position, VarMax);

            % Evaluation
            particle(i).Cost = CostFunction(particle(i).Position);

            % Update Personal Best
            if particle(i).Cost < particle(i).Best.Cost

                particle(i).Best.Position = particle(i).Position;
                particle(i).Best.Cost = particle(i).Cost;

                % Update Global Best
                if particle(i).Best.Cost < GlobalBest1(1).Cost
                    GlobalBest1(1) = particle(i).Best;
                end            

            end
			
			
		% common sub-swarm
		elseif favg1<particle(i).Cost <=favg2
			w2 = wmax-(w1-wmin)*(it/MaxIt)^2;
			% Update Velocity
            particle(i).Velocity = w2*particle(i).Velocity ...     
                + c2*rand(VarSize).*(GlobalBest1(1).Position - particle(i).Position)...
				+ c3*rand(VarSize).*(GlobalBest.Position - particle(i).Position) ...
				+ c4*rand(VarSize).*(particle(i).Best.Position + particle(i).Position);
				
            % Apply Velocity Limits
            particle(i).Velocity = max(particle(i).Velocity, MinVelocity);
            particle(i).Velocity = min(particle(i).Velocity, MaxVelocity);
            
            % Update Position
            particle(i).Position = particle(i).Position + particle(i).Velocity;
            
            % Apply Lower and Upper Bound Limits
            particle(i).Position = max(particle(i).Position, VarMin);
            particle(i).Position = min(particle(i).Position, VarMax);

            % Evaluation
            particle(i).Cost = CostFunction(particle(i).Position);

            % Update Personal Best
            if particle(i).Cost < particle(i).Best.Cost

                particle(i).Best.Position = particle(i).Position;
                particle(i).Best.Cost = particle(i).Cost;

                % Update Global Best
                if particle(i).Best.Cost < GlobalBest1(2).Cost
                    GlobalBest1(2) = particle(i).Best;
                end            

            end
			
			
		% inferior sub-swarm
		elseif particle(i).Cost > favg2
			% Update Velocity
            %particle(i).Velocity = particle(i).Velocity;

            % Apply Velocity Limits
            particle(i).Velocity = max(particle(i).Velocity, MinVelocity);
            particle(i).Velocity = min(particle(i).Velocity, MaxVelocity);
            
            % Update Position
            particle(i).Position = VarMin +rand(VarSize)*(VarMax-VarMin);
            
            % Apply Lower and Upper Bound Limits
            particle(i).Position = max(particle(i).Position, VarMin);
            particle(i).Position = min(particle(i).Position, VarMax);

            % Evaluation
            particle(i).Cost = CostFunction(particle(i).Position);

            % Update Personal Best
            if particle(i).Cost < particle(i).Best.Cost

                particle(i).Best.Position = particle(i).Position;
                particle(i).Best.Cost = particle(i).Cost;

                % Update Global Best
                if particle(i).Best.Cost < GlobalBest1(3).Cost
                    GlobalBest1(3) = particle(i).Best;
                end            

            end
        end
			
        for j =1:3
            if(GlobalBest1(j).Cost < GlobalBest.Cost)
                GlobalBest = GlobalBest1(j);
            end
        end
    
    end
    % Store the Best Cost Value
    BestCosts(it) = GlobalBest.Cost;
    % Display Iteration Information
    if ShowIterInfo
        disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCosts(it))]);
    end
    
end

%% Results

    figure;
    %plot((BestCosts,'LineWidth',2); 
    semilogy(BestCosts,'LineWidth',2);      
    xlabel('Iteration');
    ylabel('Best Cost');