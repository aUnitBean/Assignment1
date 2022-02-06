

set(0, 'DefaultFigureWindowStyle', 'docked')
close all

%Duration of simulation
num_steps = 10000;

%Silicon Temperature
T = 300;

%Constants
C.m0 = 9.11 *10 ^ (-31);
C.mn = 0.26 * C.m0;
C.k = 1.381 * 10 ^ (-23);

%Thermal velocity, mean time between collisions,
%Mean free path
v_Th = sqrt(2*C.k*T/C.mn)/10^-9;
tau = 0.2 * 10 ^(-12);
lambda = 3.74 * 10 ^(-8);

%Time Step
delta_t = tau/100;

%Number of Particles
num_part = 50;
part_present = 0;
%Silicon Dimensions
length_silicon = 200; %in nm
width_silicon = 100;
silicon = zeros(width_silicon, length_silicon);

%Boxes!
num_boxes = 4;
Boxes = {};
Box{1}.x =[80 120];
Box{1}.y =[0 40];

Box{2}.x =[80 120];
Box{2}.y =[60 100];

Box{3}.x =[140 160];
Box{3}.y =[45 55];

Box{4}.x =[40 60];
Box{4}.y =[45 55];

%Are box boundaries diffusive?
diffusive = 0;


%particles are assigned random x and y positions
part.position = zeros(num_part, 2);

%particles are assigned x and y velocities
part.velocity = zeros(num_part, 2);

% temperature_map = 300*ones(width_silicon, length_silicon);

all_x_positions = zeros(num_part, num_steps);
all_y_positions = zeros(num_part, num_steps);
all_x_positions(:,1) = part.position(:,1);
all_y_positions(:,1) = part.position(:,2);


%Array of all temperatures
temperatures = zeros(num_steps,1);

%For collisions
num_collisions = 0;
part.collisions = zeros(num_part,1);
duration = num_steps * delta_t;
%--------------------------------------------------------------------------
for i = 1:num_steps


%This is for the temperature
    v_mean = mean(part.velocity(:,1).^2 + part.velocity(:,2).^2);
    KE = (1/2) * C.mn * v_mean;
    T = KE / C.k;
    temperatures(i) = T;   
%     
    %This is the live plot
%     hold on
    
%     scatter(part.position(:,1),part.position(:,2),'.', 'b');
%      myTitle = sprintf('Electron Trajectories in Silicon, with Impurities');
%     title(myTitle)
%     ylabel('y, (nm)')
%     xlabel('x, (nm)')
%     axis([0 length_silicon 0 width_silicon])
%     pause(0.01)

    P_inject = rand();
    if(P_inject>rand())        
        if part_present < num_part
        part_present = part_present+1;
        part.phi(part_present) = (pi/2-(-pi/2))* rand()-pi/2; 
        v_rand = v_Th/3 * randn(1,1) + v_Th; 
    
        part.velocity(part_present,1) = cos(part.phi(part_present)) * v_rand;
        part.velocity(part_present,2) = sin(part.phi(part_present)) * v_rand;
    
        part.position(part_present,1) = 0;
        part.position(part_present,2) = width_silicon * rand ();
        end
    end

    %Position Updates
    part.position = part.position + part.velocity * delta_t; %time converted to nano seconds
    delta_time_step = 0;

    %Checking Boundary Conditions and scattering
    
    for n = 1:part_present

    %Scattering
%         P_scat = 1-exp(-(delta_time_step*delta_t+delta_t)/tau);
%         delta_time_step = delta_time_step + 1;
%       
%         if (P_scat > rand())
%             part.collisions (n) = part.collisions (n)+1;
%             delta_time_step = 0; %reset time step
%             %velocity reassigned
%             new_random_velocity = v_Th/3 * randn(1,1) + v_Th; 
%             part.phi(n) = 2*pi* rand(); 
%             part.velocity(n,1) = cos(part.phi(n)) * new_random_velocity;
%             part.velocity(n,2) = sin(part.phi(n)) * new_random_velocity;
%         end
        %boxes
        for b=1:num_boxes
           while part.position(n,1) > Box{b}.x (1) && part.position(n,1) < Box{b}.x (2) && part.position(n,2) > Box{b}.y (1) && part.position(n,2) < Box{b}.y (2)
              part.position (n,:) = part.position(n,:) - part.velocity(n,:) * delta_t;
               if diffusive ~=0
                new_random_velocity = v_Th/3 * randn(1,1) + v_Th; 
                part.phi(n) = 2*pi* rand(); 
                part.velocity(n,1) = cos(part.phi(n)) * new_random_velocity;
                part.velocity(n,2) = sin(part.phi(n)) * new_random_velocity;
              else
                   if part.position (n,1) < Box{b}.x (1) || part.position(n,1) > Box{b}.x (2)
                      part.velocity(n,1) =-part.velocity(n,1) ;
                   end
                   if part.position (n,2) < Box{b}.y (1) || part.position(n,2) > Box{b}.y (2)
                      part.velocity(n,2) =-part.velocity(n,2) ;
                   end  
              end
           end

        end


        if  part.position(n, 1) > length_silicon || part.position(n, 1) < 0
            if  part.position(n, 1) > length_silicon 
            part.position (n,:) = part.position(n,:) - part.velocity(n,:) * delta_t;
            part.velocity(n,1) = -part.velocity(n,1);
            end
        end
        if  part.position(n, 2) >= width_silicon || part.position(n, 2) <= 0 
            part.position (n,:) = part.position(n,:) - part.velocity(n,:) * delta_t;
            part.velocity(n,2) = -part.velocity(n,2);
        end
    end
       
    % I have arrays for storing all of the positions so I can plot them all
    % at the end
    all_x_positions(:,i) = part.position(:,1);
    all_y_positions(:,i) = part.position(:,2);    

end

  

    %Every particle gets its own colour
    colours = jet(num_part);
    
    %This plots the linear trajectories of all the particles
    
    for m =1:part_present
    scatter(all_x_positions(m,:),all_y_positions(m,:),'.', 'color', colours(m))
    hold on
    end

    
    for b =1:num_boxes
    rectangle('Position', [Box{b}.x(1) Box{b}.y(1) Box{b}.x(2) - Box{b}.x(1) Box{b}.y(2) - Box{b}.y(1)]);
    end
 
     
    allTrajTitle = sprintf('Trajectories of Injected Electrons, No Scattering and Lots of Boxes');
    title(allTrajTitle)
    axis([0 length_silicon 0 width_silicon])
    ylabel('y, (nm)')
    xlabel('x, (nm)')

 


