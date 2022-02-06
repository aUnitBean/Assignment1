set(0, 'DefaultFigureWindowStyle', 'docked')
close all

%Duration of simulation
num_steps = 2000;

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
delta_t = tau/2;

%Number of Particles
num_part = 1000;

%Silicon Dimensions
length_silicon = 200; %in nm
width_silicon = 100;
silicon = zeros(width_silicon, length_silicon);

%particles are assigned random x and y positions
part.position = zeros(num_part, 2);
part.position(:,1) = length_silicon * rand (num_part, 1);
part.position(:,2) = width_silicon* rand (num_part, 1);

%every particle has a random angle 
part.phi = 2*pi* rand(num_part, 1); 

%particles are assigned x and y velocities
part.velocity = zeros(num_part, 2);
v_rand = v_Th/3 * randn(num_part,1) + v_Th;  %random velocities from normal distributions, varriance is a quarter of the average
part.velocity(:,1) = cos(part.phi) .* v_rand;
part.velocity(:,2) = sin(part.phi) .* v_rand;

% velocities = sqrt(part.velocity(:,1).^2+part.velocity(:,2).^2);
% hist (velocities);
% VelocityTitle = sprintf('Histogram for Velocities of Electrons in Silicon,  Mean: %d m/s', mean(velocities));
% 
% title(VelocityTitle)
% ylabel('Occurence')
% xlabel('Velocity (m/s)')


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
%      myTitle = sprintf('Electron Trajectories in Silicon, with Impurities, Temperature: %d k', T);
%     title(myTitle)
%     ylabel('y, (nm)')
%     xlabel('x, (nm)')
%     axis([0 length_silicon 0 width_silicon])
%     pause(0.01)

    %Position Updates
    part.position = part.position + part.velocity * delta_t; %time converted to nano seconds
    delta_time_step = 0;

    %Checking Boundary Conditions and scattering
    
    for n = 1:num_part

        P_scat = 1-exp(-(delta_time_step*delta_t+delta_t)/tau);
        delta_time_step = delta_time_step + 1;
      
        if (P_scat > rand())
            part.collisions (n) = part.collisions (n)+1;
            delta_time_step = 0; %reset time step
            %velocity reassigned
            new_random_velocity = v_Th/3 * randn(1,1) + v_Th; 
            new_random_phi = 2*pi* rand(); 
            part.velocity(n,1) = cos(new_random_phi) * new_random_velocity;
            part.velocity(n,2) = sin(new_random_phi) * new_random_velocity;
        end

        if  part.position(n, 1) > length_silicon || part.position(n, 1) < 0
            if  part.position(n, 1) > length_silicon 
            part.position (n, 1) = 0; 
            else 
            part.position(n, 1) = length_silicon; 
            end
        end
        if  part.position(n, 2) > width_silicon || part.position(n, 2) <0 
            part.position (n,:) = part.position(n,:) - part.velocity(n,:) * delta_t;
            part.velocity(n,2) = -part.velocity(n,2);
        end
    end
       
    % I have arrays for storing all of the positions so I can plot them all
    % at the end
    all_x_positions(:,i) = part.position(:,1);
    all_y_positions(:,i) = part.position(:,2);    

end


    mean_collision_time = duration/mean(part.collisions)
    mean_free_path = mean_collision_time * v_Th
    
    
    mean_temp= mean(temperatures);



%     scatter(1:num_steps,temperatures,'.');
%     Temperature_title = sprintf('Temperature of Silicon, Mean Temperature: %d', mean_temp);
%     title(Temperature_title)
%     ylabel('Temperature (Kelvin)')
%     xlabel('Steps')
%     hold on
%     scatter(1:num_steps, mean_temp, 'r', '_')


    %Every particle gets its own colour
    colours = jet(num_part);
    
    %This plots the linear trajectories of all the particles
%     figure
%     for m =1:num_part
%     scatter(all_x_positions(m,:),all_y_positions(m,:),'.', 'color', colours(m))
%     hold on
%     end
%     allTrajTitle = sprintf("Trajectories of Electrons, with Impurities, Mean Collision Time: %d" ,mean_collision_time)
%     title(allTrajTitle)
%     axis([0 length_silicon 0 width_silicon])
%     ylabel('y, (nm)')
%     xlabel('x, (nm)')


