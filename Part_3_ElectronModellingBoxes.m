set(0, 'DefaultFigureWindowStyle', 'docked')
close all

%Duration of simulation
num_steps = 1000;

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
num_part = 50000;

%Silicon Dimensions
length_silicon = 200; %in nm
width_silicon = 100;
silicon = zeros(width_silicon, length_silicon);

%Boxes!
num_boxes = 2;
Boxes = {};
Box{1}.x =[80 120];
Box{1}.y =[0 40];

Box{2}.x =[80 120];
Box{2}.y =[60 100];

%Are box boundaries diffusive?
diffusive = 1;


%particles are assigned random x and y positions
part.position = zeros(num_part, 2);
part.position(:,1) = length_silicon * rand (num_part, 1) ;
part.position(:,2) = width_silicon * rand (num_part, 1);

% Remove particles from box region
for b = 1: num_boxes
    q = 1;
    while q <num_part
        q = q+1;
        if part.position(q,1) > Box{b}.x (1) && part.position(q,1) < Box{b}.x (2) && part.position(q,2) > Box{b}.y (1) && part.position(q,2) < Box{b}.y (2)
               part.position(q,1) = length_silicon * rand;
               part.position(q,2) = width_silicon * rand;
               q = 1;
        end
    end
end



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


temperature_map = 300*ones(width_silicon, length_silicon);

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
        my_position_y = round(part.position(n,2));
        if my_position_y <= 0
            my_position_y = 1;
        end
        if my_position_y >100
            my_position_y = 100;
        end
        my_position_x = round(part.position(n,1));
        if my_position_x <= 0
            my_position_x = 1;
        end
        if my_position_x > 200
            my_position_x = 200;
        end
    v_tot_square = (part.velocity(n,1).^2 + part.velocity(n,2).^2);
    KE_part = (1/2) * C.mn * v_tot_square;
    temperature_map (my_position_y, my_position_x) = temperature_map (my_position_y, my_position_x)+KE_part / C.k;


    %Scattering
        P_scat = 1-exp(-(delta_time_step*delta_t+delta_t)/tau);
        delta_time_step = delta_time_step + 1;
      
        if (P_scat > rand())
            part.collisions (n) = part.collisions (n)+1;
            delta_time_step = 0; %reset time step
            %velocity reassigned
            new_random_velocity = v_Th/3 * randn(1,1) + v_Th; 
            part.phi(n) = 2*pi* rand(); 
            part.velocity(n,1) = cos(part.phi(n)) * new_random_velocity;
            part.velocity(n,2) = sin(part.phi(n)) * new_random_velocity;
        end
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
            part.position (n, 1) = 0; 
            else 
            part.position(n, 1) = length_silicon; 
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

    


    %electron density
    for m =1:num_part      
        my_position_y = round(part.position(m,2));
        if my_position_y <= 0
            my_position_y = 1;
        end
        if my_position_y >100
            my_position_y = 100;
        end
        my_position_x = round(part.position(m,1));
        if my_position_x <= 0
            my_position_x = 1;
        end
        silicon(my_position_y, my_position_x)  = silicon(my_position_y, my_position_x)  + 1;
        temperature_map (my_position_y, my_position_x) = temperature_map (my_position_y, my_position_x) +(mean(temperatures))/num_part*silicon(my_position_y, my_position_x) ;
    end

    figure
    electron_density = silicon/(length_silicon*width_silicon);
     surf(electron_density);
     title('Electron Density')
    ylabel('y (nm)')
    xlabel('x(nm)')
    zlabel('Electrons per Nano Metre')

 
%     for w =1:width_silicon
%         for l =1:length_silicon
%             temperature_map(w,l)= -KE/log(silicon(w,l)/num_part) /C.k;
%         end
%     end

    figure
    surface(temperature_map);
    title('Temperature Map')
    ylabel('y (nm)')
    xlabel('x(nm)')
    zlabel('Temperature')


    
% 
%     mean_collision_time = duration/mean(part.collisions);
%     mean_free_path = mean_collision_time * v_Th
%     mean_temp= mean(temperatures);



    %Every particle gets its own colour
    colours = jet(num_part);
    
    %This plots the linear trajectories of all the particles
    figure
    for m =1:num_part
    scatter(all_x_positions(m,:),all_y_positions(m,:),'.', 'color', colours(m))
    hold on
    end

    
    for b =1:num_boxes
    rectangle('Position', [Box{b}.x(1) Box{b}.y(1) Box{b}.x(2) - Box{b}.x(1) Box{b}.y(2) - Box{b}.y(1)]);
    end
 
     
    allTrajTitle = sprintf('Trajectories of Electrons, with Diffusive Edged Boxes');
    title(allTrajTitle)
    axis([0 length_silicon 0 width_silicon])
    ylabel('y, (nm)')
    xlabel('x, (nm)')

 


