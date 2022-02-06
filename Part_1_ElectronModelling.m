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
v_Th = sqrt(2*C.k*T/C.mn);
tau = 0.2 * 10 ^(-12);
lambda = 3.74 * 10 ^(-8);

%Time Step
delta_t = 1 * 10^-5;

%Number of Particles
num_part = 10;

%Silicon Dimensions
length_silicon = 200;
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
part.velocity(:,1) = cos(part.phi) * v_Th;
part.velocity(:,2) = sin(part.phi) * v_Th;


all_x_positions = zeros(num_part, num_steps);
all_y_positions = zeros(num_part, num_steps);
all_x_positions(:,1) = part.position(:,1);
all_y_positions(:,1) = part.position(:,2);

for i = 1:num_steps

    v_mean_squared = mean(part.velocity(:,1).^2 + part.velocity(:,2).^2);
    KE = (1/2) * C.mn * v_mean_squared;
    T = KE   / C.k;
    scatter(i,T);
    title("Temperature of Silicon")
    ylabel('Temperature (Kelvin)')
    xlabel('Steps')
    hold on
%   This is the live plot
%     scatter(part.position(:,1),part.position(:,2),'.', 'b');
%     myTitle = sprintf('Electron Trajectories in Silicon, no Impurities, Temperature: %d k', T);
%     title(myTitle)
%     ylabel('y, (nm)')
%     xlabel('x, (nm)')
%     axis([0 length_silicon 0 width_silicon])
%     pause(0.01)
 
    %Position Updates
    part.position = part.position + part.velocity * delta_t;
    
    %Checking Boundary Conditions
    for n = 1:num_part
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
    %Every particle gets its own colour
    colours = jet(num_part);
    
  %  This plots the linear trajectories of all the particles
    figure
    for m =1:num_part
    scatter(all_x_positions(m,:),all_y_positions(m,:),'X', 'color', colours(m))
    hold on
    end
    title(myTitle)
    axis([0 length_silicon 0 width_silicon])
    ylabel('y, (nm)')
    xlabel('x, (nm)')


