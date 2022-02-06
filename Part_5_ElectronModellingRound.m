set(0, 'DefaultFigureWindowStyle', 'docked')
close all

%Duration of simulation
num_steps = 8000;

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
num_part = 100;

%Silicon Dimensions
length_silicon = 200; %in nm
width_silicon = 100;
silicon = zeros(width_silicon, length_silicon);



%circles!
circle_points = 1000;
circle_angle = linspace(0,2*pi,circle_points);
Circle = {};
Circle{1}.centre = [100 50];
Circle{1}.radius = 20;
Circle{1}.x = Circle{1}.radius* cos(circle_angle)'+Circle{1}.centre(1);
Circle{1}.y = Circle{1}.radius* sin(circle_angle)'+Circle{1}.centre(2);

radius_x = 0;
radius_y = 0;


% %Boxes!
% num_boxes = 2;
% Boxes = {};
% Box{1}.x =[80 120];
% Box{1}.y =[0 40];
% 
% Box{2}.x =[80 120];
% Box{2}.y =[60 100];

%Are box boundaries diffusive?
diffusive = 1;


%particles are assigned random x and y positions
part.position = zeros(num_part, 2);
part.position(:,1) = length_silicon * rand (num_part, 1) ;
part.position(:,2) = width_silicon * rand (num_part, 1);


% Remove particles from circle region
q = 1;
in = 0;
while q <num_part
    q = q+1;
    in =  inpolygon(part.position(q,1), part.position(q,2), Circle{1}.x, Circle{1}.y);
    if in 
           part.position(q,1) = length_silicon * rand;
           part.position(q,2) = width_silicon * rand;
           q = 1;
    end
end


%every particle has a random angle 
part.phi = 2*pi* rand(num_part, 1); 

%particles are assigned x and y velocities
part.velocity = zeros(num_part, 2);
v_rand = v_Th/3 * randn(num_part,1) + v_Th;  %random velocities from normal distributions, varriance is a quarter of the average
part.velocity(:,1) = cos(part.phi) .* v_rand;
part.velocity(:,2) = sin(part.phi) .* v_rand;


all_x_positions = zeros(num_part, num_steps);
all_y_positions = zeros(num_part, num_steps);
all_x_positions(:,1) = part.position(:,1);
all_y_positions(:,1) = part.position(:,2);


delta_time_step=0;

for i = 1:num_steps

%This is for the temperature
%     
  %  This is the live plot
%     scatter(Circle{1}.x , Circle{1}.y, 'b');
%     hold on
%     scatter(part.position(:,1),part.position(:,2),'.', 'b');
%     myTitle = sprintf('Electron Trajectories in Silicon, with Impurities, Temperature: %d k', T);
%     title(myTitle)
%     ylabel('y, (nm)')
%     xlabel('x, (nm)')
%     axis([0 length_silicon 0 width_silicon])
%     pause(0.01)
%     hold off
    %Position Updates
    part.position = part.position + part.velocity * delta_t; %time converted to nano seconds

    %Checking Boundary Conditions and scattering
    
    for n = 1:num_part

    %Scattering
        P_scat = 1-exp(-(delta_time_step*delta_t+delta_t)/tau);
        delta_time_step = delta_time_step + 1;
      
        if (P_scat > rand())
            delta_time_step = 0; %reset time step
            %velocity reassigned
            new_random_velocity = v_Th/3 * randn(1,1) + v_Th; 
            part.phi(n) = 2*pi* rand(); 
            part.velocity(n,1) = cos(part.phi(n)) * new_random_velocity;
            part.velocity(n,2) = sin(part.phi(n)) * new_random_velocity;
        end
          %circles!
          
          
          
          [in_circ, on_circ] =  inpolygon(part.position(n,1), part.position(n,2), Circle{1}.x, Circle{1}.y);
          
          if in_circ|| on_circ
          radius_x =  Circle{1}.centre(1)-part.position(n,1);
          radius_y =  Circle{1}.centre(2)-part.position(n,2) ;
          phi_circ = atan(radius_y/radius_x);

          while in_circ 
               part.position (n,:) = part.position(n,:) - part.velocity(n,:) * delta_t;
               in_circ =  inpolygon(part.position(n,1), part.position(n,2), Circle{1}.x, Circle{1}.y);
          end
          
                part.phi(n) = phi_circ -2*(part.phi(n)-phi_circ);
                part.velocity(n,1) = part.velocity(n,1)*cos(part.phi(n));
                part.velocity(n,2) = part.velocity(n,1)*sin(part.phi(n));
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
            part.phi(n) = part.phi(n)+pi();
            part.velocity(n,2) = part.velocity(n,1)*sin(part.phi(n));
        end
    end
       
    % I have arrays for storing all of the positions so I can plot them all
    % at the end
    all_x_positions(:,i) = part.position(:,1);
    all_y_positions(:,i) = part.position(:,2);    

end

  

    %Every particle gets its own colour
    colours = jet(num_part);
%     
%     %This plots the linear trajectories of all the particles
    figure
    for m =1:num_part
    scatter(all_x_positions(m,:),all_y_positions(m,:),'.', 'color', colours(m))
    hold on
    end

    plot(Circle{1}.x,Circle{1}.y);
 
     
    allTrajTitle = sprintf('Trajectories of Electrons, Circular Box');
    title(allTrajTitle)
    axis([0 length_silicon 0 width_silicon])
    ylabel('y, (nm)')
    xlabel('x, (nm)')

 


