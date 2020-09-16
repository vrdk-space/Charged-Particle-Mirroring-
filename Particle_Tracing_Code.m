%% Tracing Particle Trajectory in Circular Coil System

% In this work we trace the trajectory of a particle with mass = 'm' Kg and
% charge = 'q' C in a 2 or 3 circular coil system. The centre of the circular
% coils are aligned (i.e) common axis. The current in the coils, distance 
% between the coils and radii of the coils is set by the user. 

tic
clc
clear all
eps = 8.854e-12; % Permittivity of free sapce in [m^-3 kg^-1 s^4 A^2]
%% Input From User

density = 2e3;  % Density of Particle in [kg m^-3]
radius  = 1e-6; % Radius of Particle in [m]
V = 10;         % Potential of IDP in [ V = kg m^2 s^-3 A^-1]
q = 1e-12;
m = density *(4/3)*pi*(radius^3);  % Mass in [Kg]
%q = 4*pi*eps*radius*V;            % Charge of IDP in [C]

%q = 1.60217662e-19; % Charge on Electron in [C]
%m = 9.10938356e-31; % Mass of Electron in [kg]
%m = 1.6726219e-27;  % Mass of Proton in [Kg]

h = 1e-6; % Time Step
total_simulation_time = 5e-1; % Total Simulation Time Scale in [s]
total_time_steps = total_simulation_time/h;

% Coil Parameters;

I0 = 300;   % Current in Coils in [A]
L  = 5e3;  % Length of Coil in [m]
R1 = 0.600; % Radius of Coil 1 in [m]
R2 = 0.450; % Radius of Coil 2 in [m]
R3 = 0.300; % Radius of Coil 3 in [m]
R4 = 0.150; % Radius of Coil 4 in [m]

d12 = 1/3; % Distance between Coil 1 and Coil 2 in [m]
d23 = 1/3; % Distance between Coil 2 and Coil 3 in [m]
d34 = 1/3; % Distance between Coil 3 and Coil 4 in [m]

L1 = 0;
n1 = L/(2*pi*R1);
n2 = L/(2*pi*R2);
n3 = L/(2*pi*R3);
n4 = L/(2*pi*R4);

I1 = n1*I0; % Adjusted Current in coil 1 in [A]
I2 = n2*I0; % Adjusted Current in coil 2 in [A]
I3 = n3*I0; % Adjusted Current in coil 3 in [A]
I4 = n4*I0; % Adjusted Current in coil 3 in [A]

% Centre of Coil 1 is (0,0,-d12), Coil 2 is (0,0,0), 
% Coil 3 is (0,0,+d23), Coil 4 is (0,0,d23+d34)

% Current flows in anticlockwise sense in all the rings 

Coil_Input = [R1 R2 R3 R4; d12 0 d23 d34; I1 I2 I3 I4];

%% Particle Tracing

% Step 1: Input initial position and velocity (User Input)

x(1) = -0.2; % in [m]
y(1) = -0.15; % in [m]
z(1) = -0.25; % in [m]

vx(1) = 57.735; % in [m s^-1]
vy(1) = 57.735;  % in [m s^-1]
vz(1) = 57.735; % in [m s^-1]
v(1) = sqrt((vx(1)^2)+(vy(1)^2)+(vz(1)^2)); % in [m s^-1] 

% Step 2: Calculate Magnetic Field and Accleration @ Initial Position

[Bx_total(1),By_total(1),Bz_total(1)] = Magnetic_Field(Coil_Input,x(1),y(1),z(1)); % in [T]

ax(1) =  (q/m)*((vy(1)*Bz_total(1))-(vz(1)*By_total(1)));
ay(1) = -(q/m)*((vx(1)*Bz_total(1))-(vz(1)*Bx_total(1)));
az(1) =  (q/m)*((vx(1)*By_total(1))-(vy(1)*Bx_total(1)));

% Step 3: Check if its Reflection or Loss Condition

[B0_x,B0_y,B0_z] = Magnetic_Field(Coil_Input,0,0,0);    % Field components at centre of Large Ring in [T]
[Bm_x,Bm_y,Bm_z] = Magnetic_Field(Coil_Input,0,0,d23+d34); % Field components at centre of Large Ring in [T]
B0 = sqrt((B0_x^2)+(B0_y^2)+(B0_z^2));                     % Field at centre of Large Ring
Bm = sqrt((Bm_x^2)+(Bm_y^2)+(Bm_z^2));                     % Field at centre of Small Ring

B_total_initial = sqrt((Bx_total(1)^2)+(By_total(1)^2)+(Bz_total(1)^2)); % Total Field at initial point
v_initial = v(1);
Nr = (vx(1)*Bx_total(1)) + (vy(1)*By_total(1)) + (vz(1)*Bz_total(1));
Dr = B_total_initial*v_initial;
theta_0 = asind(sqrt(B0/Bm));
theta_1 = acosd(Nr/Dr);
% 
% if( theta_1 > theta_0)
%     reflection = 1
%     loss = 0
% else
%     reflection = 0
%     loss = 1
% end

v_perp = v_initial*sin(deg2rad(theta_1));
v_para = v_initial*cos(deg2rad(theta_1));
term1 = (v_perp^2)/((v_perp^2)+(v_para^2));
term2 = B0/Bm;

if(term1 > term2)
    reflection = 1
    loss = 0
else
    reflection = 0
    loss = 1
end


if (reflection == 1 || loss == 1)

    % Step 3: Calculate second point in trajectory and its velocity at that point

        x(2) = x(1) + (vx(1)*h);
        y(2) = y(1) + (vy(1)*h);
        z(2) = z(1) + (vz(1)*h);

        vx(2) = vx(1) + (ax(1)*h);
        vy(2) = vy(1) + (ay(1)*h);
        vz(2) = vz(1) + (az(1)*h);

        v(2) = sqrt((vx(2)^2)+(vy(2)^2)+(vz(2)^2)); % in [m s^-1]


    % Step 4: Calculate Magnetic Field and Accleration at second point in trajectory

        [Bx_total(2),By_total(2),Bz_total(2)] = Magnetic_Field(Coil_Input,x(2),y(2),z(2)); % in [T]

        ax(2) =  (q/m)*((vy(2)*Bz_total(2))-(vz(2)*By_total(2)));
        ay(2) = -(q/m)*((vx(2)*Bz_total(2))-(vz(2)*Bx_total(2)));
        az(2) =  (q/m)*((vx(2)*By_total(2))-(vy(2)*Bx_total(2)));

    % Step 5: Loop to calculate further points in teajectory via Second Order Central 
    % Difference Method ( Error O(x^2))
        for i = 2:total_time_steps

            k1 = (q/m)*(h/2)*Bz_total(i);
            k2 = (q/m)*(h/2)*By_total(i);
            c1 = 2*x(i) - x(i-1) - (k1*y(i-1)) + (k2*z(i-1));
    
            k3 = (q/m)*(h/2)*Bx_total(i);
            k4 = (q/m)*(h/2)*Bz_total(i);
            c2 = 2*y(i) - y(i-1) - (k3*z(i-1)) + (k4*x(i-1));

            k5 = (q/m)*(h/2)*By_total(i);
            k6 = (q/m)*(h/2)*Bx_total(i);
            c3 = 2*z(i) - z(i-1) - (k5*x(i-1)) + (k6*y(i-1));
    
            c_matrix = [c1;c2;c3];
            k_matrix = [1 -k1 k2; k4 1 -k3; -k5 k6 1];
            k_inv = inv(k_matrix);
            Sol = k_matrix\c_matrix; % Recommended by MATLAB instead of k_inv*c_matrix for 
                             % efficient calculation
            x(i+1) = Sol(1,1);
            y(i+1) = Sol(2,1);
            z(i+1) = Sol(3,1);
    
            [Bx_total(i+1),By_total(i+1),Bz_total(i+1)] = Magnetic_Field(Coil_Input,x(i+1),y(i+1),z(i+1)); % in [T]
    
        end

    % Step 6 : Calculate Velocity and Accleration of the particle with time.

        for i = 3:total_time_steps
    
            vx(i) = (x(i+1) - x(i-1))/(2*h);
            vy(i) = (y(i+1) - y(i-1))/(2*h);
            vz(i) = (z(i+1) - z(i-1))/(2*h);
            v(i) = sqrt((vx(i)^2)+(vy(i)^2)+(vz(i)^2)); % in [m s^-1]
    
            ax(i) = (x(i+1) - (2*x(i)) + x(i+1))/(h^2);
            ay(i) = (y(i+1) - (2*y(i)) + y(i+1))/(h^2);
            az(i) = (z(i+1) - (2*z(i)) + z(i+1))/(h^2);
    
        end

end
%% Plots

if (reflection == 1 || loss == 1)

        angle = 0:0.001:2*pi;

        x1 = R1*cos(angle);
        y1 = R1*sin(angle);
        z1 = ones(size(x1));
        z1 = z1.*-d12;

        x2 = R2*cos(angle);
        y2 = R2*sin(angle);
        z2 = zeros(size(x2));

        x3 = R3*cos(angle);
        y3 = R3*sin(angle);
        z3 = ones(size(x3));
        z3 = z3.*d23;

        x4 = R4*cos(angle);
        y4 = R4*sin(angle);
        z4 = ones(size(x3));
        z4 = z4.*(d23+d34);

        figure1= figure('WindowState','maximized');
        axes1 = axes('Parent',figure1);
        plot3(x,y,z,'LineWidth',2,'Color',[0 0 1],'MarkerFaceColor',[1 0 0]);hold on;
        plot3(x(1),y(1),z(1),'MarkerFaceColor',[1 0 0],'Marker','o','LineStyle','none',...
            'Color',[1 0 0]);
        plot3(x(total_time_steps+1),y(total_time_steps+1),z(total_time_steps+1),'MarkerFaceColor',[1 1 0],'Marker','o','LineStyle','none',...
            'Color',[1 0 0]);
        plot3(x1,y1,z1,'LineWidth',2,'LineStyle','--','Color',[0.24705882370472 0.24705882370472 0.24705882370472]);
        plot3(x2,y2,z2,'LineWidth',2,'LineStyle','--','Color',[0.24705882370472 0.24705882370472 0.24705882370472]);
        plot3(x3,y3,z3,'LineWidth',2,'LineStyle','--','Color',[0.24705882370472 0.24705882370472 0.24705882370472]);
        plot3(x4,y4,z4,'LineWidth',2,'LineStyle','--','Color',[0.24705882370472 0.24705882370472 0.24705882370472]);
        hold off;
        title({'Trajectory of the Charged Particle'},'FontSize',16);
        xlabel('x [m]','FontSize',14,'FontWeight','bold');
        ylabel('y [m]','FontSize',14,'FontWeight','bold');
        zlabel('z [m]','FontSize',14,'FontWeight','bold');
        set(axes1,'FontSize',12,'FontWeight','bold');
        saveas(axes1,'Trajectory.jpg')

        figure2 = figure('WindowState','maximized');
        axes2 = axes('Parent',figure2);
        plot(x,'DisplayName','x','LineWidth',2);hold on;
        plot(y,'DisplayName','y','LineWidth',2); %'LineStyle','--'
        plot(z,'DisplayName','z','LineWidth',2,'Color',[1 0 0]);hold off;
        title({'Position of Particle with Time'},'FontSize',16);
        xlabel('Time (us)','FontSize',14,'FontWeight','bold');
        ylabel('Coordinates of the Particle','FontSize',14,'FontWeight','bold');
        set(axes2,'FontSize',12,'FontWeight','bold');
        legend(axes2,'show');
        saveas(axes2,'Pos_vs_Time.jpg')

        figure3 = figure('WindowState','maximized');
        axes3 = axes('Parent',figure3);
        plot(vx,'DisplayName','vx','LineWidth',2);hold on;
        plot(vy,'DisplayName','vy','LineWidth',2);
        plot(vz,'DisplayName','vz','LineWidth',2);
        plot(v,'DisplayName','Total Velocity','LineWidth',2);hold off;
        title({'Distribution of Velocity of Particle with time'},'FontSize',16);
        xlabel('Time (us)','FontSize',14,'FontWeight','bold');
        ylabel('Velocity (m/s)','FontSize',14,'FontWeight','bold');
        set(axes3,'FontSize',12,'FontWeight','bold');
        legend(axes3,'show');
        saveas(axes3,'Vel_vs_Time.jpg')
% 
%         figure4 = figure('WindowState','maximized');
%         axes4 = axes('Parent',figure4);
%         plot(Bx_total,'DisplayName','Bx','LineWidth',2);hold on;
%         plot(By_total,'DisplayName','By','LineWidth',2);
%         plot(Bz_total,'DisplayName','Bz','LineWidth',2);hold off;
%         title({'Magnetic Field Experienced by IDP with time'},'FontSize',16);
%         xlabel('Time (us)','FontSize',14,'FontWeight','bold');
%         ylabel('Magnetic Field along Trajectory (T)','FontSize',14,'FontWeight','bold');
%         set(axes4,'FontSize',12,'FontWeight','bold');
%         legend(axes4,'show');
%         saveas(axes4,'MF_vs_Time.jpg')
% 
%         figure5 = figure('WindowState','maximized');
%         axes5 = axes('Parent',figure5);
%         plot3(Bx_total,By_total,Bz_total,'LineWidth',2,'Color',[0 0 1]);hold on;
%         plot3(Bx_total(1),By_total(1),Bz_total(1),'MarkerFaceColor',[1 0 0],'Marker','o','LineStyle','none',...
%             'Color',[1 0 0]);hold off;
%         title({'Magnetic Field Along Trajectory'},'FontSize',16);
%         xlabel('Bx (T)','FontSize',14,'FontWeight','bold');
%         ylabel('By (T)','FontSize',14,'FontWeight','bold');
%         zlabel('Bz (T)','FontSize',14,'FontWeight','bold');
%         set(axes5,'FontSize',12,'FontWeight','bold');
%         saveas(axes5,'MFOT.jpg')
    
end

toc

