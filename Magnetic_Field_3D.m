%% Magnetic Field Due to Rings

% In this work, each individual component of magnetic field is calculated
% in 3-Dimensional Space.
tic
% clc
% clear all
%% Physical Constants

mu = 4*pi*1e-7;  %[m kg s-2 A-2]

%% Input 

I0 = 30000;  % Current in Coils in [A]
L  = 5e3;  % Length of Coil in [m]
R1 = 0.50; % Radii of circular ring 1 in [m]
R2 = 0.50; % Radii of Circular Ring 2 in [m]
R3 = 0.50; % Radii of circular ring 3 in [m]
R4 = 0.25; % Radii of circular ring 4 in [m]

d12 = 1/3;  % Distance between the coils 1 and 2 in [m]
d23 = 1/3;  % Distance between the coils 2 and 3 in [m]
d34 = 1/3;  % Distance between the coils 3 and 4 in [m]
d = d12+d23+d34;

L1 = 0;
n1 = L1/(2*pi*R1);
n2 = L1/(2*pi*R2);
n3 = L1/(2*pi*R3);
n4 = L/(2*pi*R4);

I1 = n1*I0; % Adjusted Current in coil 1 in [A]
I2 = n2*I0; % Adjusted Current in coil 2 in [A]
I3 = n3*I0; % Adjusted Current in coil 3 in [A]
I4 = n4*I0; % Adjusted Current in coil 3 in [A]

% Centre of Coil 1 : (0,0,-d12), Centre of Coil 2 : (0,0,0), 
% Centre of coil 3 : (0,0,d23) , Centre of coil 4:(0,0,d23+d34)
% Coils lie in planes parallel to X-Y Plane

%% Length Scale for Calculation

% For our application, define a cube of size 'l' .We find magnitude and 
% direction of field in this space only.

                % Start - Required for certain special cases % 
                
% Length scale of importance l = c/w_gyro. l = 1.71e-3/B , where B in [T].
% If B = 1 T, l = 1.71 mm. So dividing above cube into smaller cubes of
% dimension = l will suffice.  

% B_max_1 = (mu*I1)./(2*R1); % Field at centre of ring 1 due to ring 1 in [T]
% B_max_2 = (mu*I2)./(2*R2); % Field at centre of ring 2 due to ring 2 in [T]
% B_max = B_max_1 + B_max_2; % Maximum Field Possible in the system. Not plausible
%                            % anywhere, but is a good estimate for finding length scale.                                                     
% l = 1.71e-3/B_max;         % Length Scale under consideration in [m] 

                % End - Required for certain special cases %

n = 128;    % Desired Resolution for smaller cube. Keep it even (preferably 2^x).
L = 1.5*d;  % Safe length in [m]. In magnetic mirror, particle is anyways confined between
            % the rings
l = 2*L/n;  % Shortest Distance between 2 points field is calculated in [m].

%% Coordinate Points where field is measured : 

x = zeros(n+1,n+1,n+1);
y = x;
z = x; 

for k = 1:n+1
    for j = 1:n+1
        for i= 1:n+1
            
            z(i,j,k) = L - ((i-1)*l);
            x(i,j,k) = -L + ((j-1)*l);
            y(i,j,k) = L - ((k-1)*l);
            
            if(z(i,j,k) == 0 )
                z(i,j,k) = 1e-7;
            end
            
            if(x(i,j,k) == 0 )
                x(i,j,k) = 1e-7;
            end
            
            if(y(i,j,k) == 0 )
                y(i,j,k) = 0;
            end
            
            
        end
    end
end

%% Calculation of Magnetic Field due to 3 Rings

Bx1 = zeros(n+1,n+1,n+1);   % X Component of Magnetic Field of Ring 1
By1 = Bx1;                  % Y Component of Magnetic Field of Ring 1
Bz1 = Bx1;                  % Z Component of Magnetic Field of Ring 1

Bx2 = zeros(n+1,n+1,n+1);   % X Component of Magnetic Field of Ring 2
By2 = Bx2;                  % Y Component of Magnetic Field of Ring 2
Bz2 = Bx2;                  % Z Component of Magnetic Field of Ring 2

Bx3 = zeros(n+1,n+1,n+1);   % X Component of Magnetic Field of Ring 3
By3 = Bx3;                  % Y Component of Magnetic Field of Ring 3
Bz3 = Bx3;                  % Z Component of Magnetic Field of Ring 3

Bx4 = zeros(n+1,n+1,n+1);   % X Component of Magnetic Field of Ring 3
By4 = Bx4;                  % Y Component of Magnetic Field of Ring 3
Bz4 = Bx4;                  % Z Component of Magnetic Field of Ring 3

Bx_tot = zeros(n+1,n+1,n+1); % X Component of Total Magnetic Field 
By_tot = Bx_tot;             % Y Component of Total Magnetic Field
Bz_tot = Bx_tot;             % Z Component of Total Magnetic Field

% Magnetic Field Formulas : Solving by Biot-Savart Law will yield Bx,By and
% Bz upto an integral equation similar to elliptic integrals. The numerically solved
% and simplified form of those integrals can be found in the paper "Simple 
% Analytic Expressions for the Magnetic Field of a Circular Current Loop" by James Simpson et al. 
% Link : https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20010038494.pdf

% Formula in above paper is given assuming ring centre is at origin. Shift
% origin and compute input variables for B for both rings separately.

% Loop to Calculate B
for k = 1:n+1
    for j = 1:n+1
        for i= 1:n+1
            
            % B due to Ring 1 [Centre at (0,0,-d12) => (x,y,z) -> (x,y,z + d12)]
              z0 = -1*d12;
              C = mu*I1/pi ;  
              rho_sq = x(i,j,k).^2 + y(i,j,k).^2;
              rho = sqrt(rho_sq);
              
              r_sq = x(i,j,k).^2 + y(i,j,k).^2 + (z(i,j,k)- z0).^2;
              alpha_sq = R1.^2 + r_sq - (2*R1*rho);
              beta_sq  = R1.^2 + r_sq + (2*R1*rho);
              beta = sqrt(beta_sq);
              
              k_sq = 1 - (alpha_sq/beta_sq);
              %gamma = x(i,j,k).^2 - y(i,j,k).^2;
              [K,E] = ellipke(k_sq);
              
              constant = C/(2*alpha_sq*beta);
              
              Bz1(i,j,k) = constant * (((R1^2 - r_sq)*E) + (alpha_sq*K));
              Bx1(i,j,k) = constant * ((x(i,j,k).*(z(i,j,k) - z0))/(rho_sq)) * (((R1^2 + r_sq)*E) - (alpha_sq*K));
              By1(i,j,k) = constant * ((y(i,j,k).*(z(i,j,k) - z0))/(rho_sq)) * (((R1^2 + r_sq)*E) - (alpha_sq*K));
             
%               if(isnan(Bx1(i,j,k))||isnan(By1(i,j,k)))
%                   Bx1(i,j,k) = (3*(R1^2)*mu*I1*x(i,j,k)*z(i,j,k))/(4*(R1^2+z(i,j,k)^2)^2.5);
%                   By1(i,j,k) = (3*(R1^2)*mu*I1*y(i,j,k)*z(i,j,k))/(4*(R1^2+z(i,j,k)^2)^2.5);
%               end
                  
              % B due to Ring 2 [ Centre at (0,0,0) => (x,y,z) -> (x,y,z) ]
              z0 = 0;
              C = mu*I2/pi ;  
              rho_sq = x(i,j,k).^2 + y(i,j,k).^2;
              rho = sqrt(rho_sq);
              
              r_sq = x(i,j,k).^2 + y(i,j,k).^2 + (z(i,j,k)- z0).^2;
              alpha_sq = R2.^2 + r_sq - (2*R2*rho);
              beta_sq  = R2.^2 + r_sq + (2*R2*rho);
              beta = sqrt(beta_sq);
              
              k_sq = 1 - (alpha_sq/beta_sq);
              %gamma = x(i,j,k).^2 - y(i,j,k).^2;
              [K,E] = ellipke(k_sq);
              
              constant = C/(2*alpha_sq*beta);
              
              Bz2(i,j,k) = constant * (((R2^2 - r_sq)*E) + (alpha_sq*K));
              Bx2(i,j,k) = constant * ((x(i,j,k).*(z(i,j,k) - z0))/(rho_sq)) * (((R2^2 + r_sq)*E) - (alpha_sq*K));
              By2(i,j,k) = constant * ((y(i,j,k).*(z(i,j,k) - z0))/(rho_sq)) * (((R2^2 + r_sq)*E) - (alpha_sq*K));
             
%               if(isnan(Bx1(i,j,k))||isnan(By1(i,j,k)))
%                   Bx1(i,j,k) = (3*(R1^2)*mu*I1*x(i,j,k)*z(i,j,k))/(4*(R1^2+z(i,j,k)^2)^2.5);
%                   By1(i,j,k) = (3*(R1^2)*mu*I1*y(i,j,k)*z(i,j,k))/(4*(R1^2+z(i,j,k)^2)^2.5);
%               end

              % B due to Ring 3 [Centre at (0,0,d23) => (x,y,z) -> (x,y,z - d23)]
              
              z0 = +1*d23;
              C = mu*I3/pi ;  
              rho_sq = x(i,j,k).^2 + y(i,j,k).^2;
              rho = sqrt(rho_sq);
              
              r_sq = x(i,j,k).^2 + y(i,j,k).^2 + (z(i,j,k)- z0).^2;
              alpha_sq = R3.^2 + r_sq - (2*R3*rho);
              beta_sq  = R3.^2 + r_sq + (2*R3*rho);
              beta = sqrt(beta_sq);
              
              k_sq = 1 - (alpha_sq/beta_sq);
              %gamma = x(i,j,k).^2 - y(i,j,k).^2;
              [K,E] = ellipke(k_sq);
              
              constant = C/(2*alpha_sq*beta);
              
              Bz3(i,j,k) = constant * (((R3^2 - r_sq)*E) + (alpha_sq*K));
              Bx3(i,j,k) = constant * ((x(i,j,k).*(z(i,j,k) - z0))/(rho_sq)) * (((R3^2 + r_sq)*E) - (alpha_sq*K));
              By3(i,j,k) = constant * ((y(i,j,k).*(z(i,j,k) - z0))/(rho_sq)) * (((R3^2 + r_sq)*E) - (alpha_sq*K));
             
%               if(isnan(Bx1(i,j,k))||isnan(By1(i,j,k)))
%                   Bx1(i,j,k) = (3*(R1^2)*mu*I1*x(i,j,k)*z(i,j,k))/(4*(R1^2+z(i,j,k)^2)^2.5);
%                   By1(i,j,k) = (3*(R1^2)*mu*I1*y(i,j,k)*z(i,j,k))/(4*(R1^2+z(i,j,k)^2)^2.5);
%               end
              
              % B due to Ring 4 [Centre at (0,0,d23+d34) => (x,y,z) -> (x,y,z - (d23+d34))]
              
              z0 = +1*(d23+d34);
              C = mu*I4/pi ;  
              rho_sq = x(i,j,k).^2 + y(i,j,k).^2;
              rho = sqrt(rho_sq);
              
              r_sq = x(i,j,k).^2 + y(i,j,k).^2 + (z(i,j,k)- z0).^2;
              alpha_sq = R4.^2 + r_sq - (2*R4*rho);
              beta_sq  = R4.^2 + r_sq + (2*R4*rho);
              beta = sqrt(beta_sq);
              
              k_sq = 1 - (alpha_sq/beta_sq);
              %gamma = x(i,j,k).^2 - y(i,j,k).^2;
              [K,E] = ellipke(k_sq);
              
              constant = C/(2*alpha_sq*beta);
              
              Bz4(i,j,k) = constant * (((R4^2 - r_sq)*E) + (alpha_sq*K));
              Bx4(i,j,k) = constant * ((x(i,j,k).*(z(i,j,k) - z0))/(rho_sq)) * (((R4^2 + r_sq)*E) - (alpha_sq*K));
              By4(i,j,k) = constant * ((y(i,j,k).*(z(i,j,k) - z0))/(rho_sq)) * (((R4^2 + r_sq)*E) - (alpha_sq*K));
             
%               if(isnan(Bx1(i,j,k))||isnan(By1(i,j,k)))
%                   Bx1(i,j,k) = (3*(R1^2)*mu*I1*x(i,j,k)*z(i,j,k))/(4*(R1^2+z(i,j,k)^2)^2.5);
%                   By1(i,j,k) = (3*(R1^2)*mu*I1*y(i,j,k)*z(i,j,k))/(4*(R1^2+z(i,j,k)^2)^2.5);
%               end

              % Total Magnetic Field
              Bx_tot(i,j,k) = Bx1(i,j,k) + Bx2(i,j,k)+ Bx3(i,j,k) + Bx4(i,j,k);
              By_tot(i,j,k) = By1(i,j,k) + By2(i,j,k)+ By3(i,j,k) + By4(i,j,k);
              Bz_tot(i,j,k) = Bz1(i,j,k) + Bz2(i,j,k)+ Bz3(i,j,k) + Bz4(i,j,k);
              
        end
    end
end

B_mag = sqrt(Bx_tot.^2 + By_tot.^2 + Bz_tot.^2);

%% Plots
[X1,Z1]= meshgrid(linspace(-L,L,n+1),linspace(-L,L,n+1));

%Plot 1 - 4 : Magnetic Field in y = 0 plane

% Note : B_mag(:,:,(n/2)+1)
% i)   The indice (n/2)+1 represents the mid plane. The y = 0 plane. As a
%      test see they 'By' plot. It will be blank. That is true, even
%      according to theory.
% ii)  The indice n/2 represents the plane just above (i.e) y = l plane (See Line 55)
% iii) The indice '1' will represent the top most plane (i.e) y = L plane
% Adjust the indice appropriately to see the different planes.

figure(1)
surf(X1,Z1,B_mag(:,:,(n/2)+1)) %(n/2)
xlabel('X [m]','FontWeight','bold')
ylabel('Z [m]','FontWeight','bold')
zlabel('B Total [T]','FontWeight','bold')
title({'Total Magnetic Field Due to 4 Coil System in y = 0 Plane'})
colorbar %add colorbar
set(gca, 'ZScale', 'log')
set(gca,'colorscale','log')
shading flat %Removes black lines from the mesh

figure(2)
surf(X1,Z1,abs(Bx_tot(:,:,(n/2)+1)))
xlabel('X [m]','FontWeight','bold')
ylabel('Z [m]','FontWeight','bold')
zlabel('Bx [T]','FontWeight','bold')
title({'Bx Magnetic Field Due to 4 Coil System in y = 0 Plane'})
colorbar %add colorbar
set(gca, 'ZScale', 'log')
set(gca,'colorscale','log')
shading flat %Removes black lines from the mesh

figure(3)
surf(X1,Z1,abs(By_tot(:,:,(n/2)+1)))
xlabel('X [m]','FontWeight','bold')
ylabel('Z [m]','FontWeight','bold')
zlabel('By [T]','FontWeight','bold')
title({'By Magnetic Field Due to 4 Coil System in y = 0 Plane'})
colorbar %add colorbar
set(gca, 'ZScale', 'log')
set(gca,'colorscale','log')
shading flat %Removes black lines from the mesh

figure(4)
surf(X1,Z1,abs(Bz_tot(:,:,(n/2)+1)))
xlabel('X [m]','FontWeight','bold')
ylabel('Z [m]','FontWeight','bold')
zlabel('Bz [T]','FontWeight','bold')
title({'Bz Magnetic Field Due to 4 Coil System in y = 0 Plane'})
colorbar %add colorbar
set(gca, 'ZScale', 'log')
set(gca,'colorscale','log')
shading flat %Removes black lines from the mesh

%%
% figure(5)
% axes5 = axes('Parent',figure(5));
% surf(X1,Z1,B_mag(:,:,(n/2))) %(n/2)
% xlabel('X [m]','FontWeight','bold')
% ylabel('Z [m]','FontWeight','bold')
% zlabel('B Total [T]','FontWeight','bold')
% title({'Total Magnetic Field Due to 4 Coil System in y = 1.17 cm Plane'})
% colorbar %add colorbar
% set(gca, 'ZScale', 'log')
% set(gca,'colorscale','log')
% shading flat %Removes black lines from the mesh
% view(axes5,[0.900000000000006 90]);
% 
% figure(6)
% axes6 = axes('Parent',figure(6));
% surf(X1,Z1,abs(Bx_tot(:,:,(n/2))))
% xlabel('X [m]','FontWeight','bold')
% ylabel('Z [m]','FontWeight','bold')
% zlabel('Bx [T]','FontWeight','bold')
% title({'Bx Magnetic Field Due to 4 Coil System in y = 1.17 cm Plane'})
% colorbar %add colorbar
% set(gca, 'ZScale', 'log')
% set(gca,'colorscale','log')
% shading flat %Removes black lines from the mesh
% view(axes6,[0.900000000000006 90]);
% 
% 
% figure(7)
% axes7 = axes('Parent',figure(7));
% surf(X1,Z1,abs(By_tot(:,:,(n/2))))
% xlabel('X [m]','FontWeight','bold')
% ylabel('Z [m]','FontWeight','bold')
% zlabel('By [T]','FontWeight','bold')
% title({'By Magnetic Field Due to 4 Coil System in y = 1.17 cm Plane'})
% colorbar %add colorbar
% set(gca, 'ZScale', 'log')
% set(gca,'colorscale','log')
% shading flat %Removes black lines from the mesh
% view(axes7,[0.900000000000006 90]);
% 
% figure(8)
% axes8 = axes('Parent',figure(8));
% surf(X1,Z1,abs(Bz_tot(:,:,(n/2))))
% xlabel('X [m]','FontWeight','bold')
% ylabel('Z [m]','FontWeight','bold')
% zlabel('Bz [T]','FontWeight','bold')
% title({'Bz Magnetic Field Due to 4 Coil System in y = 1.17 cm Plane'})
% colorbar %add colorbar
% set(gca, 'ZScale', 'log')
% set(gca,'colorscale','log')
% shading flat %Removes black lines from the mesh
% view(axes8,[0.900000000000006 90]);

%% 
% figure(9)
% axes9 = axes('Parent',figure(9));
% surf(X1,Z1,B_mag(:,:,1)) %(n/2)
% xlabel('X [m]','FontWeight','bold')
% ylabel('Z [m]','FontWeight','bold')
% zlabel('B Total [T]','FontWeight','bold')
% title({'Total Magnetic Field Due to 4 Coil System in y = 1.5 m Plane'})
% colorbar %add colorbar
% set(gca, 'ZScale', 'log')
% set(gca,'colorscale','log')
% shading flat %Removes black lines from the mesh
% view(axes9,[0.900000000000006 90]);
% 
% figure(10)
% axes10 = axes('Parent',figure(10));
% surf(X1,Z1,abs(Bx_tot(:,:,1)))
% xlabel('X [m]','FontWeight','bold')
% ylabel('Z [m]','FontWeight','bold')
% zlabel('Bx [T]','FontWeight','bold')
% title({'Bx Magnetic Field Due to 4 Coil System in y = 1.5 m Plane'})
% colorbar %add colorbar
% set(gca, 'ZScale', 'log')
% set(gca,'colorscale','log')
% shading flat %Removes black lines from the mesh
% view(axes10,[0.900000000000006 90]);
% 
% 
% figure(11)
% axes11 = axes('Parent',figure(11));
% surf(X1,Z1,abs(By_tot(:,:,1)))
% xlabel('X [m]','FontWeight','bold')
% ylabel('Z [m]','FontWeight','bold')
% zlabel('By [T]','FontWeight','bold')
% title({'By Magnetic Field Due to 4 Coil System in y = 1.5 m Plane'})
% colorbar %add colorbar
% set(gca, 'ZScale', 'log')
% set(gca,'colorscale','log')
% shading flat %Removes black lines from the mesh
% view(axes11,[0.900000000000006 90]);
% 
% figure(12)
% axes12 = axes('Parent',figure(12));
% surf(X1,Z1,abs(Bz_tot(:,:,1)))
% xlabel('X [m]','FontWeight','bold')
% ylabel('Z [m]','FontWeight','bold')
% zlabel('Bz [T]','FontWeight','bold')
% title({'Bz Magnetic Field Due to 4 Coil System in y = 1.5 m Plane'})
% colorbar %add colorbar
% set(gca, 'ZScale', 'log')
% set(gca,'colorscale','log')
% shading flat %Removes black lines from the mesh
% view(axes12,[0.900000000000006 90]);

%% 

% %% Some Mathematical Stuff to test approximations made by theory
% 
% % Calculating d(Bz)/dz in Y = 0 Plane
% 
% dz = l; 
% dBz = zeros(n+1,n+1,n+1);
% der_Bz_by_z = dBz;
% 
% for k = 1:n+1
%     for j = 1:n+1
%         for i= 1:n+1
%             
%             if(i == 1)
%                 
%                 dBz(i,j,k) = Bz_tot(i+1,j,k) - Bz_tot(i,j,k);
%                 der_Bz_by_z(i,j,k) = dBz(i,j,k)/dz;
%             
%             end
%             
%             if(i >= 2 && i < n+1)
%               
%                 dBz(i,j,k) = Bz_tot(i+1,j,k) - Bz_tot(i-1,j,k);
%                 der_Bz_by_z(i,j,k) = dBz(i,j,k)/(2*dz);
%                 
%             end
%             
%             if (i == n+1)
%                 
%                 dBz(i,j,k) = 0;
%                 der_Bz_by_z(i,j,k) = dBz(i,j,k)/(2*dz);
%                 
%             end
%                 
%                 
%             
%         end
%     end
% end
%             
% %% Plots 
% 
% figure(5)
% surf(X1,Z1,der_Bz_by_z(:,:,1))
% xlabel('X [m]','FontWeight','bold')
% ylabel('Z [m]','FontWeight','bold')
% zlabel('d(Bz)/dz','FontWeight','bold')
% title({'d(Bz)/dz Due to 4 Coil System in y = 0 Plane'})
% colorbar %add colorbar
% shading flat %Removes black lines from the mesh


toc
