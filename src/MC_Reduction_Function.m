clear 
close all
%% Path sampling algorithm for phonon and electron mean free path
% Algorithm found from Penelope ray tracing manual; computational physics
% tutorial youtube
%Silicon nanobridges Marconnet et al.2012; Mcgauhey

% Allocating memory 

n_particles=10000;
x=zeros(10000,1);
z=zeros(10000,1);
r=zeros(10000,1);
OP=zeros(10000,1);
mfp_eff=zeros(10000,1);
red_fn=zeros(100,1);
%xi=zeros(50,1);
%phi=zeros(1000,1);
dx=zeros(10000,1);
dz=zeros(10000,1);
lambda_intrinsic= 150e-9;
lambda_boundary=logspace(0,5,1000).*2e-9;
xi=[rand();rand();rand();rand()];

for runs=1:1000
       
    lambda_iter=lambda_boundary(runs);
    
    for count=1:n_particles


        % Defining the boundary dimensions

        rectangleWidth = 500e-9;
        rectangleHeight = 200e-9;
    
        is_scattered=0;
        iter_count=0;


        while is_scattered==0 && iter_count<1000
            
           % Initializing the source of phonon launch co-ordinates
           % For an assumption that only the MFP in z-direction varies
           
           x(count)=xi(1)*rectangleWidth;
           z(count)=xi(2)*rectangleHeight;

            % probability distributions of theta (polar) and (phi)azimuth angles

            theta =acos(1-2*xi(3)); % alternatively try asin(-1+2*xi)
            phi=2*pi*xi(4);
            
           % Computing the distance from the starting position to the
           % nearest boundary
            distancesToTop = z(count);
            distancesToBottom = rectangleHeight - z(count);
            distancesToLeft = x(count);
            distancesToRight = rectangleWidth - x(count);
            delta_x=(distancesToRight-distancesToLeft)^2;
            delta_z=(distancesToTop-distancesToBottom)^2;
            

            % Identify the nearest boundary
            
            OP(count) = sqrt(delta_x+delta_z);

            % Final direction vectors of scattering
            dx(count)=x(count)*sin(theta)*cos(phi);
            dz(count)=z(count)*cos(phi);


            % Computing distance that a phonon travels
            x(count)=x(count)+dx(count);
            z(count)=z(count)+dz(count);              
            r(count)=sqrt(x(count)^2+z(count)^2);

            % Counting the boundary scattering 

            iter_count=iter_count+1;

            if r(count) <= OP(count)
                is_scattered=1;
            end

        end

        % Calculating the effective mean free path

        if r(count) < lambda_iter
            mfp_eff(count)=r(count);
        else
            mfp_eff(count)=lambda_iter;      
        end
        
           
               
    end
    
    % Computing the reduction function
    
    red_fn(runs)=(1/(n_particles*lambda_iter))*sum(mfp_eff);
   

end
% 
% % ravg=mean(r(count));
% % disp(ravg);
% % disp(tin+tout);d
del=rectangleHeight./lambda_boundary;
figure(2)
semilogx(del,red_fn,'LineWidth',2);
xlabel('H/gamma','FontSize',16)
ylabel('Reducation function(F)','FontSize',16)
