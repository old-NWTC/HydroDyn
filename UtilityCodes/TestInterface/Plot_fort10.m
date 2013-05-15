
fileName = 'D:\DATA\Fortran\IVF Projects\HydroDyn\TestInterface\fort.10';

data = dlmread(fileName);

inpData = data(:,1:3);
invData = data(:,4:6);
eulData = data(:,7:9);

eulErr = inpData - eulData;
relEulErr = eulErr ./ inpData;
relEulErr(~isfinite(relEulErr)) = NaN;

disp( ['Max error on small angle inversion = ', num2str(norm( invData(:)  - inpData(:), inf )) ' radians.']);

disp( ['Max error on euler angle 123 inversion = ', num2str(norm( eulErr(:),inf )) ' radians.']);
disp( ['Max error on euler angle 123 inversion = ', num2str(max( abs(relEulErr) )) '%.']);

%% let's plot the errors:

theta1 = inpData(:,1);
theta2 = inpData(:,2);
theta3 = inpData(:,3);

t1 = unique(theta1);
t2 = unique(theta2);
t3 = unique(theta3);

n1 = length(t1);
n2 = length(t2);
n3 = length(t3);

theta1 = reshape(theta1     ,n3,n2,n1);
theta2 = reshape(theta2     ,n3,n2,n1);
theta3 = reshape(theta3     ,n3,n2,n1);
err{3} = reshape(eulErr(:,3),n3,n2,n1);
err{2} = reshape(eulErr(:,2),n3,n2,n1);
err{1} = reshape(eulErr(:,1),n3,n2,n1);

f(3)=figure;
f(2)=figure;
f(1)=figure;
C = ones(n3,n2,1);
lgndTxt = cell(n1,1);

for i_n1 = 1:2:n1
   
    lgndTxt{i_n1} = ['\theta_1 = ' num2str(t1(i_n1)) ' radians'];
    for ierr = 1:3
        figure(f(ierr));        
            
        surf(theta2(:,:,i_n1), theta3(:,:,i_n1), err{ierr}(:,:,i_n1),C*t1(i_n1) );
        hold on;
    end
end
%
for ierr = 1:3
    figure(f(ierr));   
    title(['Euler Angle \theta_' num2str(ierr) ' Error (radians)'])
    zlabel('Error (radians)');
    ylabel('\theta_3 (radians)');
    xlabel('\theta_2 (radians)');
    
    legend(lgndTxt(1:2:n1));
end
