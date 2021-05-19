a = 1;
L = 10;
H = 5;
dt = 1e-3;
dx = 0.1;
dy = 0.1;
UR = 0;
UL = 0;
UT = 0;
UB = 0;
m = H/dy;
n = L/dx;
T0 = 0;

U = ones(m,n)*T0;
f = zeros(m, n);
U(:,1) = UL;
U(:,n) = UR;
U(1,:) = UT;
U(m,:) = UB;
%U(10,30) = 1;
U_pr = U;
U_array = zeros(1);

%f(3, 70:75) = 1000*sin(k);


dt_array = zeros(1);

figure('Color','w')

U_n = U;

for k = 1:3:100000
    if mod(k-1,150) ==0
      clf
subplot(2,1,1)
surf(U)
axis([0 10 5 35])
title("U(k)")
xlabel('L')
ylabel('H')
colormap('winter')
shading interp
colorbar
view(0,90)

axis equal
subplot(2,1,2)
U_sum = sum(U(:));
dt_t=dt*k;

U_array(end+1) =  U_sum/(m*n); 
dt_array(end+1)= dt_t;

plot(dt_array,U_array)
title("Зависимость среднего смещения от времени")
xlabel('t')
ylabel('U')
pause(0.1)
    end
    for i = 2:(H/dy-1)
        for j = 2:(L/dx-1)
          U_n(i,j) = 2*U(i,j) - U_pr(i,j) + a*dt^2*((U(i+1,j) - 2*U(i,j)+ U(i-1,j))/dx^2+(U(i, j+1)-2*U(i,j)+ U(i, j-1))/dy^2)+ f(i, j)*dt^2;  
        end
    end
    U_pr = U;
    U = U_n;
%     U(H/dy/2, L/dx/2) = 1;
U(m*0.7,47:58) = 10*sin((k-1)/500);

 U(1,:) = 0.1;
 U(1:n) = 0.1;
 U(:, n) = 0.1;

end
