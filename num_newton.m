function num_newton()
%Adatok megadása

n = 100; h = 1/(n + 1);% osztópontok száma a tengelyek irányába, rácstávolság
xx = h*[1:n]; yy = xx;% oldalak felosztása
[x, y] = meshgrid(xx, yy); % a rácspontok koordinátáinak eltárolása
u_h = zeros(n); % 5*sin(pi*x).*sin(pi*y) zeros(n)
u_h = reshape(u_h, n^2, 1);              % a kezdõvektor
u_temp = ones(n^2, 1);% segédváltozó; mivel az iterációban u_h értékeit
      % frissítjük, ezért az n. megoldást ebben tároljuk
f = 0.9*ones(n);  % az egyenlet jobb oldala
f = reshape(f, n^2, 1);

% A diszkretizációs mátrix konstrukciója
r_egyseg = speye(n);
E = sparse(2:n,1:n-1,1,n,n);
D = E+E'-2*r_egyseg;
A = -kron(D,r_egyseg)-kron(r_egyseg,D);
log_A = logical(A);
num_iter = 0;
% Indítjuk a stopperórát
rent = tic;
while  norm(u_h - u_temp, inf) > 1e-10 %  Az iteráció
    u_temp = u_h;
    % Az alacsonyabbrendû tagok beépítése a diszkretizációs mátrixba
    em = exp(u_temp)';
    %eu = em(ones(n^2, 1), :);%  repmat(em, n^2, 1);
    A_temp = A + (h^2 / 3) * (log_A .* em(ones(n^2, 1), :));
    % az egyenlet jobb oldalának megadása
    j_o = (exp(u_temp).*(u_temp - 1) + f) * h^2;
    % Az egyenletrendszer megoldása
    u_h = A_temp \ j_o;
    num_iter = num_iter + 1;
    %sum_time = sum_time + sum(T);
end
renti= toc(rent);
disp(renti);
% ebben a változóban tároljuk az iteráció sebességét
disp(num2str(num_iter));
u_grid_in = reshape(u_h, n, n); % megoldásvektor visszalakítása mátrixsszá 
u_grid = zeros(n + 2, n + 2); % itt vesszük hozzá a peremfeltételt
u_grid(2:n + 1, 2:n + 1)=u_grid_in;

% A háromszögek konstrukciója a hozzájuk tartozó csúcsokból.
triangular_mesh=[0 0 0];
for i=1:(n + 1)*(n + 2)
    if mod(i, n + 2)~=0
       triangular_mesh(length(triangular_mesh(:, 1))+ 1,:)=[i, i + 1, i + n + 2];
       triangular_mesh(length(triangular_mesh(:,1))+1,:)=[i + 1, i + 1 + n + 2, i + n + 2];
    end
end
triangular_mesh(1,:)=[];
[x_0, y_0]=meshgrid(0:h:1,0:h:1);
[u_x, u_y] = gradient(u_grid, h);
% Az eredmények ábrázolása


subplot(1, 2, 1)
% A végeselemes megoldás ábrázolása
trisurf(triangular_mesh, x_0, y_0, u_grid)
axis([0, 1, 0, 1, -0.01, 0])
xlabel('x','FontSize',12)
ylabel('y','FontSize',12)
title(['Az egyenlet Newton-módszerrel adott megoldása.'],'FontSize',12)

subplot(1, 2, 2)
% A megoldás gradiensének ábrázolása
quiver(x_0, y_0, u_x, u_y)
axis([0, 1, 0, 1])
xlabel('x','FontSize',12)
ylabel('y','FontSize',12)
title(['A megoldás gradiensének vektormezeje.'], 'FontSize', 12)

%Az iteráció idõtartama
%X = ['Az iteráció ', num2str(T), ' másodperc alatt konvergált.'];
%disp(X)