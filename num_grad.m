function num_grad()
% Adatok megad�sa

n = 100; h = 1/(n + 1); % oszt�pontok sz�ma a tengelyek ir�ny�ba, r�cst�vols�g
xx = h*[1:n]; yy = xx; % oldalak feloszt�sa
[x, y] = meshgrid(xx, yy); % a r�cspontok koordin�t�inak elt�rol�sa
u_h = zeros(n);            %  5*sin(pi*x).*sin(pi*y)
u_h = reshape(u_h, n^2, 1);             % a kezd�vektor az azonosan nulla
u_temp = ones(n^2, 1);% seg�dv�ltoz�; mivel az iter�ci�ban u_h �rt�keit
      % friss�tj�k, ez�rt az n. megold�st ebben t�roljuk
f= 0.9*ones(n);  % az egyenlet jobb oldala
f = reshape(f, n^2, 1);

% A diszkretiz�ci�s m�trix konstrukci�ja
I = speye(n);
E = sparse(2:n,1:n-1,1,n,n);
D = E+E'-2*I;
A = -kron(D,I)-kron(I,D);
num_iter = 0;
tic % ind�tjuk a be�p�tett stopper�r�t
while  norm(u_h - u_temp, inf) > 1e-10 %Az iter�ci�
    u_temp = u_h;
    b = (exp(u_temp) - f) * h^2; % az egyenlet jobb oldal�nak megad�sa
    w = A\b; % Az egyenletrendszer megold�sa
    % Az iter�ci�s l�p�s
    u_h = (u_temp - 2 * pi^2 * w) / (2*pi^2 + 1);
    num_iter = num_iter + 1;
end
run_time = toc; % ebben a v�ltoz�ban t�roljuk az iter�ci� sebess�g�t
disp(num_iter);
u_grid_in = reshape(u_h, n, n); % megold�svektor visszalak�t�sa m�trixssz� 
u_grid = zeros(n + 2, n + 2); % itt vessz�k hozz� a peremfelt�telt
u_grid(2:n + 1, 2:n + 1) = u_grid_in;

% A h�romsz�gek konstrukci�ja a hozz�juk tartoz� cs�csokb�l.
triangular_mesh = [0 0 0];
for i=1:(n + 1)*(n + 2)
    if mod(i, n + 2)~=0
           triangular_mesh(length(triangular_mesh(:, 1))+ 1,:)=[i, i + 1, i + n + 2];
           triangular_mesh(length(triangular_mesh(:,1))+1,:)=[i + 1, i + 1 + n + 2, i + n + 2];
    end
end
triangular_mesh(1,:) = [];
[x_0, y_0] = meshgrid(0:h:1,0:h:1);
[u_x, u_y] = gradient(u_grid, h);
% Az eredm�nyek �br�zol�sa


subplot(1, 2, 1)
% A v�geselemes megold�s �br�zol�sa
trisurf(triangular_mesh, x_0, y_0, u_grid)
axis([0, 1, 0, 1, -0.01, 0])
xlabel('x','FontSize',12)
ylabel('y','FontSize',12)
title(['Az egyenlet gradines-m�dszerrel adott megold�sa.'],'FontSize',12)

subplot(1, 2, 2)
% A megold�s gradiens�nek �br�zol�sa
quiver(x_0, y_0, u_x, u_y)
axis([0, 1, 0, 1])
xlabel('x','FontSize',12)
ylabel('y','FontSize',12)
title(['A megold�s gradiens�nek vektormezeje.'], 'FontSize', 12)

%Az iter�ci� id�tartama
X = ['Az iter�ci� ', num2str(run_time), ' m�sodperc alatt konverg�lt.'];
disp(X)