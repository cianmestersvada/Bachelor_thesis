function num_newton()
%Adatok megad�sa

n = 100; h = 1/(n + 1);% oszt�pontok sz�ma a tengelyek ir�ny�ba, r�cst�vols�g
xx = h*[1:n]; yy = xx;% oldalak feloszt�sa
[x, y] = meshgrid(xx, yy); % a r�cspontok koordin�t�inak elt�rol�sa
u_h = zeros(n); % 5*sin(pi*x).*sin(pi*y) zeros(n)
u_h = reshape(u_h, n^2, 1);              % a kezd�vektor
u_temp = ones(n^2, 1);% seg�dv�ltoz�; mivel az iter�ci�ban u_h �rt�keit
      % friss�tj�k, ez�rt az n. megold�st ebben t�roljuk
f = 0.9*ones(n);  % az egyenlet jobb oldala
f = reshape(f, n^2, 1);

% A diszkretiz�ci�s m�trix konstrukci�ja
r_egyseg = speye(n);
E = sparse(2:n,1:n-1,1,n,n);
D = E+E'-2*r_egyseg;
A = -kron(D,r_egyseg)-kron(r_egyseg,D);
log_A = logical(A);
num_iter = 0;
% Ind�tjuk a stopper�r�t
rent = tic;
while  norm(u_h - u_temp, inf) > 1e-10 %  Az iter�ci�
    u_temp = u_h;
    % Az alacsonyabbrend� tagok be�p�t�se a diszkretiz�ci�s m�trixba
    em = exp(u_temp)';
    %eu = em(ones(n^2, 1), :);%  repmat(em, n^2, 1);
    A_temp = A + (h^2 / 3) * (log_A .* em(ones(n^2, 1), :));
    % az egyenlet jobb oldal�nak megad�sa
    j_o = (exp(u_temp).*(u_temp - 1) + f) * h^2;
    % Az egyenletrendszer megold�sa
    u_h = A_temp \ j_o;
    num_iter = num_iter + 1;
    %sum_time = sum_time + sum(T);
end
renti= toc(rent);
disp(renti);
% ebben a v�ltoz�ban t�roljuk az iter�ci� sebess�g�t
disp(num2str(num_iter));
u_grid_in = reshape(u_h, n, n); % megold�svektor visszalak�t�sa m�trixssz� 
u_grid = zeros(n + 2, n + 2); % itt vessz�k hozz� a peremfelt�telt
u_grid(2:n + 1, 2:n + 1)=u_grid_in;

% A h�romsz�gek konstrukci�ja a hozz�juk tartoz� cs�csokb�l.
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
% Az eredm�nyek �br�zol�sa


subplot(1, 2, 1)
% A v�geselemes megold�s �br�zol�sa
trisurf(triangular_mesh, x_0, y_0, u_grid)
axis([0, 1, 0, 1, -0.01, 0])
xlabel('x','FontSize',12)
ylabel('y','FontSize',12)
title(['Az egyenlet Newton-m�dszerrel adott megold�sa.'],'FontSize',12)

subplot(1, 2, 2)
% A megold�s gradiens�nek �br�zol�sa
quiver(x_0, y_0, u_x, u_y)
axis([0, 1, 0, 1])
xlabel('x','FontSize',12)
ylabel('y','FontSize',12)
title(['A megold�s gradiens�nek vektormezeje.'], 'FontSize', 12)

%Az iter�ci� id�tartama
%X = ['Az iter�ci� ', num2str(T), ' m�sodperc alatt konverg�lt.'];
%disp(X)