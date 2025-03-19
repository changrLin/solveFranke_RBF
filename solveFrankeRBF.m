% 定义Franke函数
Franke = @(x,y) 0.75*exp(-((9*x-2).^2 + (9*y-2).^2)/4) + ...
              0.75*exp(-((9*x+1).^2)/49 - ((9*y+1).^2)/10) + ...
              0.5 *exp(-((9*x-7).^2 + (9*y-3).^2)/4) - ...
              0.2 *exp(-((9*x-4).^2 + (9*y-7).^2));
%取一组形状参数
ep = [0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0];
N = [10,20,30,50,70,90,110,130];
conj = zeros(length(N),1);
for m =1: length(N)
 %取10个插值点
    X_inter = linspace(0, 1, N(m));
    Y_inter = linspace(0,1, N(m));
%取20个点计算点
    x = linspace(0,1,20);
    y = linspace(0,1,20);
    A = [x',y'];
    coord = [X_inter',Y_inter'];
    F_values = Franke(X_inter,Y_inter);
    distances = pdist2(coord,coord);
%什么是径向基函数(RBF)插值技术？简单来说就是使用一组径向基的线性组合逼近一个难以计算的函数
%我们考虑两个径向基：Gauss 和 inverse quadratic
%-------Gauss kernel 
    RMES_1 = zeros(10,1);
    for j = 1:10
        temp = 0;
        basis_func_1 = exp(-(ep(j)*distances).^2);
        coefficients = pinv(basis_func_1)*F_values';
        X = zeros(length(X_inter),1);
        for k = 1:20
            for i =1:length(X_inter)
                X(i) = norm(A(k,:)-coord(i,:));
            end
        basis_func_1_val = exp(-(ep(j)*X).^2);
        f_val = coefficients'*basis_func_1_val;
        f_ture = Franke(A(k,1),A(k,2)); 
        t = (f_ture - f_val)^2;%误差
        temp = temp + t; %残差累计
        end
        RMES_1(j) = sqrt(temp/20);%均分根
    end

    disp(RMES_1);
    figure
    plot(ep,RMES_1,'r-*')
    title(['插值点数量 = ', num2str(N(m))]);
    
%------------inverse quadratic kernel
    RMES_2 = zeros(10,1);
    for j = 1:10
        temp = 0;
        basis_func_2 = 1./(1+(ep(j)*distances).^2);
        coefficients = pinv(basis_func_2)*F_values';
        X = zeros(length(X_inter),1);
        for k = 1:20
             for i =1:length(X_inter)
                 X(i) = norm(A(k,:)-coord(i,:));
             end
        basis_func_2_val = 1./(1+(ep(j)*X).^2);
        f_val = coefficients'*basis_func_2_val;
        f_ture = Franke(A(k,1),A(k,2)); 
        t = (f_ture - f_val)^2;
        temp = temp + t; 
        end
        RMES_2(j) = sqrt(temp/20);
    end
    hold on
    plot(ep,RMES_2,'b-*')
    legend('Guass kernel','IQ kernel')
    xlabel('形状参数ep之取值')
    ylabel('误差')
    hold off
    conj(m) = sum(RMES_1(:)-RMES_2(:));
end
figure
plot(N,conj,'b-*');
xlabel('插值点个数')
ylabel('两种径向核的均方差之差')