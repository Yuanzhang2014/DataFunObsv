function [Phi_x, Phi_y, centers_x, centers_y] = lift_koopman(x, y, num_centers_x, num_centers_y, num_steps, centers_x, centers_y)
% 新增输入参数centers_x, centers_y，训练阶段不传入，测试阶段传入训练得到的中心点
% 新增输出参数centers_x, centers_y，训练阶段保存后用于测试阶段


% 定义薄板样条径向基函数:∥x − x0 ∥2 log(∥x − x0 ∥)
%thin_plate_spline = @(x, center0) norm(x - center0)^2 * log(max(norm(x - center0), eps));
thin_plate_spline = @(x, c) max(norm(x-c), 1e-3)^2 * log(max(norm(x-c), 1e-3));


centers_x = 1.5*lhsdesign(num_centers_x,2) - 0.75;
% 创建一个包含 100 个基函数x的函数句柄数组
basis_functions_x = cell(num_centers_x, 1);
for i = 1:num_centers_x
     center0_x = centers_x(i, :);
    basis_functions_x{i} = @(x) thin_plate_spline(x, center0_x);
end

centers_y = 1.5*lhsdesign(num_centers_y,2) - 0.75;
% 创建一个包含 100 个基函数y的函数句柄数组
basis_functions = cell(num_centers_x, 1);
for i = 1:num_centers_y
     center0_y = centers_y(i, :);
    basis_functions_y{i} = @(x) thin_plate_spline(x, center0_y);
end

phix = [];
phiy = [];
% x提升维度，基函数phix生成
for i = 1:num_centers_x %num_center个函数
    for j = 1:num_steps
        phix(i,j) = basis_functions_x{i}(x(:,j));
    end
end

%y提升维度，基函数phiy生成
for i = 1:num_centers_y %100个函数
    for j = 1:num_steps
        phiy(i,j) = basis_functions_y{i}(y(:,j));
    end
end

%定义phix，phiy
    Phi_x = [x;phix];
    Phi_y = [y;phiy];
    
    %数据标准化
    Phi_x = (Phi_x - mean(Phi_x,2)) ./ std(Phi_x,0,2);
    Phi_y = (Phi_y - mean(Phi_y,2)) ./ std(Phi_y,0,2);

end
