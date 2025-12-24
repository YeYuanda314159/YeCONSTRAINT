%%%% PENALTY METHOD USING IMGAUSSFILT %%%%
%** 最小单元尺度固定为1
function [y, loop, loop_k, c, x, energies, energies_k]= ...
    topthr_penlity_general(nelx, nely, lambda, r0, volfrac, frac, g, sd, objectfunc, bc, w,continuation, x, fileID,logtype, V_constrain, descent_type)
% nelx: number of elements on x axis
% nely: number of elements on y axis
% volfrac: volume fraction of material to total area
% frac: fraction of Emin to E0
% sd1: filter size for thresholding
% sd2: filter size for solving equation
% sd3: filter size for the perimeter term
% bc:  takes values in {'cantilever_rb', 'cantilever_central', 'mbb'},
%      enforces BCs (boundary condition)
% g: coefficient of the perimeter term
% continuation: 1 indicates to use x as the initial guess
%               0 indicates to use constant density as the initial guess
% x: initial guess, if continuaton is 0, can be anything.
% fileID: the opened file to log outputs
% V_constrain: 
% If you don't want to save the results, just let fileID to be negative 
% and the results will be displayed in the command window
% sameEmin: calculate the compliance each iteration using this specified value as Emin
% return c: compliance computed with filtered chi
%% OPTIMIZATION PARAMETERS
tol = 1;   %容许变化
maxtimes = 100;  %最大迭代次数
r = r0;
gamma1 = 3;
gamma2 = 1.5;
M = floor(nelx*nely*volfrac); %向下取整，\Omega1的元的数量
loop = 1;
loop_k = 1;
change = 10; %两次迭代中分布场xPhys的变化
gamma = g*sqrt(2*pi)/(sd*1/nely); %sd*1 = sd, 若想对分量求和应使用norm(sd,1) sd/nely = tau
energies = []; %存储每次迭代的目标函数值
energies_k = []; %存储每次搜索的目标函数值
%% MATERIAL PROPERTIES
E0 = 5000*8/3;  %设置杨氏模量
Emin = E0*frac; %设置人工材料的杨氏模量

%% DEFINE LOADS AND SUPPORTS assuming the mesh size is 1/nely
Fin = zeros(2*(nely+1)*(nelx+1),1);
Hout = zeros(2*(nely+1)*(nelx+1),1);
switch bc
    case 'left_bdc'
        fixeddofs = 1:2*(nely+1); %左边界
    case 'left_down_bdc_right_down_qin'
        force_length = floor(nelx/10); %1/10宽的边受力
        %sparse(参数1，参数2，参数3，参数4，参数5);
        %参数1填入行坐标，参数2填入列坐标，参数3填入值，参数4、5填入矩阵大小
        Fin((2*(nelx+1)*(nely+1)-2*force_length-1):2:(2*(nelx+1)*(nely+1)-1)) = -2; 
        fixeddofs = union(1 :2: 2*(nely+1)-1,(2*(nely+1)*(1:(nelx+1)))); 
        %fixeddofs = union(fixeddofs,2*(nely+1)*(1:(nelx+1))-1);
    case 'left_up_down_bdc_left_down_qin'
        force_length = floor(nely/10);
        Fin(2*(nely + 1-force_length)+1:2:(2*(nely+1)-1)) = 1;
        fixeddofs = union(1 : 2*force_length,(2*(nely+1)*(1:(nelx+1))));  
end
switch objectfunc %伴随问题的Dirichlet边界同元问题相同
    case 'right_down_central'
        %右边中点的单位位移
        Hout(2*(nely + 1) * nelx + 2*floor((nely+1)/2), 1) = -1;
    case 'left_up'
        force_length = floor(nelx/10);
        Hout(2+(2*(nely+1)*(0:force_length))) = -1; %H(2,1) = -1, 表示1号节点受到向下大小为1的点力
    case 'right_down'
        force_length = floor(nely/10); %1/8宽的边受力
        %下边界右侧1/8宽的边受到向下大小250的面力
        %sparse(参数1，参数2，参数3，参数4，参数5);
        %参数1填入行坐标，参数2填入列坐标，参数3填入值，参数4、5填入矩阵大小
        Hout((2*(nelx+1)*(nely+1)-2*force_length-1):2:(2*(nelx+1)*(nely+1)-1)) = 1; 
end
F = Fin + Hout;
H = Fin + w*Hout;
alldofs = [1:2*(nely+1)*(nelx+1)]; %所有自由度编号
freedofs = setdiff(alldofs,fixeddofs);  %Dirichlet 边界外的自由度编号
%% INITIALIZE ITERATION
if continuation == 0 %采用均匀密度作为初值
    x = repmat(volfrac,nely,nelx);
    xPhys = x; %matrix_dim = nely*nelx
else
    if sd > 0 %磨光因子大于零，可以进行磨光
        xPhys = imgaussfilt(x, sd); %基于内部启发式方法，在空间域或频域中执行卷积
        %重复边界元素填充矩阵
    else
        xPhys = x;
    end
end

%% START ITERATION
fid = fopen(fileID,logtype);
print_to_file= true;
try 
    fprintf(fid, 'Displaying\n'); %若此行执行失败则执行catch后的内容
    fprintf(fid, 'bc:%s | objectfunc:%s | Vforce:%3.5f | delta_x:%3.5f | volfrac:%3.5f\n',...
        bc,objectfunc,V_constrain,1/nely,volfrac);
    fprintf(fid, 'Emin:%5.5f | E0:%5.5f | gamma:%3f | tau:%5.5f | lambda:%5.5f | w:%5.5f | gamma1:%5.5f | gamma2:%5.5f\n',...
        Emin, E0, gamma, sd/nely, lambda,w, gamma1,gamma2);
catch err
    disp(err.message)
    disp('Now display the output in the command window.')
    print_to_file=false;
end
%创建固定的图像窗口
ax = axes('Parent', gcf);
while 1
    if change < tol  %只要变化大于0.01, 迭代不停
        fprintf('收敛原因2：自变量无法更新！\n');
        break;
    elseif loop > maxtimes %最多计算maxtimes次
        fprintf('收敛原因3：达到最大迭代次数！\n');
        break;
    end
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    if loop == 1
        [UV,c] = solver_elasticity_Q1(xPhys,F,H,freedofs,nelx,nely,E0,Emin);
        PG = g*sum(sum((1-x).*xPhys));
        energies(loop) = c+PG; %记录目标函数值
        energies_k(loop_k) = energies(loop);
    end
    
    %% Penalty Method--calculate Phi = (1/E0-1/Emin)G_\epsilon*(E_0^{-1}\sigma:\varrho)+ gamma/\epsilon*(x - xPhys)
    gk = 1./(1/Emin+xPhys*(1/E0-1/Emin)).*UV; %E_0^{-1}\sigma:\varrho = 1/A(\chi)^2 Ee(u):e(v)
    % filtering
    if sd > 0 
        gk = imgaussfilt(gk, sd, 'Padding', 'symmetric');
    end 
    gk = ((1/E0-1/Emin)*gk+ g*(x - xPhys));
    %共轭梯度法
    if loop == 1 && strcmp(descent_type, 'conjugate')
        gk0 = gk;
        Phi = gk;
    elseif loop > 1 && strcmp(descent_type, 'conjugate')
        betak = max(sum(sum(gk.*(gk-gk0)))/norm(gk0,2)^2,0);
        gk0 = gk;
        Phi = gk + betak*Phi;
    elseif strcmp(descent_type, 'gradient')
        Phi = gk;
    end
    if loop == 1
        r_min = 0;
    else
        %找出x中取值为1的索引
        mask = x > 0;
        indx1 = find(mask);
        indx0 = setdiff(1:numel(x), indx1); % 索引0（补集）

        % 最高效方法：逻辑索引
        max_val_indx1 = max(Phi(indx1));  % indx1中最大值
        min_val_indx0 = min(Phi(indx0));  % indx0中最小值
        if max_val_indx1 <= min_val_indx0 %满足一阶最优性条件，迭代终止
            fprintf('收敛原因1：达到一阶最优性条件！\n');
            break;
        end
        r_min = 1/(max_val_indx1 - min_val_indx0);
    end
    sorted_A = sort(Phi(:));  % 从大到小排序成向量
    min_spacing = min(diff(sorted_A)); % 相邻元素最小间距
    if r <= r_min + 1;
        r = min(1/min_spacing,r0);
        r_max = r;
    end
    %% Penalty Method--linear research
    while 1
        if loop_k > maxtimes
            change = 0;
            fprintf('线搜索结束4：超出最大迭代次数！\n');
            break;
        end
        loop_k = loop_k + 1;
        bar_Phi = x - r*(Phi + g*(x - xPhys)); %计算L^\infty中的 局部极小Phi
        [~,I] = sort(bar_Phi(:),'descend'); %由大到小快速排序
        % Project
        xnew = zeros(nelx*nely, 1);
        if (bar_Phi(I(M)) > 0.5) || (V_constrain == 0)
            xnew(I(1:M)) = 1; %最大的M个元是新的最优区域
        else
            xnew(bar_Phi > 0.5) = 1;
        end
        xnew = reshape(xnew, nely, nelx);
        %调整邻近因子
        if sd > 0
            xPhys_new = imgaussfilt(xnew, sd, 'Padding', 'symmetric');%	用自身的镜面反射填充图像。
        else
            xPhys_new = xnew;
        end
        [UV,c] = solver_elasticity_Q1(xPhys_new,F,H,freedofs,nelx,nely,E0,Emin); %计算试探的目标函数值
        PG = g*sum(sum((1-xnew).*xPhys_new)); %计算新的周长约束
        til_c = c + PG;
        energies_k(loop_k) = til_c;
        change = norm(xnew-x,1);%计算更新前后的区域的无穷范数
        %输出新图像
        cla(ax);  % 清除当前axes内容（保留axes设置）
        imshow(1-xnew, [], 'Parent', ax);  % 在同一axes显示新图像
        title(ax, sprintf('迭代 %d/%d', loop_k, maxtimes));
        colormap(ax, gray);
        drawnow;
        % PRINT RESULTS
        if print_to_file
            fprintf(fid, ' It.:%5i ||Obj.:%10.6f ||Vol.:%7.3f ||ch.:%7.3f ||r.:%7.3f \n', loop_k, til_c, ...
            mean(xPhys(:)),change,r); %mean: 计算xPhys的品均值
        else
            fprintf(' It.:%5i Obj.:%10.6f Vol.:%7.3f ch.:%7.3f r.%7.3f:\n', loop_k, til_c,...
            mean(xPhys(:)),change,r);
        end
        if til_c > energies(loop)
          r_max = r;
          if change < 2
              change = 0;
              fprintf('线搜索结束2：自变量不能更新！\n');
              break;
          end
        elseif (til_c < energies(loop)) % && ((energies(loop) - til_c)/energies(loop) > 0.001)
          x = xnew;
          xPhys = xPhys_new;
          loop = loop + 1;
          energies(loop) = til_c;
          fprintf('线搜索结束1：目标函数下降！\n');
          break;
        elseif til_c == energies(loop) && change > tol
            r_max = r;
        end
        if change == 0;
            r_min = r;
        end
        r = (r_max + r_min)/2;
        if abs(r_max-r_min) < 10^-3
            change = 0;
            fprintf('线搜索结束3：搜索区间过小！\n');
            break;
        end
    end
end
[~,c] = solver_elasticity_Q1(x,F,H,freedofs,nelx,nely,E0,Emin);
loop = loop + 1;
energies(loop) = c+g*(sum(sum((1-x).*xPhys))); %计算最后的目标函数值
%% FINAL OBJECTIVE FUNCTION WITHOUT SMOOTHING
y = energies(loop);
if print_to_file
    fprintf(fid, ' sharp interface energy: %11.6f\n', y);
    fprintf(fid, '--------------------------------------\n');
    fprintf(fid, '--------------------------------------\n');
else
    fprintf(' sharp interface energy: %11.8f\n', y);
end
fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written based on 88 line MATLAB code written by 
% E. Andreassen, A. Clausen, M. Schevenels,
% B. S. Lazarov and O. Sigmund,  Department of Solid  Mechanics,           %
% Technical University of Denmark,                                         %
% DK-2800 Lyngby, Denmark.                                                 %
%                                                                          %
% Disclaimer:                                                              %
% The author reserves all rights but do not guaranty that the code is      %
% free from errors. Furthermore, I shall not be liable in any event        %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%