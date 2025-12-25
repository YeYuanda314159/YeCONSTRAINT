function [sigma_square] = solver_elasticity_correction_Q1(xPhys,lambda,xi,C)

global F H freedofs nelx nely E0 Emin  mu lambdae 
%% 单位单元刚度矩阵
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
AA = [A11 A12;A12' A11];
BB = [B11 B12;B12' B11]; 
ak = 1/(1+lambda)*(1./(1/Emin+xPhys*(1/E0-1/Emin))).*max(xi-C,0).*xPhys;
A = 1./(1+2*mu*ak);
B = lambdae.*ak./A./((3*lambda+2*mu)*ak+1);
%% PREPARE FINITE ELEMENT ANALYSIS
% Ae = 2\mu e + \lambda tr(e) I
% \mu = E /2/ (1+nu) \lambda = E*nu/(1+nu)/(1+nu*(1-d))
%matrix_dim = 8*8
nodenrs = reshape(1:(1+nelx)*(1+nely), 1+nely, 1+nelx); %按列取数据排成(1+nely)行,(1+nelx)列
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1, nelx*nely, 1); %按列取数据排成(nelx*mely)行,1列
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1], nelx*nely, 1); %给出每个元的节点的编号
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);%遍历每个单元依次经过的节点编号
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);%用于装配刚度矩阵
%% 求解
U = zeros(2*(nely + 1)*(nelx + 1), 1); %约束的解向量
%V = zeros(2*(nely + 1)*(nelx + 1), 1); %目标函数对应的解向量
sK = reshape(AA(:)*A(:)+BB(:)*B(:), 64*nelx*nely, 1);%混合问题
%重排前的矩阵每一列表示一个物理单元的刚度矩阵
%xPhys(:),表示按列取数据将矩阵变成一列向量
K = sparse(iK,jK,sK); K = (K+K')/2; %装配整体刚度矩阵并保证对称性
U(freedofs) = K(freedofs,freedofs)\F(freedofs); %U是状态的位移向量
%V(freedofs) = K(freedofs,freedofs)\H(freedofs); %V是伴随的位移向量
  
%% 计算相关函数和泛函值
AU = reshape(sum((U(edofMat)*AA).*ones(nelx*nely,1),2),nely,nelx);
BU = reshape(sum((U(edofMat)*BB).*ones(nelx*nely,1),2),nely,nelx);
sigma_square = A.*AU+B.*BU; %柔度函数: 1/A(\chi)E_0e(u):e(v)
end