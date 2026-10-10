function [vis_data, range_data] = isl_visibility_range(pos_data, cfg)
%ISL_VISIBILITY_RANGE  星间几何可见性（地球遮挡）与星间距离计算
%
%   [vis_data, range_data] = ISL_VISIBILITY_RANGE(pos_data, cfg)
%
%   输入：
%       pos_data : [N_slot, N_sat, 3]，J2000 惯性系卫星位置，单位 m
%       cfg      : 参数结构体，使用字段
%           .R_occ          遮挡用地球半径 [m]（平均半径 6371 km）
%           .occTol_m       遮挡安全余量 [m]
%           .maxIslRange_m  最大星间链路距离 [m]，inf 表示不限制
%
%   输出：
%       vis_data   : [N_slot, N_sat, N_sat] logical
%                    vis_data(k,i,j)=1 表示第 k 个时隙卫星 i 与卫星 j
%                    的连线不被地球遮挡（几何可视），否则为 0；
%                    对角线（卫星自身）恒为 0。
%       range_data : [N_slot, N_sat, N_sat] double，单位 m
%                    两星之间的“全几何直线距离”（无论是否可见均计算），
%                    对角线为 0；研究建链距离时可用 vis_data 作为掩码：
%                    linkRange = range_data .* double(vis_data)
%
%   ---------------------------------------------------------------------
%   地球遮挡判断逻辑（直线视距法 LOS，与数据集“直线视距法”一致）：
%
%       设两星位置矢量为 r_i、r_j（均为地心指向卫星），连线为线段
%       L: q(t) = r_i + t*(r_j-r_i),  t ∈ [0,1]。
%       地心 O 到该线段的最近点对应参数
%           t* = - r_i · (r_j-r_i) / |r_j-r_i|^2
%       并将 t* 截断到 [0,1]（最近点必须落在线段内部，而非延长线上）。
%       地心到星间连线的最近距离
%           d_min = | r_i + t* (r_j-r_i) |
%       判定：
%           d_min > R_earth + tol  ⇒ 连线不穿地球，可见（vis=1）
%           d_min <= R_earth + tol ⇒ 连线穿过地球，不可见（vis=0）
%
%   等高度圆轨道下的等价解析判据（便于论文中书写，代码中已用上面的
%   通用线段法实现，二者严格等价）：设两星地心夹角为 θ，则连线中点
%   即最近点，d_min = a cos(θ/2)，可见条件
%           θ < θcrit = 2 arccos(R_earth / a)
%   对 h=600 km、R=6371 km：θcrit ≈ 48.2°。
%   即同一轨道面内相邻卫星（间隔 30°）互相可见，隔一颗（60°）不可见。
%
%   全部时隙采用矩阵化运算（隐式扩展），避免逐星双层 for 循环。

[Nt, N, ~] = size(pos_data);

% 预分配：vis 用 logical（1 字节/元素），range 用 double
vis_data   = false(Nt, N, N);
range_data = zeros(Nt, N, N);

Rocc = cfg.R_occ + cfg.occTol_m;          % 遮挡判定半径（含安全余量）
maxR = cfg.maxIslRange_m;
eyeMask = logical(eye(N));                % 对角线掩码

for k = 1:Nt
    p = squeeze(pos_data(k,:,:));         % 当前时隙所有卫星位置 [N×3]

    % 配对位置矢量：A(i,1,:)=r_i，B(1,j,:)=r_j
    A = reshape(p, N, 1, 3);
    B = reshape(p, 1, N, 3);
    w = B - A;                            % [N×N×3]，w(i,j)=r_j-r_i

    % ---- 星间直线距离 |r_j-r_i| ----
    rg = sqrt(sum(w.^2, 3));              % [N×N]

    % ---- 地心到星间线段的最近距离 ----
    w2 = sum(w.^2, 3);                    % |w|^2，对角线为 0
    tStar = -sum(A .* w, 3) ./ w2;        % 最近点参数（对角线为 NaN）
    tStar = min(max(tStar, 0), 1);        % 截断到线段 [0,1] 内
    q = A + tStar .* w;                   % 线段上距地心最近的点
    dMin = sqrt(sum(q.^2, 3));            % [N×N]

    % ---- 可见性：不穿地球（且满足最大链路距离约束）----
    vis = dMin > Rocc;
    if isfinite(maxR)
        vis = vis & (rg <= maxR);
    end

    % ---- 对角线与数值对称化处理 ----
    vis(eyeMask) = false;                 % 自身不可见
    rg(eyeMask)  = 0;                     % 自身距离为 0
    vis = vis | transpose(vis);           % 消除浮点误差导致的不对称
    rg  = (rg + transpose(rg)) / 2;       % 距离矩阵强制对称

    vis_data(k,:,:)   = vis;
    range_data(k,:,:) = rg;

    if mod(k, 200) == 0
        fprintf('      已完成时隙 %d / %d\n', k, Nt);
    end
end

end
