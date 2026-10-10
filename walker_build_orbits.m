function orb = walker_build_orbits(cfg)
%WALKER_BUILD_ORBITS  构建 Walker-delta 星座的初始开普勒轨道根数
%
%   orb = WALKER_BUILD_ORBITS(cfg) 根据 cfg 中的 Walker 参数
%   （T/P/F、高度、倾角）生成每颗卫星的经典轨道六根数（圆轨道情形）。
%
%   Walker 星座记号 T/P/F：
%       T = 卫星总数，P = 轨道面数，F = 相位因子（0 <= F < P）
%       S = T/P = 每个轨道面内的卫星数
%
%   对第 p 个轨道面（p = 0,1,...,P-1）、面内第 s 颗卫星（s = 0,...,S-1）：
%       (1) 升交点赤经 RAAN：Omega_p = 2*pi * p / P
%           —— 各轨道面升交点在赤道上均匀分布；
%       (2) 面内相位：u_s = 2*pi * s / S
%           —— 同一轨道面内卫星均匀分布；
%       (3) 相邻轨道面相位偏移：dPhi_p = 2*pi * F * p / T
%           —— 相位因子 F 决定相邻轨道面卫星之间的错位，
%              使星座覆盖与星间链路几何更均匀。
%       故第 (p,s) 颗卫星的初始纬度幅角（圆轨道下即从升交点起算的
%       在轨相位）为：u0 = u_s + dPhi_p。
%
%   卫星全局编号规则（与后续 vis_data / range_data 的行列号、
%   EXata 节点编号一一对应）：
%       satId = p*S + (s+1)，即第 1 个轨道面排满后再排第 2 个轨道面。
%
%   输出 orb 字段（角度单位均为 rad，长度单位 m）：
%       .a          半长轴（标量，所有卫星相同）
%       .e          偏心率（标量，0）
%       .inc        轨道倾角（标量，rad）
%       .omega      近地点幅角（标量，rad，圆轨道取 0）
%       .Omega[N×1] 每颗卫星的升交点赤经 RAAN
%       .u0[N×1]    每颗卫星的初始纬度幅角
%       .planeIdx[N×1]   轨道面编号（1..P）
%       .inPlaneIdx[N×1] 面内编号（1..S）
%
%   适用对象：通信方向研究生。理解 RAAN / 纬度幅角 / 相位因子是阅读
%   Walker 星座拓扑类论文的基础。

% ---------------- 基本参数 ----------------
P = cfg.N_plane;        % 轨道面数
S = cfg.N_perPlane;     % 每面卫星数
T = cfg.N_sat;          % 总卫星数 = P*S
F = cfg.phaseF;         % 相位因子

% ---------------- 轨道形状参数 ----------------
orb.a     = cfg.Re_equ + cfg.h_m;        % 半长轴 = 地球赤道半径 + 轨道高度 [m]
orb.e     = cfg.ecc;                     % 偏心率（圆轨道 = 0）
orb.inc   = deg2rad(cfg.inc_deg);        % 倾角 [rad]
orb.omega = deg2rad(cfg.argPer_deg);     % 近地点幅角 [rad]

% ---------------- 逐星生成轨道面与相位 ----------------
orb.Omega       = zeros(T, 1);
orb.u0          = zeros(T, 1);
orb.planeIdx    = zeros(T, 1, 'int32');
orb.inPlaneIdx  = zeros(T, 1, 'int32');

idx = 0;
for p = 0:P-1                  % 轨道面索引（从 0 起，便于套公式）
    for s = 0:S-1              % 面内卫星索引（从 0 起）
        idx = idx + 1;

        % 升交点赤经：各轨道面在赤道上均匀间隔 2π/P
        orb.Omega(idx) = 2*pi * p / P;

        % 初始纬度幅角 = 面内均匀相位 + 相邻面 Walker 相位偏移
        u_inplane = 2*pi * s / S;             % 面内相位
        dPhi      = 2*pi * F * p / T;         % 第 p 面相对第 0 面的相位偏移
        orb.u0(idx) = mod(u_inplane + dPhi, 2*pi);

        % 编号记录（从 1 起，便于与 MATLAB 矩阵索引、EXata 节点号对应）
        orb.planeIdx(idx)   = p + 1;
        orb.inPlaneIdx(idx) = s + 1;
    end
end

end
