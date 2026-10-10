function [pos_data, info, vel_data] = walker_propagate_eci(orb, cfg)
%WALKER_PROPAGATE_ECI  二体开普勒轨道传播，输出 J2000 惯性系三维位置/速度
%
%   [pos_data, info] = WALKER_PROPAGATE_ECI(orb, cfg)
%
%   物理模型：
%       二体问题（地球视为均匀球体、仅受中心引力）下，开普勒方程
%           M(t) = M0 + n*t = E - e*sin(E)
%       其中 n = sqrt(mu/a^3) 为平均运动角速度；解出偏近点角 E 后，
%       卫星在“近焦点坐标系”（轨道面内，x 轴指向近地点）中的位置为
%           x_p = a*(cos(E)-e),   y_p = a*sqrt(1-e^2)*sin(E)
%       再经经典 3-1-3 欧拉旋转（R3(-Omega)*R1(-i)*R3(-omega)）得到
%       J2000 地心惯性系（ECI）坐标。
%
%   圆轨道简化（本星座 e=0）：
%       e=0 时 E=M，纬度幅角 u = omega + nu 随时间线性变化：
%           u(t) = u0 + n*t
%       3-1-3 旋转可直接展开为（omega=0）：
%           X = a*( cos(Omega)*cos(u) - sin(Omega)*cos(i)*sin(u) )
%           Y = a*( sin(Omega)*cos(u) + cos(Omega)*cos(i)*sin(u) )
%           Z = a*( sin(i)*sin(u) )
%       上式即本函数主分支，对“时隙×卫星”全矩阵向量化计算，速度快。
%       椭圆轨道通用分支（解开普勒方程）保留在本文件末尾，供扩展使用。
%
%   坐标系说明（J2000 ECI）：
%       原点：地心；XY 平面：J2000.0 平赤道面；X 轴：指向 J2000.0
%       平春分点；Z 轴：指向平北极。坐标轴相对恒星惯性固定，地球在
%       该系中自转。星间相对几何（距离、遮挡、可见性）在惯性系中
%       计算最为直接，且与在 ECEF 中计算结果严格等价（ECEF 只是 ECI
%       绕 Z 轴整体旋转 GMST，刚体旋转不改变星间内积与遮挡关系）。
%
%   输出：
%       pos_data : [N_slot, N_sat, 3]，J2000 坐标 (X,Y,Z)，单位 m
%       info     : 附带平均运动 n、轨道周期等信息的结构体
%       vel_data : [N_slot, N_sat, 3]，J2000 速度 (Vx,Vy,Vz)，单位 m/s
%                  （圆轨道由位置对纬度幅角解析求导得到，供 STK .e 星历）
%
%   注：本函数为理想二体模型，不含 J2 摄动（轨道面在惯性系中固定）。
%       24 h 尺度上 J2 主要造成 RAAN 漂移与近地点旋转；若后续需要
%       SGP4 级精度，可在本函数中加入 J2 平均根数修正，接口保持不变。

N  = cfg.N_sat;                 % 卫星数
Nt = cfg.N_slot;                % 时隙数
t  = cfg.time_s;                % 时刻列向量 [Nt×1]，单位 s

a   = orb.a;                    % 半长轴 [m]
e   = orb.e;                    % 偏心率
inc = orb.inc;                  % 倾角 [rad]

% 平均运动与轨道周期
n_mean = sqrt(cfg.mu / a^3);    % [rad/s]
info.n_mean  = n_mean;
info.T_orbit = 2*pi / n_mean;   % [s]

% 每颗卫星的 RAAN、初始纬度幅角（转成行向量以便向量化扩展）
Omega = orb.Omega(:).';         % [1×N]
u0    = orb.u0(:).';            % [1×N]

if e == 0
    % ============ 圆轨道解析分支（Walker 星座使用）============
    % 纬度幅角随时间线性推进：u(t) = u0 + n*t
    % t 为 [Nt×1]、u0 为 [1×N]，借助隐式扩展得到 [Nt×N]
    u = mod(u0 + n_mean .* t, 2*pi);

    cosu = cos(u);  sinu = sin(u);
    cosO = cos(Omega);  sinO = sin(Omega);      % [1×N]
    cosi = cos(inc);  sini = sin(inc);          % 标量

    % J2000 ECI 坐标（3-1-3 旋转在圆轨道、omega=0 下的展开式）
    X = a * ( cosu .* cosO - sinu .* (sinO * cosi) );
    Y = a * ( cosu .* sinO + sinu .* (cosO * cosi) );
    Z = a * ( sinu * sini );

    % J2000 ECI 速度：对 u 求导后乘 du/dt = n（供 STK .e 星历使用）
    VX = n_mean*a * ( -cosO .* sinu - sinO*cosi .* cosu );
    VY = n_mean*a * ( -sinO .* sinu + cosO*cosi .* cosu );
    VZ = n_mean*a * ( sini * cosu );
else
    % ============ 椭圆轨道通用分支（保留，供扩展）============
    % 若设置非零偏心率，需要初始平近点角 M0；圆轨道配置下不会进入
    M0 = mod(u0 - orb.omega, 2*pi);                  % 近似：以纬度幅角反推
    M  = mod(M0 + n_mean .* t, 2*pi);                % [Nt×N]
    E  = kepler_solve(M, e);                         % 解 Kepler 方程
    nu = atan2(sqrt(1-e^2).*sin(E), cos(E)-e);       % 真近点角
    u  = mod(orb.omega + nu, 2*pi);

    xp = a*(cos(E) - e);                             % 近焦点系位置
    yp = a*sqrt(1-e^2).*sin(E);

    cosO = cos(Omega); sinO = sin(Omega);
    cosi = cos(inc);  sini = sin(inc);
    cosw = cos(orb.omega); sinw = sin(orb.omega);
    % 先在轨道面内旋转 -omega（近地点->升交点方向），再转倾角、RAAN
    x1 = cosw.*xp - sinw.*yp;
    y1 = sinw.*xp + cosw.*yp;
    X = cosO.*x1 - sinO.*cosi.*y1;
    Y = sinO.*x1 + cosO.*cosi.*y1;
    Z = sini.*y1;

    % 近焦点系速度（dE/dt = n/(1-e cosE)），再做相同的 3-1-3 旋转
    vxp = -a*n_mean*sinE ./ (1 - e*cosE);
    vyp =  a*sqrt(1-e^2)*n_mean*cosE ./ (1 - e*cosE);
    vx1 = cosw.*vxp - sinw.*vyp;
    vy1 = sinw.*vxp + cosw.*vyp;
    VX = cosO.*vx1 - sinO.*cosi.*vy1;
    VY = sinO.*vx1 + cosO.*cosi.*vy1;
    VZ = sini.*vy1;
end

% 组装为 [N_slot, N_sat, 3]
pos_data = zeros(Nt, N, 3);
pos_data(:,:,1) = X;
pos_data(:,:,2) = Y;
pos_data(:,:,3) = Z;

vel_data = zeros(Nt, N, 3);       % J2000 ECI 速度 [m/s]
vel_data(:,:,1) = VX;
vel_data(:,:,2) = VY;
vel_data(:,:,3) = VZ;

end


% =========================================================================
function E = kepler_solve(M, e)
%KEPLER_SOLVE  牛顿迭代法求解开普勒方程 M = E - e*sin(E)
%   M : 平近点角 [rad]（可为矩阵）
%   e : 偏心率（0 <= e < 1）
%   E : 偏近点角 [rad]
%
%   迭代格式：E_{k+1} = E_k - (E_k - e*sin(E_k) - M)/(1 - e*cos(E_k))
%   圆轨道 e=0 时直接 E=M；椭圆轨道一般 3~5 次迭代即收敛到 1e-12。

if e == 0
    E = M;
    return;
end

E = M;                         % 初值取 M（e 较小时收敛很快）
for k = 1:20
    f  = E - e*sin(E) - M;
    fp = 1 - e*cos(E);
    dE = f ./ fp;
    E  = E - dE;
    if max(abs(dE(:))) < 1e-12
        break;
    end
end
end
