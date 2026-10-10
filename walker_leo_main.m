%% walker_leo_main.m
% =========================================================================
%  Walker 12/12/1 低轨星座星间可见性 / 星间距离仿真
%  面向 EXata 天地一体化网络仿真的数据生成脚本（纯原生 MATLAB 实现）
% -------------------------------------------------------------------------
%  星座参数（Walker 记号 T/P/F = 144/12/1）：
%     轨道面数 P = 12，每个轨道面卫星数 S = 12，总卫星数 T = 144
%     轨道高度 h = 600 km，轨道倾角 i = 55°，相位因子 F = 1
%
%  本脚本完成的工作：
%   (1) 按 Walker-delta 规则生成 144 颗卫星的初始开普勒轨道根数；
%   (2) 用二体开普勒轨道方程解析传播，计算 24 h、步长 60 s 内每颗卫星
%       在 J2000 地心惯性系（ECI）下的三维坐标 (X,Y,Z)，单位 m；
%   (3) 直线视距法（LOS）判断两星连线是否被地球遮挡，生成可见矩阵
%         vis_data   : [时隙数, 144, 144]，logical，1=可见，0=不可见；
%   (4) 生成星间距离矩阵
%         range_data : [时隙数, 144, 144]，double，单位 m；
%   (5) 保存 .mat 文件，对齐国家基础学科公共科学数据中心（西安电子科技
%       大学）《大规模低轨卫星网络运行参数数据集》的“固定时隙 + 三维
%       位置 + 可见性 + 实时距离”结构化组织方式；
%   (6) 导出 EXata 可直接导入的节点运动轨迹：J2000 惯性系 CSV、
%       ECEF 地固系 CSV、经纬高 CSV、EXata/QualNet 块格式 mobility trace。
%
%  坐标系、遮挡逻辑、运行前提与 EXata 导入方法详见同目录
%       README_使用说明.md
%
%  依赖：仅需 MATLAB 基础环境（建议 R2016b 及以上），不需要 STK、
%        Aerospace Toolbox、Mapping Toolbox 等任何工具箱。
%
%  运行方式：在 MATLAB 中打开本文件所在目录，直接按 F5，或命令行输入
%       >> walker_leo_main
% =========================================================================

clear; clc;
tRunStart = tic;                         % 记录总运行时间

%% 0. 路径处理：把本脚本所在目录加入搜索路径，保证同目录函数文件可调用
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)                    % 兼容“新建->运行未保存文件”的情况
    scriptDir = pwd;
end
addpath(scriptDir);
outDir = fullfile(scriptDir, 'output');          % mat 输出目录
csvDir = fullfile(outDir, 'traces_csv');         % 轨迹文本输出目录
stkDir = fullfile(outDir, 'stk_ephemeris');      % STK .e 星历输出目录
if ~exist(outDir, 'dir'); mkdir(outDir); end
if ~exist(csvDir, 'dir'); mkdir(csvDir); end
if ~exist(stkDir, 'dir'); mkdir(stkDir); end

%% 1. 仿真参数配置 =======================================================
cfg = struct();

% --- 1.1 Walker 星座参数（T/P/F = 144/12/1）---
cfg.N_plane    = 12;                     % P：轨道面数
cfg.N_perPlane = 12;                     % S：每个轨道面的卫星数
cfg.N_sat      = cfg.N_plane * cfg.N_perPlane;   % T = P*S = 144 总卫星数
cfg.phaseF     = 1;                      % F：Walker 相位因子
cfg.walkerType = 'delta';                % delta 型：RAAN 在 0~360° 均匀铺开

% --- 1.2 轨道参数 ---
cfg.h_m        = 600e3;                  % 轨道高度 [m]
cfg.inc_deg    = 55;                     % 轨道倾角 [deg]
cfg.ecc        = 0;                      % 偏心率（Walker 星座取圆轨道 e=0）
cfg.argPer_deg = 0;                      % 近地点幅角 [deg]（圆轨道取 0）

% --- 1.3 地球常数（WGS84 / IERS 常用值）---
cfg.mu         = 3.986004418e14;         % 地球引力常数 GM [m^3/s^2]
cfg.Re_equ     = 6378137.0;              % WGS84 赤道半径 [m]，半长轴 a=Re+h
cfg.R_occ      = 6371000.0;              % 遮挡判断用的地球半径（平均半径）[m]
cfg.f_wgs      = 1/298.257223563;        % WGS84 扁率（ECEF->经纬高使用）

% --- 1.4 仿真时间 ---
cfg.T_total_s  = 24*3600;                % 仿真总时长：24 h = 86400 s
cfg.dt_s       = 60;                     % 时间步长：60 s
cfg.time_s     = (0:cfg.dt_s:cfg.T_total_s-cfg.dt_s).';  % 时隙时刻列向量 [s]
cfg.N_slot     = numel(cfg.time_s);      % 时隙数 = 86400/60 = 1440

% --- 1.5 仿真起始历元（UTC）---
% 星间几何（位置、可见性、距离）在惯性系中计算，与起始历元无关；
% 起始历元仅在把 J2000 惯性系坐标旋转到 ECEF 地固系（供 EXata 使用）
% 计算格林尼治恒星时（GMST）时使用。默认取 J2000.0 历元。
cfg.epochUtc   = [2000 1 1 12 0 0];      % [年 月 日 时 分 秒]，UTC

% --- 1.6 星间链路（ISL）几何约束 ---
cfg.maxIslRange_m = inf;                 % 最大建链距离 [m]，inf=不限制
cfg.occTol_m      = 1000;               % 遮挡安全余量 [m]：连线到地心最近
                                         % 距离需 > R_occ + tol，避免临界抖动

% --- 1.7 导出开关 ---
cfg.saveMat           = true;           % 保存汇总 .mat
cfg.saveSeparateMat   = true;           % 额外按“位置/可见/距离”拆分保存 3 个 mat
cfg.exportCsvEci      = true;           % 导出 J2000 惯性系位置长表 CSV
cfg.exportCsvEcef     = true;           % 导出 ECEF 地固系位置长表 CSV（EXata 推荐）
cfg.exportLlh         = true;           % 导出经度/纬度/高度 CSV
cfg.exportExataTrace  = true;           % 导出 EXata/QualNet 块格式 mobility trace
cfg.exportPerSatTrace = true;           % 同时为每颗卫星导出单独的 trace 文件
cfg.exportStkE        = true;           % 导出 STK 外部星历文件（.e，含位置与速度）

%% 2. 构建 Walker 星座初始轨道根数 =======================================
fprintf('============================================================\n');
fprintf(' Walker 12/12/1 LEO 星座仿真（h=600km, i=55°, 144 颗卫星）\n');
fprintf('============================================================\n');
orb = walker_build_orbits(cfg);
fprintf('半长轴 a = %.3f km，平均运动 n = %.8f rad/s\n', ...
        orb.a/1e3, sqrt(cfg.mu/orb.a^3));
fprintf('轨道周期 T = %.2f s = %.2f min\n', ...
        2*pi/sqrt(cfg.mu/orb.a^3), 2*pi/sqrt(cfg.mu/orb.a^3)/60);

%% 3. 开普勒轨道传播：J2000 ECI 三维位置时序 ============================
fprintf('[1/4] 开普勒轨道传播，计算 J2000 ECI 位置...\n');
[pos_data, orbInfo, vel_data] = walker_propagate_eci(orb, cfg);
% pos_data : [N_slot, N_sat, 3]，第 3 维依次为 X,Y,Z，单位 m，J2000 惯性系
% vel_data : [N_slot, N_sat, 3]，J2000 惯性系速度 Vx,Vy,Vz，单位 m/s

%% 4. 星间可见性（地球遮挡）与星间距离 ==================================
fprintf('[2/4] 星间可见性 / 距离计算，共 %d 个时隙...\n', cfg.N_slot);
[vis_data, range_data] = isl_visibility_range(pos_data, cfg);
% vis_data  : [N_slot, N_sat, N_sat] logical
% range_data: [N_slot, N_sat, N_sat] double，全几何距离（单位 m）

%% 5. 保存 .mat 数据文件 ================================================
fprintf('[3/4] 保存 mat 数据文件...\n');
sat_id   = int32(1:cfg.N_sat).';              % 卫星全局编号 1..144
plane_id = int32(orb.planeIdx);               % 所属轨道面编号 1..12
slot_s   = int32((1:cfg.N_slot).');           % 时隙编号 1..1440
time_s   = cfg.time_s;                        % 各时隙对应时刻 [s]（列向量）

if cfg.saveMat
    % 5.1 汇总 mat：一个文件包含全部数据
    matAll = fullfile(outDir, 'walker12121_leo_dataset.mat');
    save(matAll, '-v7.3', ...
        'pos_data', 'vel_data', 'vis_data', 'range_data', ...
        'time_s', 'slot_s', 'sat_id', 'plane_id', 'cfg');
    fprintf('      汇总文件: %s\n', matAll);

    % 5.2 按数据类型拆分保存（对齐数据集多文件组织，便于按需加载）
    if cfg.saveSeparateMat
        save(fullfile(outDir, 'position_data.mat'),   '-v7.3', ...
             'pos_data', 'vel_data', 'time_s', 'slot_s', 'sat_id', 'plane_id', 'cfg');
        save(fullfile(outDir, 'visibility_data.mat'), '-v7.3', ...
             'vis_data', 'time_s', 'slot_s', 'sat_id', 'cfg');
        save(fullfile(outDir, 'range_data.mat'),      '-v7.3', ...
             'range_data', 'time_s', 'slot_s', 'sat_id', 'cfg');
        fprintf('      拆分文件: position_data.mat / visibility_data.mat / range_data.mat\n');
    end
end

%% 6. 导出 EXata 节点运动轨迹文本 =======================================
fprintf('[4/4] 导出 EXata 轨迹 CSV / mobility trace');
if cfg.exportStkE; fprintf(' 与 STK .e 星历'); end
fprintf('...\n');
export_exata_csv(pos_data, cfg, csvDir);
if cfg.exportStkE
    export_stk_e(pos_data, vel_data, cfg, stkDir);
end

%% 7. 内置自检（数值正确性快速核验）=====================================
fprintf('------------------------------------------------------------\n');
fprintf('内置自检：\n');

% 7.1 维度检查
assert(isequal(size(pos_data),  [cfg.N_slot cfg.N_sat 3]), 'pos_data 维度错误');
assert(isequal(size(vis_data),  [cfg.N_slot cfg.N_sat cfg.N_sat]), 'vis_data 维度错误');
assert(isequal(size(range_data),[cfg.N_slot cfg.N_sat cfg.N_sat]), 'range_data 维度错误');
fprintf('  [PASS] 矩阵维度：pos[%d,%d,3], vis/range[%d,%d,%d]\n', ...
        cfg.N_slot, cfg.N_sat, cfg.N_slot, cfg.N_sat, cfg.N_sat);

% 7.2 卫星到地心距离应恒等于半长轴（圆轨道）
rNorm = squeeze(sqrt(sum(pos_data.^2, 3)));          % [N_slot, N_sat]
errR  = max(abs(rNorm(:) - orb.a));
assert(errR < 1e-6, '卫星地心距离不等于半长轴');
fprintf('  [PASS] 地心距离恒为 a=%.3f km（最大偏差 %.2e m）\n', orb.a/1e3, errR);

% 7.2b 速度模长应恒等于圆轨道速度 sqrt(mu/a)
vNorm = squeeze(sqrt(sum(vel_data.^2, 3)));
vTheory = sqrt(cfg.mu/orb.a);
errV = max(abs(vNorm(:) - vTheory));
assert(errV < 1e-6, '卫星速度模长异常');
fprintf('  [PASS] 速度模长恒为 %.3f m/s（最大偏差 %.2e m/s）\n', vTheory, errV);

% 7.3 vis / range 对称性与对角线
v0 = squeeze(vis_data(1,:,:));
g0 = squeeze(range_data(1,:,:));
assert(isequal(v0, v0') && norm(g0-g0','fro') < 1e-6, '矩阵不对称');
assert(~any(diag(v0)) && max(abs(diag(g0)))==0, '对角线应为不可见/零距离');
fprintf('  [PASS] vis/range 对称，对角线 vis=0、range=0\n');

% 7.4 临界可见地心角理论值：theta_crit = 2*arccos(R_occ/a)
thetaCrit = 2*acosd(cfg.R_occ/orb.a);
fprintf('  [INFO] 星间几何可见临界地心角 = %.2f°（连线与地球相切）\n', thetaCrit);

% 7.5 同轨道面相邻卫星：面内间隔 30°，应可见；隔一颗间隔 60°，应不可见
%     相邻星理论距离 = 2 a sin(15°)
dAdjTheory = 2*orb.a*sind(15);
p0 = squeeze(pos_data(1,:,:));
dAdj  = norm(p0(1,:)-p0(2,:));
dSkip = norm(p0(1,:)-p0(3,:));
assert(abs(dAdj-dAdjTheory) < 1, '同轨相邻星距离与理论值不符');
assert(v0(1,2)==true  && v0(1,3)==false, '同轨可见性与遮挡理论不符');
fprintf('  [PASS] 同轨相邻星(30°)可见，距离 %.1f km；隔一颗(60°)不可见\n', dAdj/1e3);

% 7.6 可见邻居数统计（首个时隙）
nNeigh = sum(v0, 2);
fprintf('  [INFO] 首时隙每星可见邻居数：最小 %d，最大 %d，平均 %.1f\n', ...
        min(nNeigh), max(nNeigh), mean(nNeigh));

fprintf('------------------------------------------------------------\n');
fprintf('全部完成，总耗时 %.1f s。输出目录：%s\n', toc(tRunStart), outDir);
clear v0 g0 p0 rNorm nNeigh vNorm                       % 仅释放自检临时变量
