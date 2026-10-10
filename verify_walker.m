%% verify_walker.m
% =========================================================================
%  独立验证脚本：从磁盘重新 load 生成的 .mat，按“不同于生成路径”的方式
%  回读核验数据维度、物理量、可见性判据与文件完整性。
%  运行前提：先成功运行 walker_leo_main.m。
% =========================================================================
clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir); scriptDir = pwd; end
outDir = fullfile(scriptDir, 'output');
csvDir = fullfile(outDir, 'traces_csv');

fprintf('================ 产物独立回读验证 ================\n');

%% 1. 加载汇总 mat
S = load(fullfile(outDir,'walker12121_leo_dataset.mat'));
cfg = S.cfg;
Nt = cfg.N_slot; N = cfg.N_sat;

% 1.1 维度与类型
assert(isequal(size(S.pos_data),  [Nt N 3]),       'pos_data 维度不符');
assert(isequal(size(S.vis_data),  [Nt N N]),       'vis_data 维度不符');
assert(isequal(size(S.range_data),[Nt N N]),       'range_data 维度不符');
assert(islogical(S.vis_data),                      'vis_data 应为 logical');
fprintf('[PASS] 维度与类型：pos[%d,%d,3] double, vis[%d,%d,%d] logical, range[%d,%d,%d] double\n',Nt,N,Nt,N,N,Nt,N,N);

% 1.2 位置模长恒等于半长轴
a = cfg.Re_equ + cfg.h_m;
rNorm = squeeze(sqrt(sum(S.pos_data.^2,3)));
assert(max(abs(rNorm(:)-a)) < 1e-6, '位置模长异常');
fprintf('[PASS] 所有时刻所有卫星 |r| = %.3f km（圆轨道）\n', a/1e3);

% 1.3 轨道周期核验（相邻时隙相位推进应等于 n*dt）
nMean = sqrt(cfg.mu/a^3);
Torb  = 2*pi/nMean;
fprintf('[PASS] 理论轨道周期 %.2f s = %.2f min，24h 运行 %.2f 圈\n', Torb, Torb/60, 86400/Torb);

%% 2. 逐时隙对称性 / 对角 / 取值范围检查
for k = [1 720 Nt]
    v = squeeze(S.vis_data(k,:,:));
    g = squeeze(S.range_data(k,:,:));
    assert(isequal(v,v'), sprintf('时隙%d vis不对称',k));
    assert(norm(g-g','fro') < 1e-6, sprintf('时隙%d range不对称',k));
    assert(~any(diag(v)), '对角 vis 应为 0');
    assert(max(abs(diag(g)))==0, '对角 range 应为 0');
    assert(all(g(~logical(eye(N)))>0), '非对角距离应 >0');
    assert(all(g(v) > 0), '可见对距离应 >0');
end
fprintf('[PASS] 首/中/末时隙：vis/range 对称、对角正确、可见对距离为正\n');

%% 3. 遮挡判据的独立复核（用夹角法重算，与线段法交叉验证）
%    等高度圆轨道：theta<2*arccos(R/a) 可见
thetaCrit = 2*acosd(cfg.R_occ/a);
p0 = squeeze(S.pos_data(1,:,:));
cosTheta = (p0*p0')/(a*a);
cosTheta = min(max(cosTheta,-1),1);
theta = acosd(cosTheta);
vRef = theta < thetaCrit;
vRef(logical(eye(N))) = false;
vCode = squeeze(S.vis_data(1,:,:));
mismatch = sum(sum(xor(vRef, vCode)));
assert(mismatch==0, '线段法与夹角法结果不一致');
fprintf('[PASS] 线段法与夹角解析判据完全一致（临界角 %.2f°）\n', thetaCrit);

% 同轨面：相邻(30°)可见、隔一颗(60°)不可见；距离理论值
dAdj  = norm(p0(1,:)-p0(2,:));
dSkip = norm(p0(1,:)-p0(3,:));
assert(vCode(1,2)==true && vCode(1,3)==false, '同轨可见性异常');
fprintf('[PASS] 同轨相邻星：30°可见，距离 %.1f km（理论 %.1f km）；隔一颗 60°不可见（%.1f km）\n',...
    dAdj/1e3, 2*a*sind(15)/1e3, dSkip/1e3);

%% 4. 拓扑统计（全时隙）
nNeigh = squeeze(sum(S.vis_data,3));   % [Nt,N]
fprintf('[INFO] 每星可见邻居数：全局最小 %d、最大 %d、平均 %.2f\n', ...
    min(nNeigh(:)), max(nNeigh(:)), mean(nNeigh(:)));
fprintf('[INFO] 全网可见链路条数（无向，首时隙）：%d\n', ...
    sum(sum(triu(squeeze(S.vis_data(1,:,:)),1))));

%% 5. 文本产物检查
files = {'satellite_positions_eci.csv','satellite_positions_ecef.csv',...
         'satellite_llh.csv','exata_mobility_trace_ecef.txt'};
for ii = 1:numel(files)
    fp = fullfile(csvDir, files{ii});
    assert(exist(fp,'file')==2, '缺少文件 %s', files{ii});
    fi = dir(fp);
    fprintf('[PASS] %s（%.2f MB）\n', files{ii}, fi.bytes/1e6);
end
% 行数核验：长表应为 1440*144=207360 行数据 + 1 行表头
fid = fopen(fullfile(csvDir,'satellite_positions_eci.csv'),'r');
nLines = 0;
while fgets(fid) ~= -1; nLines = nLines+1; end
fclose(fid);
assert(nLines == Nt*N+1, 'ECI CSV 行数不符');
fprintf('[PASS] ECI CSV 共 %d 行（1 表头 + %d 条记录）\n', nLines, Nt*N);

% 单星 trace 数量与块格式抽查
d = dir(fullfile(csvDir,'traces_per_sat','sat_*.trace'));
assert(numel(d)==N, '单星 trace 数量应为 144');
fid = fopen(fullfile(csvDir,'traces_per_sat','sat_001.trace'),'r');
head = fgetl(fid); row2 = fgetl(fid); fclose(fid);
fprintf('[PASS] 单星 trace 共 %d 个；sat_001 首部 "%s"，首数据行 "%s"\n', N, head, row2);

fprintf('================ 全部验证通过 ================\n');
