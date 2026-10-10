function export_exata_csv(pos_data, cfg, csvDir)
%EXPORT_EXATA_CSV  导出 EXata 可直接导入的节点运动轨迹文本
%
%   EXPORT_EXATA_CSV(pos_data, cfg, csvDir) 将 J2000 惯性系位置时序
%   pos_data [N_slot,N_sat,3] 导出为以下文本（单位：坐标 m，时间 s）：
%
%   (1) satellite_positions_eci.csv
%         长表，表头 time_s,sat_id,X_m,Y_m,Z_m（J2000 惯性系）
%         —— 即题目要求的“每个时刻三维位置 CSV”；
%   (2) satellite_positions_ecef.csv
%         长表，time_s,sat_id,X_m,Y_m,Z_m（WGS84 地固系，EXata 推荐）；
%   (3) satellite_llh.csv
%         长表，time_s,sat_id,lat_deg,lon_deg,alt_m（经纬高，便于核查）；
%   (4) exata_mobility_trace_ecef.txt
%         EXata/QualNet 外部移动轨迹（MOBILITY-TRACE）块格式，ECEF 坐标：
%             <节点编号> <航点数>
%             <时间s> <x m> <y m> <z m>
%             ...
%         每个节点一个数据块，节点编号 1..144 与卫星全局编号一致；
%   (5) traces_per_sat/sat_001.trace ... sat_144.trace
%         每颗卫星单独一个同格式 trace 文件（部分 EXata 版本按节点
%         分别指定轨迹文件时使用），由 cfg.exportPerSatTrace 控制。
%
%   说明：EXata 场景坐标为随地球固定的地固系，故 mobility trace
%   统一使用 ECEF 坐标；ECI CSV 用于与 .mat 数据集对齐和通用分析。

[Nt, N, ~] = size(pos_data);
t = cfg.time_s;

% ---- 逐时隙把 ECI 旋转为 ECEF，并换算经纬高（预分配）----
ecef_data = zeros(Nt, N, 3);
llh_data  = zeros(Nt, N, 3);
for k = 1:Nt
    eci = squeeze(pos_data(k,:,:));                  % [N×3] J2000
    ecef = eci2ecef(eci, cfg.epochUtc, t(k));       % [N×3] ECEF
    ecef_data(k,:,:) = ecef;
    llh_data(k,:,:)  = ecef2llh(ecef, cfg.Re_equ, cfg.f_wgs);
end

% ---------------------------------------------------------------------
% (1) J2000 惯性系位置长表 CSV
% ---------------------------------------------------------------------
if cfg.exportCsvEci
    fp = fullfile(csvDir, 'satellite_positions_eci.csv');
    fid = fopen(fp, 'w');
    assert(fid > 0, '无法写入文件：%s', fp);
    fprintf(fid, 'time_s,sat_id,X_m,Y_m,Z_m\n');
    for k = 1:Nt
        blk = [repmat(t(k), N, 1), (1:N).', squeeze(pos_data(k,:,:))];
        fprintf(fid, '%.0f,%.0f,%.3f,%.3f,%.3f\n', blk.');
    end
    fclose(fid);
    fprintf('      %s\n', fp);
end

% ---------------------------------------------------------------------
% (2) ECEF 地固系位置长表 CSV
% ---------------------------------------------------------------------
if cfg.exportCsvEcef
    fp = fullfile(csvDir, 'satellite_positions_ecef.csv');
    fid = fopen(fp, 'w');
    assert(fid > 0, '无法写入文件：%s', fp);
    fprintf(fid, 'time_s,sat_id,X_ecef_m,Y_ecef_m,Z_ecef_m\n');
    for k = 1:Nt
        blk = [repmat(t(k), N, 1), (1:N).', squeeze(ecef_data(k,:,:))];
        fprintf(fid, '%.0f,%.0f,%.3f,%.3f,%.3f\n', blk.');
    end
    fclose(fid);
    fprintf('      %s\n', fp);
end

% ---------------------------------------------------------------------
% (3) 经纬高 CSV
% ---------------------------------------------------------------------
if cfg.exportLlh
    fp = fullfile(csvDir, 'satellite_llh.csv');
    fid = fopen(fp, 'w');
    assert(fid > 0, '无法写入文件：%s', fp);
    fprintf(fid, 'time_s,sat_id,lat_deg,lon_deg,alt_m\n');
    for k = 1:Nt
        blk = [repmat(t(k), N, 1), (1:N).', squeeze(llh_data(k,:,:))];
        fprintf(fid, '%.0f,%.0f,%.6f,%.6f,%.3f\n', blk.');
    end
    fclose(fid);
    fprintf('      %s\n', fp);
end

% ---------------------------------------------------------------------
% (4) EXata/QualNet 统一块格式 mobility trace（ECEF）
% ---------------------------------------------------------------------
if cfg.exportExataTrace
    fp = fullfile(csvDir, 'exata_mobility_trace_ecef.txt');
    fid = fopen(fp, 'w');
    assert(fid > 0, '无法写入文件：%s', fp);
    for id = 1:N
        fprintf(fid, '%d %d\n', id, Nt);             % 节点号 + 航点数
        blk = [t, squeeze(ecef_data(:,id,:))];       % [Nt×4]
        fprintf(fid, '%.0f %.3f %.3f %.3f\n', blk.');
    end
    fclose(fid);
    fprintf('      %s\n', fp);
end

% ---------------------------------------------------------------------
% (5) 每颗卫星单独的 trace 文件（部分 EXata 版本按节点分别加载）
% ---------------------------------------------------------------------
if cfg.exportPerSatTrace
    perDir = fullfile(csvDir, 'traces_per_sat');
    if ~exist(perDir, 'dir'); mkdir(perDir); end
    for id = 1:N
        fp = fullfile(perDir, sprintf('sat_%03d.trace', id));
        fid = fopen(fp, 'w');
        assert(fid > 0, '无法写入文件：%s', fp);
        fprintf(fid, '%d %d\n', id, Nt);
        blk = [t, squeeze(ecef_data(:,id,:))];
        fprintf(fid, '%.0f %.3f %.3f %.3f\n', blk.');
        fclose(fid);
    end
    fprintf('      %s%s*.trace（共 %d 个单星 trace）\n', perDir, filesep, N);
end

end
