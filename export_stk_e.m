function export_stk_e(pos_data, vel_data, cfg, outDir)
%EXPORT_STK_E  导出 STK 外部星历文件（.e，External Ephemeris）
%
%   EXPORT_STK_E(pos_data, vel_data, cfg, outDir) 为每颗卫星生成一个
%   STK 12 可直接加载的 .e 星历文件（sat_001.e … sat_144.e）。
%
%   ---------------------------------------------------------------------
%   生成的头部严格遵循 STK 外部星历（Ephemeris File Format）语法：
%
%       stk.v.12.0
%       BEGIN Ephemeris
%       NumberOfEphemerisPoints 1440
%       ScenarioEpoch           1 Jan 2000 12:00:00.000000
%       InterpolationMethod     Lagrange
%       InterpolationOrder      7
%       CentralBody             Earth
%       CoordinateSystem        J2000
%       DistanceUnit            Kilometers
%
%       EphemerisTimePosVel
%        0.000000000000000e+00   6978.137000   0.000000   0.000000   0.000000   7.557865  -0.000000
%        ...
%       END Ephemeris
%
%   ---------------------------------------------------------------------
%   关键点（此前版本失败的根因，务必不要再改回去）：
%
%   1) 【最重要】数据行的时间列是“相对 ScenarioEpoch 的秒数偏移”，
%      不是绝对日期！STK 外部星历的数据块格式为“时间 + 位置(+速度)”，
%      其中“时间”一律解释为自 ScenarioEpoch 起的流逝秒数（TimeFormat
%      EpSec 语义）。写成 "1 Jan 2000 12:00:00.000" 这类绝对日期串时，
%      STK 无法把它们当作合法时间读入 → 文件头能读出第一个点时刻、
%      但星历实际为空，于是 3D/2D 不显示卫星、Classical Orbit Elements
%      报告报 "The reference object has no ephemeris data"、对象起止时间
%      退化为零长度。
%
%   2) 必须有 ScenarioEpoch 行，且格式为 "d mmm yyyy hh:mm:ss.ssssss"
%      （日不补零，英文三字母月份）。
%
%   3) 插值用 STK 的关键字 InterpolationSamplesM1，不是
%      InterpolationOrder / InterpolationType。n 阶 Lagrange 插值对应
%      写 InterpolationSamplesM1 n（即 n 阶 = 每侧 n/2 个采样点）。
%
%   4) 单位：STK 外部星历默认米（Meters）；本函数写千米，因此必须显式
%      声明 DistanceUnit Kilometers，并把 m→km、m/s→km/s 换算后再写。
%
%   5) 位置精度：km 下保留 9 位小数（≈1e-6 m，远优于 60 s 步长的插值误差）；
%      速度保留 12 位小数。切勿降到 3~6 位小数，否则 60 s 采样点之间的
%      Hermite/Lagrange 插值会出现可见抖动。
%
%   参考：STK 12 Help → "Ephemeris File Format"；AGI STK Components
%         StkEphemerisFile.EphemerisTimePosVel；Orekit STKEphemerisFileParser
%         实现的 STK 关键字集合（NumberOfEphemerisPoints / ScenarioEpoch /
%         InterpolationMethod / InterpolationSamplesM1 / CentralBody /
%         CoordinateSystem / DistanceUnit / EphemerisTimePos[Vel[Acc]]）。
%   ---------------------------------------------------------------------

[Nt, N, ~] = size(pos_data);
if ~exist(outDir, 'dir'); mkdir(outDir); end

% ---- 1. 时间列：相对 ScenarioEpoch 的秒偏移 ----------------------------
if isfield(cfg, 'time_s')
    tOff = double(cfg.time_s(:));                 % 列向量 [Nt×1]，单位 s
else
    tOff = (0:Nt-1).' * cfg.dt_s;
end
assert(numel(tOff) == Nt, 'cfg.time_s 长度(%d)与星历点数(%d)不一致。', numel(tOff), Nt);

% ---- 2. 历元字符串（英文月份，日不补零，与 STK 写法一致）--------------
ek = epoch_vector(cfg.epochUtc);                  % [y m d H M S]
mon = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
epochStr = sprintf('%d %s %d %02d:%02d:%09.6f', ...
                   ek(3), mon{ek(2)}, ek(1), ek(4), ek(5), ek(6));

% ---- 3. 预生成时间列字符串（只做一次，避免逐行格式化）----------------
tLine = cell(Nt, 1);
for k = 1:Nt
    tLine{k} = sprintf('%.15e', tOff(k));         % 例如 '6.000000000000000e+01'
end

useVel = ~isempty(vel_data) && size(vel_data, 1) == Nt && size(vel_data, 2) == N;
if useVel
    dataTag = 'EphemerisTimePosVel';
else
    dataTag = 'EphemerisTimePos';                 % 仅位置（Lagrange 插值即可）
end

% ---- 4. 逐星写文件 ----------------------------------------------------
for id = 1:N
    fp  = fullfile(outDir, sprintf('sat_%03d.e', id));
    fid = fopen(fp, 'wt');
    assert(fid > 0, '无法写入文件：%s', fp);

    fprintf(fid, 'stk.v.12.0\n');
    fprintf(fid, 'BEGIN Ephemeris\n\n');
    fprintf(fid, 'NumberOfEphemerisPoints %d\n', Nt);
    fprintf(fid, 'ScenarioEpoch           %s\n', epochStr);
    fprintf(fid, 'InterpolationMethod     Lagrange\n');
    fprintf(fid, 'InterpolationSamplesM1  7\n');
    fprintf(fid, 'CentralBody             Earth\n');
    fprintf(fid, 'CoordinateSystem        J2000\n');
    fprintf(fid, 'DistanceUnit            Kilometers\n\n');
    fprintf(fid, '%s\n', dataTag);

    if useVel
        % m -> km，m/s -> km/s
        P = pos_data(:, id, :) / 1000;
        V = vel_data(:, id, :) / 1000;
        P = reshape(P, Nt, 3);
        V = reshape(V, Nt, 3);
        for k = 1:Nt
            fprintf(fid, '%s %16.9f %16.9f %16.9f %16.12f %16.12f %16.12f\n', ...
                    tLine{k}, P(k,1), P(k,2), P(k,3), V(k,1), V(k,2), V(k,3));
        end
    else
        P = reshape(pos_data(:, id, :) / 1000, Nt, 3);
        for k = 1:Nt
            fprintf(fid, '%s %16.9f %16.9f %16.9f\n', ...
                    tLine{k}, P(k,1), P(k,2), P(k,3));
        end
    end

    fprintf(fid, '\nEND Ephemeris\n');
    fclose(fid);
end

fprintf('      %s%ssat_001.e … sat_%03d.e（共 %d 个 STK .e 星历）\n', ...
        outDir, filesep, N, N);
fprintf('        格式：stk.v.12.0 / %s / DistanceUnit Kilometers / J2000\n', dataTag);
fprintf('        历元：%s（数据行时间 = 相对该历元的秒偏移，0 … %.6f s）\n', ...
        epochStr, tOff(end));

end

% ------------------------------------------------------------------------
function ek = epoch_vector(e)
%EPOCH_VECTOR  把 cfg.epochUtc 规范化为 [年 月 日 时 分 秒]
if ischar(e) || isstring(e)
    d  = datetime(char(e), 'InputFormat', 'dd MMM yyyy HH:mm:ss');
    ek = [year(d), month(d), day(d), hour(d), minute(d), second(d)];
elseif isstruct(e)
    ek = [e.year, e.month, e.day, e.hour, e.minute, e.second];
else
    e  = double(e(:)).';
    assert(numel(e) >= 6, 'cfg.epochUtc 需要 [年 月 日 时 分 秒] 六个元素。');
    ek = e(1:6);
end
end
