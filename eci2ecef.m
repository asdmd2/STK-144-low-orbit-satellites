function ecef = eci2ecef(eci, epochUtc, tSec)
%ECI2ECEF  J2000 地心惯性系（ECI）-> WGS84 地心地固系（ECEF）
%
%   ecef = ECI2ECEF(eci, epochUtc, tSec)
%
%   输入：
%       eci      : [M×3]，某一时刻所有卫星（或任意点）的 J2000 坐标，单位 m
%       epochUtc : 仿真起始历元 [年 月 日 时 分 秒]（UTC）
%       tSec     : 该时刻相对起始历元的秒数（标量）
%
%   输出：
%       ecef     : [M×3]，WGS84 地固系坐标 (X,Y,Z)，单位 m
%
%   原理：
%       ECEF 坐标系随地球以格林尼治平恒星时 GMST 相对 J2000 惯性系
%       绕 Z 轴旋转。忽略岁差/章动/极移（对网络拓扑级仿真影响在
%       角秒量级，完全可接受），有
%           ECEF = Rz(-GMST) * ECI
%       即
%           X_ecef =  cos(g)*X_eci + sin(g)*Y_eci
%           Y_ecef = -sin(g)*X_eci + cos(g)*Y_eci
%           Z_ecef =  Z_eci
%       其中 g = GMST（rad）。
%
%   为什么需要 ECEF：
%       EXata/QualNet 的场景坐标是“随地球固定”的地固笛卡尔坐标，
%       地面节点（关口站、用户终端）在该系中静止。卫星轨迹若直接用
%       ECI 导入，相对地面节点会多出一圈/天的虚假自转，因此导入
%       EXata 的轨迹必须使用 ECEF（或经纬高）。

% 当前时刻儒略日 = 起始历元儒略日 + 经过秒数/86400
% （直接在 JD 上加天数，可正确处理跨日/跨月，避免手动进位出错）
JD0 = datetime2jd(epochUtc);
JD  = JD0 + tSec/86400;
g   = gmst_rad(JD);                          % 格林尼治平恒星时 [rad]

cg = cos(g); sg = sin(g);
ecef = zeros(size(eci));
ecef(:,1) =  cg*eci(:,1) + sg*eci(:,2);
ecef(:,2) = -sg*eci(:,1) + cg*eci(:,2);
ecef(:,3) =  eci(:,3);

end


% =========================================================================
function JD = datetime2jd(utc)
%DATETIME2JD  公历日期时间（UTC）转儒略日 JD
%   utc = [年 月 日 时 分 秒]
%   采用通用 Gregorian 历法公式。校验：2000-01-01 12:00 UTC -> 2451545.0

y = utc(1); m = utc(2); d = utc(3);
H = utc(4) + utc(5)/60 + utc(6)/3600;

a = floor((14 - m)/12);
y = y + 4800 - a;
m = m + 12*a - 3;
JD = d + floor((153*m + 2)/5) + 365*y ...
     + floor(y/4) - floor(y/100) + floor(y/400) - 32045 ...
     + (H - 12)/24;
end


% =========================================================================
function g = gmst_rad(JD)
%GMST_RAD  给定儒略日（UT1≈UTC）计算格林尼治平恒星时，单位 rad
%   采用 Vallado《Fundamentals of Astrodynamics and Applications》
%   中的 IERS 常用多项式（度），再取模 360° 并转弧度。

T  = (JD - 2451545.0) / 36525;                         % 儒略世纪数
gmstDeg = 280.46061837 ...
        + 360.98564736629*(JD - 2451545.0) ...
        + 0.000387933*T.^2 ...
        - T.^3/38710000;
gmstDeg = mod(gmstDeg, 360);
g = deg2rad(gmstDeg);
end
