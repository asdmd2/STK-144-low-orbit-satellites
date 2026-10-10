function llh = ecef2llh(ecef, Re, f)
%ECEF2LLH  WGS84 地心地固坐标 (X,Y,Z) -> 大地纬度、经度、椭球高
%
%   llh = ECEF2LLH(ecef, Re, f)
%
%   输入：
%       ecef : [M×3]，ECEF 坐标，单位 m
%       Re   : WGS84 椭球长半轴（赤道半径）[m]
%       f    : WGS84 扁率
%   输出：
%       llh  : [M×3]，依次为 纬度[deg]、经度[deg]、椭球高[m]
%
%   采用对 h 的不动点迭代（Bowring 迭代形式），纯 MATLAB 实现，
%   不依赖 Mapping Toolbox；初值 h=0，通常 4~6 次迭代收敛到 1e-6 m。
%
%   公式：
%       lon = atan2(Y, X)
%       p   = sqrt(X^2+Y^2)
%       N(lat) = a / sqrt(1 - e2*sin(lat)^2)      （卯酉圈曲率半径）
%       lat = atan2(Z, p*(1 - e2*N/(N+h)))
%       h   = p/cos(lat) - N
%       e2 = f*(2-f)

e2 = f*(2 - f);                              % 第一偏心率平方
X = ecef(:,1); Y = ecef(:,2); Z = ecef(:,3);

lon = atan2(Y, X);                           % 经度可直接解析求出
p   = sqrt(X.^2 + Y.^2);

lat = atan2(Z, p*(1 - e2));                  % 纬度初值（取 h=0）
h   = zeros(size(Z));
for it = 1:10
    N  = Re ./ sqrt(1 - e2*sin(lat).^2);
    hNew = p./max(cos(lat),eps) - N;
    latNew = atan2(Z, p.*(1 - e2.*N./(N + hNew)));
    if max(abs(hNew - h)) < 1e-6 && max(abs(latNew - lat)) < 1e-12
        lat = latNew; h = hNew;
        break;
    end
    lat = latNew; h = hNew;
end
N = Re ./ sqrt(1 - e2*sin(lat).^2);
h = p./max(cos(lat),eps) - N;                % 再刷新一次高度，保证末值一致

llh = [rad2deg(lat), rad2deg(lon), h];
end
