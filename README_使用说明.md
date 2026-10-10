# Walker 12/12/1 低轨星座仿真（EXata 天地一体化网络）使用说明

本工程用**纯原生 MATLAB**（不使用 STK、Aerospace Toolbox、Mapping Toolbox）复现
Walker 12/12/1 低轨星座，输出 24 小时、60 s 步长的卫星 J2000 惯性系位置、
星间可见性矩阵、星间距离矩阵，并导出 EXata/QualNet 可直接导入的节点运动轨迹。

- 星座：轨道面数 P=12，每面 S=12 颗，共 **144 颗**；轨道高度 **600 km**，
  倾角 **55°**，Walker 相位因子 **F=1**（Walker 记号 **T/P/F = 144/12/1**，delta 型）
- 仿真时长：24 h（86400 s），步长 60 s，共 **1440 个时隙**
- 已在 MATLAB R2026a 实测通过（语法兼容 R2016b 及以上），整机运行约 **10 s**

---

## 一、文件清单

```
Walker_LEO_EXata\
├── walker_leo_main.m         主脚本：参数配置 + 全流程编排 + 内置自检（直接运行它）
├── walker_build_orbits.m     函数：按 Walker-delta 规则生成初始开普勒轨道根数
├── walker_propagate_eci.m    函数：二体开普勒传播 -> J2000 ECI 三维坐标（含开普勒方程求解）
├── isl_visibility_range.m    函数：直线视距法地球遮挡判断 + 星间距离（向量化）
├── eci2ecef.m                函数：J2000 ECI -> WGS84 ECEF（含儒略日、GMST 计算）
├── ecef2llh.m                函数：ECEF -> 大地纬度/经度/椭球高（Bowring 迭代）
├── export_exata_csv.m        函数：导出 ECI/ECEF/经纬高 CSV 与 EXata mobility trace
├── export_stk_e.m            函数：导出 STK 外部星历文件（.e，J2000 位置）
├── verify_walker.m           独立验证脚本：从磁盘回读 mat，交叉核验全部产物
├── run_log.txt               最近一次运行的完整日志
├── README_使用说明.md         本文件
└── output\                   运行后自动生成
    ├── walker12121_leo_dataset.mat   汇总 mat（全部变量，含 vel_data 速度）
    ├── position_data.mat            仅位置/速度时序
    ├── visibility_data.mat          仅可见性矩阵
    ├── range_data.mat               仅距离矩阵
    ├── stk_ephemeris\sat_001.e … sat_144.e   STK 可直接导入的星历（144 个）
    └── traces_csv\
        ├── satellite_positions_eci.csv    J2000 惯性系位置长表（题目要求的 CSV）
        ├── satellite_positions_ecef.csv   ECEF 地固系位置长表（EXata 推荐）
        ├── satellite_llh.csv              经度/纬度/高度长表
        ├── exata_mobility_trace_ecef.txt  EXata/QualNet 统一块格式轨迹
        └── traces_per_sat\sat_001.trace … sat_144.trace   每星单独轨迹（共144个）
```

---

## 二、坐标系说明

| 坐标系 | 是否随地球自转 | 原点 / 坐标轴 | 在本工程中的用途 |
|---|---|---|---|
| **J2000 地心惯性系 ECI** | 否（相对恒星固定） | 原点地心；XY 平面为 J2000.0 平赤道面；X 轴指向 J2000.0 平春分点；Z 轴指向平北极 | `pos_data`、星间距离与遮挡判断均在此系完成 |
| **WGS84 地心地固系 ECEF** | 是（随地球自转） | 原点地心；X 轴指向赤道/本初子午线交点，Z 轴指向协议地极 | **EXata 节点轨迹统一用此系**（地面节点在该系静止） |
| **大地坐标系 LLH** | — | 纬度、经度、WGS84 椭球高 | 结果核查、支持经纬高航点的 EXata 卫星扩展 |

要点：

1. **星间几何在 ECI 中计算**。ECEF 只是 ECI 绕 Z 轴整体旋转一个格林尼治恒星时
   （GMST），刚体旋转不改变卫星之间的内积、距离与遮挡关系，因此
   `vis_data / range_data` 在 ECI 中计算与在 ECEF 中结果严格等价，但 ECI 中
   公式最简单（轨道面固定）。
2. **导入 EXata 的轨迹必须用 ECEF**。EXata/QualNet 的场景坐标是地固笛卡尔坐标，
   地面终端/关口站在该系中静止；若直接导入 ECI 轨迹，卫星相对地面会多出
   “一圈/天”的虚假自转。本工程导出的 trace 已全部转换为 ECEF。
3. ECI→ECEF 采用 `ECEF = Rz(−GMST)·ECI`，GMST 由儒略日按 IERS 常用多项式计算，
   忽略岁差/章动/极移（角秒量级，对网络拓扑仿真无影响）。起始历元默认
   J2000.0（2000-01-01 12:00 UTC），可在 `cfg.epochUtc` 修改；历元只影响
   ECEF/经纬度，不影响星间矩阵。
4. 本工程为**理想二体开普勒模型**（轨道面在惯性系固定），不含 J2 摄动；
   24 h 尺度上 J2 主要造成升交点赤经漂移，需要 SGP4 级精度时可在
   `walker_propagate_eci.m` 中加入 J2 平均根数修正，对外接口不变。

---

## 三、轨道与星座数学模型

### 3.1 Walker-delta 初始轨道根数

记号 T/P/F：T=144 总星数，P=12 轨道面数，F=1 相位因子，S=T/P=12 每面星数。
对第 p 个轨道面（p=0…P−1）、面内第 s 颗（s=0…S−1）：

```
半长轴       a   = Re + h = 6378.137 + 600 = 6978.137 km
升交点赤经   Ω_p = 2π·p / P                       （各面在赤道上均匀分布）
初始纬度幅角 u0  = 2π·s/S + 2π·F·p/T              （面内均匀 + 相邻面相位错位）
偏心率 e = 0，倾角 i = 55°，近地点幅角 ω = 0
```

卫星全局编号 `satId = p·S + (s+1)`（第 1 面排满再排第 2 面），与矩阵行列号、
EXata 节点号 1…144 一一对应。

### 3.2 开普勒传播（二体）

平均运动角速度与轨道周期：

```
n = sqrt(μ/a³) = 0.00108308 rad/s，   T_orb = 2π/n = 5801.23 s ≈ 96.69 min
```

圆轨道下纬度幅角线性推进 `u(t) = u0 + n·t`，解开普勒方程 M=E−e·sinE 在 e=0
时退化为 E=M。经 3-1-3 欧拉旋转（R3(−Ω)R1(−i)R3(−ω)）得到 J2000 坐标：

```
X = a·( cosΩ cosu − sinΩ cosi sinu )
Y = a·( sinΩ cosu + cosΩ cosi sinu )
Z = a·( sini sinu )
```

椭圆轨道通用分支（牛顿迭代解开普勒方程）已保留在
`walker_propagate_eci.m` 中，把 `cfg.ecc` 改为非零即可启用。

### 3.3 星间距离

```
range(i,j) = | r_j − r_i |  （单位 m，全几何直线距离，无论是否可见均计算）
```

---

## 四、地球遮挡（星间可见性）判断逻辑

采用与数据集一致的**直线视距法（LOS）**：判断两星连线线段是否与地球球体相交。

设两星位置矢量为 r_i、r_j，连线参数方程 `q(t) = r_i + t·(r_j−r_i)，t∈[0,1]`。
地心到该线段的最近点参数为

```
t* = − r_i·(r_j−r_i) / |r_j−r_i|² ，  并截断到 [0,1]（保证最近点在线段内部）
d_min = | r_i + t*·(r_j−r_i) |        （地心到星间连线的最近距离）
```

判定（`cfg.R_occ=6371 km` 为地球平均半径，另加 `cfg.occTol_m=1000 m` 安全余量）：

```
d_min > R_earth + tol   ⇒ 可见 vis = 1
d_min ≤ R_earth + tol   ⇒ 连线穿地，不可见 vis = 0
```

**等高度圆轨道等价解析判据**（代码用通用线段法实现，验证脚本证明两者逐对一致）：
设两星地心夹角为 θ，连线中点即最近点，`d_min = a·cos(θ/2)`，故

```
θ < θcrit = 2·arccos(R_earth/a)
```

代入 a=6978.137 km、R=6371 km：**θcrit ≈ 48.16°**。物理含义与实测核验：

- 同一轨道面相邻卫星间隔 30°（< 48.16°）→ **可见**，距离
  `2a·sin15° = 3612.1 km`；
- 同轨隔一颗间隔 60°（> 48.16°）→ **不可见**，距离 6978.1 km；
- 仿真统计：每颗星可见邻居数 20～25，平均 22.77；首时隙全网无向可见链路 1644 条。

如需限制星间激光/微波链路的最大作用距离，在主脚本中设置
`cfg.maxIslRange_m`（例如 5016e3），可见性会额外要求 `range ≤ 该值`。

---

## 五、运行前提与步骤

### 前提
- MATLAB R2016b 或更高版本（使用了结构体、隐式扩展；本机 R2026a 实测通过）；
- **无需任何工具箱**：轨道、坐标转换、遮挡判断全部为手写基础数学；
- 无需联网、无需 STK；
- 磁盘剩余空间建议 ≥ 400 MB（mat 用 double 存储，汇总 mat 约 136 MB）。

### 步骤
1. 将整个 `Walker_LEO_EXata` 文件夹复制到本地，在 MATLAB 中“切换到该文件夹”；
2. 打开 `walker_leo_main.m`，按 F5，或命令行执行 `>> walker_leo_main`；
3. 控制台会依次打印轨道周期、进度、内置自检结果（全部 `[PASS]`）与输出路径；
4. （可选）运行 `>> verify_walker`，从磁盘重新加载产物做独立交叉验证；
5. 所有产物在 `output\` 及其 `traces_csv\` 子目录下。

常用修改（均在 `walker_leo_main.m` 第 1 节参数区）：

| 需求 | 修改项 |
|---|---|
| 改仿真时长/步长 | `cfg.T_total_s`、`cfg.dt_s` |
| 改星座（如 24×66、StarLink 壳层） | `cfg.N_plane`、`cfg.N_perPlane`、`cfg.phaseF`、`cfg.h_m`、`cfg.inc_deg` |
| 限制最大建链距离 | `cfg.maxIslRange_m` |
| 更保守的遮挡（用赤道半径） | `cfg.R_occ = 6378137.0` |
| 距离矩阵改单精度节省一半空间 | 主脚本预分配处 `zeros(...)` 改为 `zeros(...,'single')` |
| 改起始历元（影响 ECEF/经纬度） | `cfg.epochUtc` |

---

## 六、mat 数据格式（与西电数据集对齐）

对齐对象：国家基础学科公共科学数据中心、西安电子科技大学发布的
《大规模低轨卫星网络运行参数数据集》
（CSTR:16666.11.nbsdc.mVw79QFM；其方法学为“直线视距法 + 开普勒轨道模型”、
固定时隙结构化存储）。本工程与其对齐的要素：**固定 60 s 时隙、1440 时隙/天、
J2000 三维位置、星间可见性 0/1 矩阵、星间实时距离矩阵、米/秒单位**。

汇总文件 `walker12121_leo_dataset.mat`（同时拆分为 position / visibility /
range 三个 mat，便于按需加载）：

| 变量名 | 类型 | 维度 | 含义与单位 |
|---|---|---|---|
| `pos_data` | double | **[1440, 144, 3]** | 各时隙卫星 J2000 坐标 (X,Y,Z)，单位 m |
| `vis_data` | logical | **[1440, 144, 144]** | 星间可见性，1=可见，0=不可见，对角线为 0 |
| `range_data` | double | **[1440, 144, 144]** | 星间直线距离，单位 m，对角线为 0 |
| `time_s` | double | [1440, 1] | 各时隙时刻 0,60,…,86340 s |
| `slot_s` | int32 | [1440, 1] | 时隙编号 1…1440 |
| `sat_id` | int32 | [144, 1] | 卫星全局编号 1…144 |
| `plane_id` | int32 | [144, 1] | 每颗卫星所属轨道面 1…12 |
| `cfg` | struct | — | 全部星座/地球/时间参数与坐标系配置，随数据自描述 |

说明：
- `range_data` 保存的是**全几何距离**（不可见的卫星对也有距离值），做路由/建链
  分析时用可见性掩码即可：`linkRange = range_data .* double(vis_data);`
- 第 1 维为时隙、第 2/3 维为卫星编号；矩阵对称，即 `vis(k,i,j)=vis(k,j,i)`。
- mat 以 `-v7.3` 保存，可被 MATLAB、Python（`scipy.io.loadmat(..., verify_decompressed_data=True)`
  对 v7.3 需用 `h5py`，或在 MATLAB 中另存 v7）读取。Python 读 v7.3 示例：

```python
import h5py, numpy as np
f = h5py.File('walker12121_leo_dataset.mat','r')
pos   = np.array(f['pos_data'])      # MATLAB v7.3 按列存储，注意维度转置
vis   = np.array(f['vis_data'])
rng   = np.array(f['range_data'])
```

MATLAB 加载与取数示例：

```matlab
S = load('walker12121_leo_dataset.mat');
k = 100;                 % 第 100 个时隙（t = 99*60 s）
p = squeeze(S.pos_data(k,:,:));      % 144×3，该时隙所有卫星 J2000 位置
V = squeeze(S.vis_data(k,:,:));      % 144×144 可见性
D = squeeze(S.range_data(k,:,:));    % 144×144 距离(m)
```

> 注：该数据集原始文件需在平台注册申请下载。若下载后发现其变量命名与本工程
> 存在细微差异（例如位置变量名不同），在 MATLAB 中按其名字重命名保存即可，
> 维度顺序、单位、0/1 语义均已对齐，无需重新计算。

---

## 七、把生成数据导入 EXata

### 7.1 该用哪个文件

| EXata 接入方式 | 使用文件 |
|---|---|
| 节点外部移动轨迹（MOBILITY-TRACE，地固笛卡尔坐标，**推荐**） | `exata_mobility_trace_ecef.txt`（全部节点一个块文件）或 `traces_per_sat\sat_001.trace …`（每节点一个文件） |
| EXata 卫星扩展包按“经纬度+高度”航点导入 | `satellite_llh.csv` |
| 自行二次开发/与 mat 拓扑联合仿真 | `satellite_positions_ecef.csv`、各 `.mat` |

**务必使用 ECEF/经纬高文件，不要用 ECI 文件**（理由见第二节）。

### 7.2 EXata/QualNet trace 文件格式

`exata_mobility_trace_ecef.txt` 内部按节点分块，每块第一行为
“节点编号 + 航点个数”，后续每行一个航点“时间(s) X(m) Y(m) Z(m)”：

```
1 1440
0 1266948.132 6862159.895 0.000
60 ......
...（共1440行）
2 1440
0 ......
...（直到节点 144）
```

坐标为 WGS84 ECEF、单位米；时间为仿真相对秒、0 对应起始历元；节点编号 1…144
与 `sat_id` 完全一致。单星文件 `sat_001.trace` 只含一个节点块。

### 7.3 在 EXata Architect 中的操作（通用流程，不同版本菜单位置略有差异）

1. **批量建立 144 个卫星节点**：在 Architect 中用脚本/批量复制（Nodes 面板设置
   数量与编号 1…144），节点类型选卫星/空中节点（Satellite / Aerial Node）；
   地面关口站、用户终端按常规地面节点另行放置（它们在 ECEF 场景中静止）。
2. **指定移动模型**：选中卫星节点 → 节点属性 → Mobility，将 Mobility Model
   设为 **Trace Based / MOBILITY-TRACE（外部轨迹）**；
   - 若版本支持一个总文件：轨迹文件选 `exata_mobility_trace_ecef.txt`，
     解析器按“节点号 + 航点数”分块读取；
   - 若版本要求每节点一个文件：节点 n 对应 `traces_per_sat\sat_00n.trace`
     （可在节点批量属性中按编号规则填充路径）。
3. **坐标与单位**：EXata 场景全局坐标为米、仿真时间为秒，与 trace 完全一致，
   无需缩放；ECEF 坐标原点为地心，节点高度量级为 6.98×10⁶ m，若视图中看不到
   节点，请把视图缩放/跟随对象设为卫星节点，并关闭地形对空中节点的遮挡。
4. **时间对齐**：令 EXata 场景仿真起始时刻对应 trace 的 t=0；步长可设为 60 s 的
   约数（如 1 s），EXata 会在航点之间线性插值，60 s 内卫星运动约 460 km，
   插值误差对网络拓扑无影响（可见性在 60 s 边界已按最严几何给出）。
5. **星间链路**：在无线/光链路模型中按 `vis_data` 配置邻居或开启基于几何
   通视的链路判定；也可直接用 mat 中的拓扑时序驱动自定义链路上下线
   （`vis_data(k,i,j)=1` 建链、0 拆链）。
6. **运行**：运行仿真，卫星节点将按 24 h 真实轨道运动，地面节点固定，
   星地/星间几何与 mat 数据一致。

> 若所用 EXata 版本的 trace 首行格式不同（例如不含节点号、或要求
> `time, x, y, z` 逗号分隔），直接使用 `satellite_positions_ecef.csv`
> 删除表头后按其格式微调即可；所有数值已按米/秒生成，不需要 STK 中转。

---

## 八、在 STK 中查看星座效果

运行主脚本后，`output\stk_ephemeris\` 下已生成 **sat_001.e … sat_144.e**
（STK External Ephemeris 外部星历：J2000 惯性系、位置 km、速度 km/s，
历元 `1 Jan 2000 12:00:00` UTCG，数据行时间为**相对历元的秒偏移**，
Lagrange 7 阶插值）。下面给两种方式，**方式 A 看本工程数据，
方式 B 用 STK 原生工具交叉核对**。

> **为什么数据行时间不是绝对日期？**
> STK 外部星历（`*.e`）的数据块里，时间列一律是“自 `ScenarioEpoch` 起的
> 流逝秒数”，即 `TimeFormat EpSec` 语义；绝对日期串（如
> `1 Jan 2000 12:00:00.000`）不是合法的时间列。早期版本把绝对日期写进
> 数据行，结果 STK 只读出了文件头（属性页能显示文件名和第一个点时刻），
> 星历本体为空 —— 3D/2D 不显示卫星、`Classical Orbit Elements` 报告报
> *"The reference object has no ephemeris data"*、对象起止时间退化成
> 零长度。现版本已改为标准写法，详见 `export_stk_e.m` 头部注释。

### 方式 A：导入 MATLAB 生成的 .e 星历（推荐，用于验证本工程）

1. **新建场景并设置时间**（星历覆盖 2000-01-01 12:00 → 2000-01-02 12:00
   UTCG，场景时段必须覆盖它）：
   - 打开 STK → New Scenario；
   - Start = `1 Jan 2000 12:00:00.000`（UTCG），
     Stop = `2 Jan 2000 12:00:00.000`（与 `cfg.epochUtc` + 24 h 对应）。
2. **批量导入 144 个星历**：
   - 菜单 `File → Import → STK Data File...`（或工具栏 Import 向导）；
   - 文件类型选 **STK External Ephemeris (*.e)**；
   - 进入 `output\stk_ephemeris\`，Ctrl+A 全选 144 个 `.e` → Open；
   - 导入向导中坐标系保持 **J2000/ICRF**（文件头已声明，STK 会自动识别），
     完成后对象浏览器中出现 144 个 Satellite（名称取自文件名）。
3. **建星座分组（可选但推荐）**：`Insert → Constellation`，把 144 颗卫星
   全部加入，可在 Constellation 的 3D Graphics 里统一颜色、轨道、点亮/熄灭，
   便于整体观察 12 个轨道面与 F=1 的相位错位。
4. **播放动画**：时间轴回到 1 Jan 2000 12:00，点播放，即可看到 24 h 星座
   运行（每颗星 24 h 绕地 14.89 圈）。
5. **与 vis_data 交叉验证星间可见性**：
   - 选中两颗卫星 → `Analysis → Access...`（或先建 Chain）；
   - Access 的约束中默认已包含**地球遮挡（Central Body Line of Sight）**；
   - Compute 后 Access 时间段就是该星对“可见”的时段，应与
     `vis_data(:,i,j)` 为 1 的时隙一一对应（临界地心角 48.16°）；
   - 批量分析可用一颗星对 Constellation 做 Access，导出 Access 区间与
     mat 矩阵逐时隙比对。
6. 机器性能一般时，可先只导入 12 个 `.e`（一个轨道面）确认效果，再全量导入。

7. **加载成功后 self-check（30 秒确认星历真的生成了）**：

   | 检查项 | 正确表现 | 星历为空的错误表现 |
   |---|---|---|
   | StkExternal 属性页 Time of first ephemeris point | `1 Jan 2000 12:00:00.000000 UTCG` | 同样显示（**不能**用它判断成功） |
   | 3D 窗口 | 出现卫星点 + 轨道线 | 什么都没有，Zoom To 无效 |
   | 对象 `Properties → Basic → Orbit` 的 Start/Stop | 12:00:00 → 次日 11:59:00 | 起止都等于 12:00:00（零长度） |
   | `Classical Orbit Elements` 报告 | 正常出表：a≈6978.137 km、e≈0、i≈55°、T≈96.69 min | 报 *The reference object has no ephemeris data* |

   只要第 3、4 行是右列的“错误表现”，就说明 STK 没有读到数据行——**先检查
   `.e` 文件是不是旧版本**（数据行时间列是绝对日期就一定是旧的），重新运行
   `walker_leo_main.m` 生成新文件即可。

   也可以完全不打开 STK 先自检：`output\stk_ephemeris\sat_001.e` 的前 12 行
   应当长这样（注意 `ScenarioEpoch` 行存在、`InterpolationSamplesM1 7`、
   数据行以秒偏移 `0` 开头）：

   ```
   stk.v.12.0
   BEGIN Ephemeris

   NumberOfEphemerisPoints 1440
   ScenarioEpoch           1 Jan 2000 12:00:00.000000
   InterpolationMethod     Lagrange
   InterpolationSamplesM1  7
   CentralBody             Earth
   CoordinateSystem        J2000
   DistanceUnit            Kilometers

   EphemerisTimePosVel
   0  6978.137000000  0.000000000  0.000000000  0  7.557865000000  0
   60  6963.407824000  259.917774000  371.201051000  ...
   ```

### 方式 B：STK 原生 Walker 工具（约 30 秒，核对几何布局）

1. 新建 Scenario（时间任意，二体星座的惯性几何与日期无关）。
2. 建一颗**种子卫星**：`Insert → Satellite`，Propagator 选 **Two-body**，设置
   - 半长轴/轨道高度：Perigee = Apogee = **600 km**（a = 6978.137 km）；
   - Inclination = **55°**，RAAN = 0°，Arg of Perigee = 0°，True/Mean Anomaly = 0°。
3. 选中该种子星 → 菜单 `Satellite → Walker...`（部分版本在
   `Utilities → Walker Tool`）：
   - Type = **Delta**；RAAN Spread = **360°**；
   - Number of Planes = **12**；Satellites per Plane = **12**；
   - Inter-plane Spacing（相位因子）= **1**；
   - Apply → 自动展开为 144 颗。
4. 同样建议建 Constellation 分组观察。方式 B 与方式 A 的星座**整体几何完全
   一致**；若个别卫星编号对不上，只是 Walker 展开顺序/编号约定不同。

### 常见问题

| 现象 | 原因与处理 |
|---|---|
| 属性页能读出文件名和第一个点时刻，但 3D/2D 没有卫星；报告提示 *The reference object has no ephemeris data*；对象起止时间都是 12:00:00 | **数据行时间列用了绝对日期**（旧版本 bug）。必须是从 `ScenarioEpoch` 起的**秒偏移**。重新运行 `walker_leo_main.m`，或运行 `tools\fix_stk_e.ps1` 就地改写现有 `.e` 文件 |
| 导入后提示时间不在场景范围内、卫星不显示 | 场景 Start/Stop 未设为 2000-01-01 12:00 ~ 2000-01-02 12:00（UTCG） |
| 导入向导询问坐标系 | 选 J2000（或 ICRF）；文件头已写 `CoordinateSystem J2000` |
| 卫星显示但轨迹有轻微抖动 | 星历点数与步长不匹配，或位置列精度过低；本工程位置保留 9 位小数（km）、60 s 步长，正常无抖动 |
| 144 个对象动画卡顿 | 3D Globe 关闭云层/高分辨率纹理，或先用一个轨道面 12 颗调试 |
| 看不到星间连线 | 对 Access/Chain 对象在其 `3D Graphics` 中打开连线显示 |
| STK 版本兼容性 | `.e` 文件头写 `stk.v.12.0`（STK 12 原生版本号）；STK 11 及更早版本请把首行改成 `stk.v.11.0`，其余内容不变 |

### 想让 STK 直接读取生成的文件（可选）

`.e` 文件本身与路径无关，放在任何目录都能加载。若希望每次重新生成后
STK 里指向的仍是同一批文件，只需把 STK 的 StkExternal 属性页 Filename
指向 `output\stk_ephemeris\sat_00x.e`，之后点 **Reload Ephemeris** +
**Apply** 即可刷新（不必重新建对象）。

> 注意：`Reload Ephemeris` 只是重新读盘。如果磁盘上还是旧格式文件，
> 刷新多少次都不会出现星历——务必先确认文件的头部符合上一节的样式。


---

## 九、实测验证结论（run_log.txt）

- 半长轴 a = 6978.137 km，平均运动 n = 1.08308×10⁻³ rad/s，轨道周期 5801.23 s
  （96.69 min），24 h 运行 14.89 圈，与理论值一致；
- 全部 1440×144 颗次位置模长恒为 a（最大偏差 1.86×10⁻⁹ m）；
- MATLAB 计算的速度模长恒为 7557.865 m/s = √(μ/a)（最大偏差 3.6×10⁻¹² m/s）；.e 仅写入位置，STK 从插值生成速度；
- vis/range 逐时隙对称、对角线正确；线段法与夹角解析判据逐对完全一致；
- 同轨相邻星（30°）可见、距离 3612.1 km；隔一颗（60°）不可见；
- 产物维度：pos_data [1440,144,3]、vis_data [1440,144,144] logical、
  range_data [1440,144,144] double；
- 3 个长表 CSV 各 207360 条记录（1440×144），144 个单星 trace 格式正确；
- 经纬高抽查：t=0、卫星 1 位于纬度 0°、经度 79.5394°E、椭球高 600 km，
  与 J2000 历元春分点地理经度的理论值一致。
