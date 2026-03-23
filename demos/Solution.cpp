#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include <string>

using namespace std;

const double EPS = 1e-8; // 稍微减小精度阈值以提高稳定性

struct Vec2 {
    double x, y;
    Vec2 operator+(const Vec2& o) const { return {x + o.x, y + o.y}; }
    Vec2 operator-(const Vec2& o) const { return {x - o.x, y - o.y}; }
    Vec2 operator*(double s) const { return {x * s, y * s}; }
    Vec2 operator/(double s) const { return {x / s, y / s}; }
    double dot(const Vec2& o) const { return x * o.x + y * o.y; }
    double cross(const Vec2& o) const { return x * o.y - y * o.x; }
    double normSq() const { return x * x + y * y; }
    double norm() const { return sqrt(normSq()); }
};

// 确保多边形点是逆时针顺序，并从最左下角的点开始
void reorder(vector<Vec2>& p) {
    if (p.empty()) return;
    
    // 1. 计算有向面积判断方向
    double area = 0;
    for (size_t i = 0; i < p.size(); ++i) {
        area += p[i].cross(p[(i + 1) % p.size()]);
    }
    // 如果面积为负，说明是顺时针，需要反转
    if (area < 0) {
        reverse(p.begin(), p.end());
    }
    
    // 2. 找到最左下角的点作为起点 (y 最小，若相同则 x 最小)
    int start = 0;
    for (int i = 1; i < (int)p.size(); ++i) {
        if (p[i].y < p[start].y - EPS || 
            (abs(p[i].y - p[start].y) < EPS && p[i].x < p[start].x - EPS)) {
            start = i;
        }
    }
    rotate(p.begin(), p.begin() + start, p.end());
}

struct SatAxis {
    Vec2 axis;
    double minProj;
    double maxProj;
};

vector<SatAxis> mAxes;
Vec2 mDiffCenter = {0, 0};
int n1, n2, m;

// 辅助函数：获取多边形的所有边向量（已排序）
vector<Vec2> getSortedEdges(const vector<Vec2>& p) {
    vector<Vec2> edges;
    if (p.empty()) return edges;
    for (size_t i = 0; i < p.size(); ++i) {
        Vec2 edge = p[(i + 1) % p.size()] - p[i];
        edges.push_back(edge);
    }
    return edges;
}

void PreProcess(vector<Vec2>& p1, vector<Vec2>& p2) {
    // 1. 标准化多边形：逆时针，起点一致
    reorder(p1);
    
    // 构建 -P2 (将 P2 关于原点对称)
    vector<Vec2> negP2;
    negP2.reserve(p2.size());
    for (const auto& v : p2) {
        negP2.push_back({-v.x, -v.y});
    }
    reorder(negP2);
    
    // 2. 构建闵可夫斯基差 A - B = A + (-B)
    // 由于 A 和 -B 都是凸多边形且已按逆时针排序，它们的边向量也是按极角排序的。
    // 我们可以使用双指针法合并边向量，从而在线性时间内构建结果多边形的顶点。
    
    vector<Vec2> e1 = getSortedEdges(p1);
    vector<Vec2> e2 = getSortedEdges(negP2);
    
    vector<Vec2> mDiffVerts;
    if (!p1.empty() && !negP2.empty()) {
        // 起始点：P1 的第一个点 + (-P2) 的第一个点
        Vec2 currentPt = p1[0] + negP2[0];
        mDiffVerts.push_back(currentPt);
        
        size_t i = 0, j = 0;
        while (i < e1.size() || j < e2.size()) {
            Vec2 nextEdge;
            bool takeE1 = false;
            
            if (i >= e1.size()) {
                nextEdge = e2[j++];
            } else if (j >= e2.size()) {
                nextEdge = e1[i++];
            } else {
                // 比较两个边向量的极角
                // 使用叉积判断：如果 e1 x e2 > 0，说明 e1 在 e2 的逆时针方向（更小或更大取决于坐标系，这里我们要按逆时针顺序合并）
                // 在标准数学坐标系中，叉积 > 0 表示 e1 到 e2 是逆时针旋转，即 e1 的极角小于 e2。
                // 我们希望按极角从小到大（逆时针）添加边。
                double cr = e1[i].cross(e2[j]);
                if (cr > EPS) { 
                    // e1 极角更小，先取 e1
                    nextEdge = e1[i++];
                } else if (cr < -EPS) {
                    // e2 极角更小，先取 e2
                    nextEdge = e2[j++];
                } else {
                    // 极角相同（平行），合并长度或直接取一个（凸多边形合并通常取长度和，或者直接依次添加）
                    // 这里为了保持顶点数最少，可以将同向边合并，也可以依次添加。
                    // 依次添加更简单且不会出错
                    nextEdge = e1[i++]; 
                    // 注意：这里如果不合并，下一轮循环还会比较剩下的那个，逻辑依然成立。
                    // 但为了严谨，如果方向完全一致，应该把两个向量都加进去，或者加起来加一次。
                    // 简单的做法是：如果叉积为 0，任意取一个，下一次循环另一个还会被比较。
                    // 但为了防止死循环或漏掉，我们这里采取：如果平行，两个都取（先取 e1，下次循环 e2 还会和新的 e1 比，或者和旧的 e1 比？不对）
                    // 修正：如果平行，应该同时前进吗？不，闵可夫斯基和的边是两组边的并集排序。
                    // 如果角度相同，它们在排序序列中是相邻的。
                    // 所以这里如果相等，取 e1，i++。下一轮循环，e1[i] (新) 和 e2[j] (旧) 比。
                    // 这样没问题。
                }
            }
            
            currentPt = currentPt + nextEdge;
            mDiffVerts.push_back(currentPt);
        }
        
        // 去除首尾重复点（如果闭合）
        if (mDiffVerts.size() > 1) {
            if (abs(mDiffVerts.front().x - mDiffVerts.back().x) < EPS && 
                abs(mDiffVerts.front().y - mDiffVerts.back().y) < EPS) {
                mDiffVerts.pop_back();
            }
        }
    }

    // 3. 计算闵可夫斯基差的中心
    mDiffCenter = {0, 0};
    if (!mDiffVerts.empty()) {
        for(const auto& v : mDiffVerts) {
            mDiffCenter = mDiffCenter + v;
        }
        mDiffCenter = mDiffCenter / (double)mDiffVerts.size();
    }

    // 4. 预处理所有分离轴（法线）及其投影区间
    mAxes.clear();
    int n = mDiffVerts.size();
    if (n == 0) return;

    for (int k = 0; k < n; ++k) {
        Vec2 p1_v = mDiffVerts[k];
        Vec2 p2_v = mDiffVerts[(k + 1) % n];
        Vec2 edge = p2_v - p1_v;
        
        // 跳过极短的边
        if (edge.normSq() < EPS * EPS) continue;
        
        // 法线指向多边形外部还是内部？
        // 对于逆时针多边形，边向量 (p1->p2) 的左手法线 (-dy, dx) 指向多边形外部。
        // 但 SAT 测试中，我们需要的是投影轴。方向不重要，只要统一即可。
        // 这里我们取垂直于边的单位向量。
        Vec2 axis = {-edge.y, edge.x};
        double len = axis.norm();
        if (len < EPS) continue;
        axis = axis / len;
        
        // 计算整个闵可夫斯基差多边形在该轴上的投影区间
        double minP = 1e18;
        double maxP = -1e18;
        for (const auto& v : mDiffVerts) {
            double proj = v.dot(axis);
            if (proj < minP) minP = proj;
            if (proj > maxP) maxP = proj;
        }
        mAxes.push_back({axis, minP, maxP});
    }
}

Vec2 GenSolution(Vec2 v) {
    double minOverlap = 1e18;
    Vec2 bestAxis = {0, 0};
    bool foundOverlap = false;
    
    for (const auto& ax : mAxes) {
        double projV = v.dot(ax.axis);
        
        // 如果点在任何轴的投影区间之外，说明没有重叠，MTV 为 0
        // 考虑到浮点误差，使用 EPS 缓冲
        if (projV > ax.maxProj + EPS || projV < ax.minProj - EPS) {
            return {0, 0};
        }
        
        // 计算点在区间内的重叠深度
        // 点在区间内，距离左边界的距离是 projV - minP，距离右边界是 maxP - projV
        // 要将点移出区间，需要的最小位移是 min(到左边界的距离，到右边界的距离)
        // 注意：这里的“重叠”是指点侵入多边形的深度。
        // 如果点在多边形内，projV 必然在 [minP, maxP] 之间。
        // 最小推出距离 = min(maxP - projV, projV - minP)
        
        double distToMax = ax.maxProj - projV;
        double distToMin = projV - ax.minProj;
        
        double overlap = 0.0;
        Vec2 candidateAxis;
        
        if (distToMax < distToMin) {
            overlap = distToMax;
            candidateAxis = ax.axis; // 向右推 (沿 axis 正方向)
        } else {
            overlap = distToMin;
            candidateAxis = ax.axis * (-1.0); // 向左推 (沿 axis 负方向)
        }
        
        // 更新最小重叠
        if (overlap < minOverlap) {
            minOverlap = overlap;
            bestAxis = candidateAxis;
            foundOverlap = true;
        }
    }
    
    if (!foundOverlap || minOverlap < EPS) {
        return {0, 0};
    }
    
    // 方向修正：
    // 理论上，上述逻辑已经根据 distToMax/distToMin 选择了将点推向外部的最短方向。
    // 但是，为了确保物理意义上的正确性（即从多边形 A 指向多边形 B 的推力，或者说将 B 移出 A 的方向），
    // 我们通常需要验证方向是否指向多边形外部。
    // 在闵可夫斯基差空间中，原点代表重合。点 v 代表相对位置。
    // 如果 v 在 M(A-B) 内部，我们需要将 v 移到边界。
    // 上述逻辑计算的 bestAxis 是将 v 推向最近边界的方向。
    // 这个方向就是 MTV 的方向。
    // 之前的代码有一个额外的检查：if (bestAxis.dot(dir) < 0) flip。
    // dir = v - center。如果 bestAxis 和 (v - center) 夹角大于 90 度，说明 bestAxis 指向中心，这是不对的，应该背离中心。
    // 让我们保留这个安全检查，因为它能纠正某些极端情况下的法线方向选择错误。
    
    Vec2 dir = v - mDiffCenter;
    // 如果计算出的最佳轴方向指向多边形内部（与从中心到点的向量相反），则翻转
    // 注意：如果 v 就在中心附近，dir 接近 0，点积可能不准，但此时 overlap 也很小，影响不大
    if (bestAxis.dot(dir) < 0.0) {
        bestAxis = bestAxis * (-1.0);
    }
    
    return bestAxis * minOverlap;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    
    if (!(cin >> n1 >> n2)) return 0;
    
    vector<Vec2> p1(n1), p2(n2);
    for (int i = 0; i < n1; ++i) cin >> p1[i].x >> p1[i].y;
    for (int i = 0; i < n2; ++i) cin >> p2[i].x >> p2[i].y;
    
    string ok;
    cin >> ok; 
    
    PreProcess(p1, p2);
    
    cout << "OK" << endl;
    
    if (!(cin >> m)) return 0;
    
    vector<Vec2> queries(m);
    for (int i = 0; i < m; ++i) {
        cin >> queries[i].x >> queries[i].y;
    }
    cin >> ok; 
    
    cout << fixed << setprecision(5);
    for (int i = 0; i < m; ++i) {
        Vec2 res = GenSolution(queries[i]);
        cout << res.x << " " << res.y << "\n";
    }
    cout << "OK" << endl;
    
    return 0;
}