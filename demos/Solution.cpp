#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

struct Vector2D;
struct Polygon;

using namespace std;
const double EPS = 1e-6;

struct Vector2D {
    double x, y;

    Vector2D(double x = 0, double y = 0) : x(x), y(y) {}
    
    inline double Dot(const Vector2D& other) const { return x * other.x + y * other.y; }
    
    inline Vector2D operator+(const Vector2D& other) const { return Vector2D(x + other.x, y + other.y); }
    inline Vector2D operator-(const Vector2D& other) const { return Vector2D(x - other.x, y - other.y); }
    inline Vector2D operator*(double s) const { return Vector2D(x * s, y * s); }
    
    inline Vector2D Normalize() const {
        double len = std::sqrt(x * x + y * y);
        return (len > 0) ? Vector2D(x / len, y / len) : Vector2D(0, 0);
    }
    inline Vector2D Perp() const { return Vector2D(-y, x); }
};

struct Polygon {
    std::vector<Vector2D> vertices;
    Vector2D center; // 简化：直接使用 center，去掉 mutable 和 computed 标志位

    Polygon() = default;
    Polygon(std::initializer_list<Vector2D> vts) : vertices(vts) {}
};

struct Projection {
    double min, max;
};

struct ProjectionData {
    double centerDiff; // (CenterB - CenterA) 静态部分
    double sumRadii;   // RadiusA + RadiusB
    Vector2D axis;
};

// 优化：使用单一结构体数组存储预处理数据，提高缓存命中率
vector<ProjectionData> axesData; 

int n1 = 0, n2 = 0, m = 0;
Polygon polygon1;
Polygon polygon2;
vector<Vector2D> testCases;

// 优化：拆分 Projection 为独立的数组，提高缓存局部性
vector<Vector2D> poly1Axes;
vector<double> poly1MinOwn, poly1MaxOwn;       // Polygon1 在自己轴上的投影
vector<double> poly2MinOnPoly1, poly2MaxOnPoly1; // Polygon2 在 Polygon1 轴上的投影

vector<Vector2D> poly2Axes;
vector<double> poly2MinOwn, poly2MaxOwn;       // Polygon2 在自己轴上的投影
vector<double> poly1MinOnPoly2, poly1MaxOnPoly2; // Polygon1 在 Polygon2 轴上的投影

Vector2D GenSolution(const Vector2D& vec)
{
    double minOverlap = std::numeric_limits<double>::infinity();
    Vector2D smallestAxis(0, 0);
    
    const size_t totalAxes = axesData.size();
    const ProjectionData* data = axesData.data();

    // 预计算中心差值的静态部分已在 data[i].centerDiff 中
    // 运行时只需加上 vec 在轴上的投影
    
    // 优化：手动展开或使用指针遍历减少开销
    for (size_t i = 0; i < totalAxes; ++i) {
        double shift = vec.Dot(data[i].axis);
        double dist = data[i].centerDiff + shift;
        
        // 快速绝对值计算 (通常比 std::abs 快或内联)
        double absDist = (dist < 0.0) ? -dist : dist;
        
        if (absDist >= data[i].sumRadii) {
            return {0.0, 0.0}; // 发现分离轴
        }

        double overlap = data[i].sumRadii - absDist;
        if (overlap < minOverlap) {
            minOverlap = overlap;
            smallestAxis = data[i].axis;
        }
    }

    if (minOverlap == std::numeric_limits<double>::infinity()) {
         // 理论上不会发生，除非多边形退化且完全重合且无轴？
         return {0.0, 0.0};
    }

    // 优化：直接使用预存的 center 差值基础向量
    Vector2D dir = (polygon2.center - polygon1.center) + vec;
    
    if (smallestAxis.Dot(dir) < 0.0) {
        smallestAxis = smallestAxis * -1.0;
    }

    return smallestAxis * minOverlap;
}

void PreProcess()
{
    size_t n1Verts = polygon1.vertices.size();
    size_t n2Verts = polygon2.vertices.size();

    // 预计算中心
    if (n1Verts > 0) {
        polygon1.center = Vector2D(0, 0);
        for (const auto& v : polygon1.vertices) {
            polygon1.center = polygon1.center + v;
        }
        polygon1.center = polygon1.center * (1.0 / n1Verts);
    }

    if (n2Verts > 0) {
        polygon2.center = Vector2D(0, 0);
        for (const auto& v : polygon2.vertices) {
            polygon2.center = polygon2.center + v;
        }
        polygon2.center = polygon2.center * (1.0 / n2Verts);
    }

    const Vector2D* verts1 = polygon1.vertices.data();
    const Vector2D* verts2 = polygon2.vertices.data();

    axesData.clear();
    axesData.reserve(n1Verts + n2Verts);

    auto processAxis = [&](const Vector2D& axis) {
        // 计算 Polygon1 投影
        double val1 = verts1[0].Dot(axis);
        double minA = val1, maxA = val1;
        for (size_t k = 1; k < n1Verts; ++k) {
            val1 = verts1[k].Dot(axis);
            if (val1 < minA) minA = val1;
            else if (val1 > maxA) maxA = val1;
        }
        double centerA = (minA + maxA) * 0.5;
        double radiusA = (maxA - minA) * 0.5;

        // 计算 Polygon2 投影
        double val2 = verts2[0].Dot(axis);
        double minB = val2, maxB = val2;
        for (size_t k = 1; k < n2Verts; ++k) {
            val2 = verts2[k].Dot(axis);
            if (val2 < minB) minB = val2;
            else if (val2 > maxB) maxB = val2;
        }
        double centerB = (minB + maxB) * 0.5;
        double radiusB = (maxB - minB) * 0.5;

        ProjectionData pd;
        pd.axis = axis;
        pd.centerDiff = centerB - centerA; // 预计算静态中心差
        pd.sumRadii = radiusA + radiusB;   // 预计算半径和
        axesData.push_back(pd);
    };

    // --- 处理多边形 1 的轴 ---
    for (size_t i = 0; i < n1Verts; ++i) {
        Vector2D edge = verts1[(i + 1) % n1Verts] - verts1[i];
        processAxis(edge.Perp().Normalize());
    }

    // --- 处理多边形 2 的轴 ---
    for (size_t i = 0; i < n2Verts; ++i) {
        Vector2D edge = verts2[(i + 1) % n2Verts] - verts2[i];
        processAxis(edge.Perp().Normalize());
    }
    
    // 清理旧的大向量以释放内存，虽然它们会被重新分配，但显式清理有助于峰值内存控制
    poly1Axes.clear(); poly1Axes.shrink_to_fit();
    poly1MinOwn.clear(); poly1MinOwn.shrink_to_fit();
    poly1MaxOwn.clear(); poly1MaxOwn.shrink_to_fit();
    poly2MinOnPoly1.clear(); poly2MinOnPoly1.shrink_to_fit();
    poly2MaxOnPoly1.clear(); poly2MaxOnPoly1.shrink_to_fit();
    
    poly2Axes.clear(); poly2Axes.shrink_to_fit();
    poly2MinOwn.clear(); poly2MinOwn.shrink_to_fit();
    poly2MaxOwn.clear(); poly2MaxOwn.shrink_to_fit();
    poly1MinOnPoly2.clear(); poly1MinOnPoly2.shrink_to_fit();
    poly1MaxOnPoly2.clear(); poly1MaxOnPoly2.shrink_to_fit();
}

int main()
{
    // 优化：关闭同步，加速 IO
    std::ios::sync_with_stdio(false);
    cin.tie(nullptr);

    cin >> n1 >> n2;

    if (n1 <= 0 || n2 <= 0) {
        cerr << "Input data error: the number of vertices of both polygons should be greater than 2" << endl;
        return 1;
    }

    polygon1.vertices.resize(n1);
    for (int i = 0; i < n1; ++i) {
        cin >> polygon1.vertices[i].x >> polygon1.vertices[i].y;
    }

    polygon2.vertices.resize(n2);
    for (int i = 0; i < n2; ++i) {
        cin >> polygon2.vertices[i].x >> polygon2.vertices[i].y;
    }

    string okResp;
    cin >> okResp;
    if (okResp != "OK") {
        cerr << "Input data error: waiting for OK after obtaining polygons but I get " << okResp << endl;
        return 0;
    }

    PreProcess();
    
    cout << "OK" << "\n"; // 优化：使用 \n 代替 endl 避免频繁 flush
    cout.flush();

    cin >> m;
    testCases.resize(m);

    for (int i = 0; i < m; ++i) {
        cin >> testCases[i].x >> testCases[i].y;
    }

    cin >> okResp;
    if (okResp != "OK") {
        cerr << "Input data error: waiting for OK after that I have received all test points but I get " << okResp
             << endl;
        return 0;
    }

    for (int i = 0; i < m; ++i) {
        const Vector2D& res = GenSolution(testCases[i]);
        cout << fixed << setprecision(5) << res.x << " " << res.y << "\n";
    }

    cout << "OK" << "\n";
    cout.flush();

    return 0;
}