#pragma GCC optimize("O3")
#include <cstdio>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <algorithm>

using namespace std;

const double EPS = 1e-9;
const double INF = 1e100;

int n1, n2, m;
double* A;   // 存储顶点坐标，顺序：x1,y1,x2,y2,...
double* B;   // 存储原始 B 顶点坐标

// 新增：MTV 计算计数器
int mtv_count = 0;

// 包围盒
double a_minx, a_maxx, a_miny, a_maxy;
double b_minx, b_maxx, b_miny, b_maxy;

// 轴数据（最多 n1+n2 条轴）
int axis_cnt;
double axis_dir_x[2000], axis_dir_y[2000];
double axis_a_min[2000], axis_a_max[2000];
double axis_b_min[2000], axis_b_max[2000];

// 投影函数（内联）
static inline void projectPolygon(const double* poly, int n, double dir_x, double dir_y,
                                  double& minProj, double& maxProj) {
    minProj = INF;
    maxProj = -INF;
    for (int i = 0; i < n; ++i) {
        double proj = poly[2*i] * dir_x + poly[2*i+1] * dir_y;
        if (proj < minProj) minProj = proj;
        if (proj > maxProj) maxProj = proj;
    }
}

// 判断两个单位向量是否平行（方向相同或相反）
static inline bool isParallel(double x1, double y1, double x2, double y2) {
    double dot = x1*x2 + y1*y2;
    return fabs(fabs(dot) - 1.0) < EPS;
}

void PreProcess() {
    // 计算A包围盒
    a_minx = a_miny = INF;
    a_maxx = a_maxy = -INF;
    for (int i = 0; i < n1; ++i) {
        double x = A[2*i], y = A[2*i+1];
        if (x < a_minx) a_minx = x;
        if (x > a_maxx) a_maxx = x;
        if (y < a_miny) a_miny = y;
        if (y > a_maxy) a_maxy = y;
    }
    // 计算B包围盒
    b_minx = b_miny = INF;
    b_maxx = b_maxy = -INF;
    for (int i = 0; i < n2; ++i) {
        double x = B[2*i], y = B[2*i+1];
        if (x < b_minx) b_minx = x;
        if (x > b_maxx) b_maxx = x;
        if (y < b_miny) b_miny = y;
        if (y > b_maxy) b_maxy = y;
    }

    // 收集唯一轴方向
    double unique_x[2000], unique_y[2000];
    int uniq_cnt = 0;

    auto addUnique = [&](double dx, double dy) {
        // 归一化
        double len = sqrt(dx*dx + dy*dy);
        if (len < EPS) return;
        dx /= len; dy /= len;
        for (int i = 0; i < uniq_cnt; ++i) {
            if (isParallel(unique_x[i], unique_y[i], dx, dy)) return;
        }
        unique_x[uniq_cnt] = dx;
        unique_y[uniq_cnt] = dy;
        uniq_cnt++;
    };

    // A的边法线
    for (int i = 0; i < n1; ++i) {
        int i1 = i, i2 = (i+1) % n1;
        double dx = A[2*i2] - A[2*i1];
        double dy = A[2*i2+1] - A[2*i1+1];
        // 法线是 (-dy, dx) 或 (dy, -dx)，这里取 (-dy, dx)
        addUnique(-dy, dx);
    }

    // B的边法线
    for (int i = 0; i < n2; ++i) {
        int i1 = i, i2 = (i+1) % n2;
        double dx = B[2*i2] - B[2*i1];
        double dy = B[2*i2+1] - B[2*i1+1];
        addUnique(-dy, dx);
    }

    // 为每个唯一轴计算投影区间
    axis_cnt = uniq_cnt;
    for (int i = 0; i < axis_cnt; ++i) {
        double dx = unique_x[i], dy = unique_y[i];
        axis_dir_x[i] = dx;
        axis_dir_y[i] = dy;

        double minA, maxA, minB, maxB;
        projectPolygon(A, n1, dx, dy, minA, maxA);
        projectPolygon(B, n2, dx, dy, minB, maxB);
        axis_a_min[i] = minA;
        axis_a_max[i] = maxA;
        axis_b_min[i] = minB;
        axis_b_max[i] = maxB;
    }
}

// 计算最小平移向量
static inline void genSolution(double vec_x, double vec_y, double& out_x, double& out_y) {
    // 快速包围盒检测
    double b_minx_t = b_minx + vec_x;
    double b_maxx_t = b_maxx + vec_x;
    double b_miny_t = b_miny + vec_y;
    double b_maxy_t = b_maxy + vec_y;
    if (a_maxx < b_minx_t - EPS || a_minx > b_maxx_t + EPS ||
        a_maxy < b_miny_t - EPS || a_miny > b_maxy_t + EPS) {
        out_x = out_y = 0.0;
        return;
    }

    double bestDist = INF;
    double bestVecX = 0.0, bestVecY = 0.0;

    for (int i = 0; i < axis_cnt; ++i) {
        double dx = axis_dir_x[i];
        double dy = axis_dir_y[i];
        double dot = dx * vec_x + dy * vec_y;
        double b_min = axis_b_min[i] + dot;
        double b_max = axis_b_max[i] + dot;

        if (axis_a_max[i] < b_min - EPS || b_max < axis_a_min[i] - EPS) {
            out_x = out_y = 0.0;
            return;
        }

        double dist_right = axis_a_max[i] - b_min;
        double dist_left  = b_max - axis_a_min[i];
        double dist = dist_right < dist_left ? dist_right : dist_left;
        if (dist < bestDist - EPS) {
            bestDist = dist;
            if (dist_right <= dist_left) {
                bestVecX = dx * dist;
                bestVecY = dy * dist;
            } else {
                bestVecX = dx * (-dist);
                bestVecY = dy * (-dist);
            }
        }
    }
    out_x = bestVecX;
    out_y = bestVecY;
    
    // 如果计算出了有效的移动向量（距离大于 0），则计数加一
    if (bestDist < INF && bestDist > EPS) {
        mtv_count++;
    }
}

int main() {
    // 读入多边形
    scanf("%d%d", &n1, &n2);
    A = (double*)malloc(2 * n1 * sizeof(double));
    B = (double*)malloc(2 * n2 * sizeof(double));
    for (int i = 0; i < n1; ++i) {
        scanf("%lf%lf", &A[2*i], &A[2*i+1]);
    }
    for (int i = 0; i < n2; ++i) {
        scanf("%lf%lf", &B[2*i], &B[2*i+1]);
    }

    // 读取OK
    char ok[4];
    scanf("%s", ok);
    if (strcmp(ok, "OK") != 0) return 0;

    PreProcess();
    printf("OK\n");
    fflush(stdout);


    // 读取测试用例数量
    scanf("%d", &m);
    double* testX = (double*)malloc(m * sizeof(double));
    double* testY = (double*)malloc(m * sizeof(double));
    for (int i = 0; i < m; ++i) {
        scanf("%lf%lf", &testX[i], &testY[i]);
    }
    scanf("%s", ok);
    if (strcmp(ok, "OK") != 0) return 0;

    // 输出结果
    for (int i = 0; i < m; ++i) {
        double ans_x, ans_y;
        genSolution(testX[i], testY[i], ans_x, ans_y);
        printf("%.5f %.5f\n", ans_x, ans_y);
        fflush(stdout);
    }

    printf("OK\n");
    fflush(stdout);

    printf("%d\n", mtv_count);
    fflush(stdout);

    free(A); free(B); free(testX); free(testY);
    return 0;
}