#pragma GCC optimize("O3")
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

using namespace std;

const double GEOM_EPS = 1e-12;
const double QUERY_EPS = 1e-10;
const double PRINT_EPS = 5e-8;
const double INF = 1e100;
const double OUT_STEP = 1e-7;
const double OUTPUT_MARGIN = 3e-5;

struct Vec {
    double x, y;

    Vec() : x(0.0), y(0.0) {}
    Vec(double x_, double y_) : x(x_), y(y_) {}

    Vec operator + (const Vec& other) const { return Vec(x + other.x, y + other.y); }
    Vec operator - (const Vec& other) const { return Vec(x - other.x, y - other.y); }
    Vec operator * (double k) const { return Vec(x * k, y * k); }
};

struct Segment {
    Vec a, b;
};

struct BoundarySegment {
    Segment seg;
    double minx, maxx, miny, maxy;
};

struct ConvexRegion {
    vector<Vec> pt;
    vector<Segment> edges;
    vector<double> nx, ny, c;
    double minx, maxx, miny, maxy;
};

static inline double dot(const Vec& a, const Vec& b) {
    return a.x * b.x + a.y * b.y;
}

static inline double cross(const Vec& a, const Vec& b) {
    return a.x * b.y - a.y * b.x;
}

static inline double norm2(const Vec& v) {
    return dot(v, v);
}

static inline double length(const Vec& v) {
    return sqrt(norm2(v));
}

static inline bool samePoint(const Vec& a, const Vec& b) {
    return fabs(a.x - b.x) <= GEOM_EPS && fabs(a.y - b.y) <= GEOM_EPS;
}

static inline bool sameDirection(const Vec& a, const Vec& b) {
    double la2 = norm2(a);
    double lb2 = norm2(b);
    if (la2 <= GEOM_EPS || lb2 <= GEOM_EPS) {
        return true;
    }
    double limit = GEOM_EPS * sqrt(la2 * lb2);
    return fabs(cross(a, b)) <= limit && dot(a, b) >= 0.0;
}

static inline double signedArea2(const vector<Vec>& poly) {
    double area2 = 0.0;
    int n = static_cast<int>(poly.size());
    for (int i = 0; i < n; ++i) {
        area2 += cross(poly[i], poly[(i + 1) % n]);
    }
    return area2;
}

static inline void normalizeCCW(vector<Vec>& poly) {
    if (signedArea2(poly) < 0.0) {
        reverse(poly.begin(), poly.end());
    }
}

static inline bool isConvexPolygon(const vector<Vec>& poly) {
    int n = static_cast<int>(poly.size());
    if (n < 3) {
        return false;
    }
    int sign = 0;
    for (int i = 0; i < n; ++i) {
        Vec e1 = poly[(i + 1) % n] - poly[i];
        Vec e2 = poly[(i + 2) % n] - poly[(i + 1) % n];
        double cr = cross(e1, e2);
        if (fabs(cr) <= GEOM_EPS) {
            continue;
        }
        int cur = cr > 0.0 ? 1 : -1;
        if (sign == 0) {
            sign = cur;
        } else if (sign != cur) {
            return false;
        }
    }
    return true;
}

static inline bool pointInTriangleInclusive(const Vec& p, const Vec& a, const Vec& b, const Vec& c) {
    double c1 = cross(b - a, p - a);
    double c2 = cross(c - b, p - b);
    double c3 = cross(a - c, p - c);
    return c1 >= -GEOM_EPS && c2 >= -GEOM_EPS && c3 >= -GEOM_EPS;
}

static vector<vector<Vec>> triangulateSimplePolygon(const vector<Vec>& polyInput) {
    vector<Vec> poly = polyInput;
    normalizeCCW(poly);
    int n = static_cast<int>(poly.size());

    vector<vector<Vec>> triangles;
    if (n == 3) {
        triangles.push_back(poly);
        return triangles;
    }

    vector<int> prev(n), next(n);
    vector<char> alive(n, 1);
    for (int i = 0; i < n; ++i) {
        prev[i] = (i - 1 + n) % n;
        next[i] = (i + 1) % n;
    }

    auto isEar = [&](int idx) -> bool {
        if (!alive[idx]) {
            return false;
        }
        int a = prev[idx];
        int c = next[idx];
        if (!alive[a] || !alive[c]) {
            return false;
        }
        if (cross(poly[idx] - poly[a], poly[c] - poly[idx]) <= GEOM_EPS) {
            return false;
        }
        for (int j = next[c]; j != a; j = next[j]) {
            if (!alive[j] || j == idx) {
                continue;
            }
            if (pointInTriangleInclusive(poly[j], poly[a], poly[idx], poly[c])) {
                return false;
            }
        }
        return true;
    };

    int aliveCount = n;
    int start = 0;
    while (aliveCount > 3) {
        bool clipped = false;
        int idx = start;
        for (int steps = 0; steps < n; ++steps) {
            if (alive[idx] && isEar(idx)) {
                int a = prev[idx];
                int c = next[idx];
                triangles.push_back(vector<Vec>{poly[a], poly[idx], poly[c]});
                alive[idx] = 0;
                next[a] = c;
                prev[c] = a;
                start = c;
                --aliveCount;
                clipped = true;
                break;
            }
            idx = next[idx];
        }
        if (!clipped) {
            vector<int> remain;
            remain.reserve(aliveCount);
            int first = -1;
            for (int i = 0; i < n; ++i) {
                if (alive[i]) {
                    first = i;
                    break;
                }
            }
            if (first == -1) {
                break;
            }
            int cur = first;
            do {
                remain.push_back(cur);
                cur = next[cur];
            } while (cur != first);
            for (size_t i = 1; i + 1 < remain.size(); ++i) {
                triangles.push_back(vector<Vec>{poly[remain[0]], poly[remain[i]], poly[remain[i + 1]]});
            }
            return triangles;
        }
    }

    int first = -1;
    for (int i = 0; i < n; ++i) {
        if (alive[i]) {
            first = i;
            break;
        }
    }
    if (first != -1) {
        int second = next[first];
        int third = next[second];
        triangles.push_back(vector<Vec>{poly[first], poly[second], poly[third]});
    }
    return triangles;
}

static int lowestVertexIndex(const vector<Vec>& poly) {
    int best = 0;
    for (int i = 1, n = static_cast<int>(poly.size()); i < n; ++i) {
        if (poly[i].y < poly[best].y - GEOM_EPS ||
            (fabs(poly[i].y - poly[best].y) <= GEOM_EPS && poly[i].x < poly[best].x - GEOM_EPS)) {
            best = i;
        }
    }
    return best;
}

static vector<Vec> reorderPolygon(const vector<Vec>& poly) {
    int n = static_cast<int>(poly.size());
    int start = lowestVertexIndex(poly);
    vector<Vec> res;
    res.reserve(n);
    for (int i = 0; i < n; ++i) {
        res.push_back(poly[(start + i) % n]);
    }
    return res;
}

static vector<Vec> compressConvexPolygon(vector<Vec> poly) {
    while (poly.size() > 1 && samePoint(poly.front(), poly.back())) {
        poly.pop_back();
    }
    if (poly.size() <= 2) {
        return poly;
    }

    bool changed = true;
    while (changed && poly.size() > 2) {
        changed = false;
        vector<Vec> next;
        next.reserve(poly.size());
        int n = static_cast<int>(poly.size());
        for (int i = 0; i < n; ++i) {
            const Vec& prev = poly[(i - 1 + n) % n];
            const Vec& cur = poly[i];
            const Vec& nxt = poly[(i + 1) % n];
            if (samePoint(prev, cur) || samePoint(cur, nxt)) {
                changed = true;
                continue;
            }
            Vec e1 = cur - prev;
            Vec e2 = nxt - cur;
            if (sameDirection(e1, e2)) {
                changed = true;
                continue;
            }
            next.push_back(cur);
        }
        poly.swap(next);
    }
    return poly;
}

static vector<Vec> cleanupPolygon(vector<Vec> poly) {
    vector<Vec> clean;
    clean.reserve(poly.size());
    for (size_t i = 0; i < poly.size(); ++i) {
        if (clean.empty() || !samePoint(clean.back(), poly[i])) {
            clean.push_back(poly[i]);
        }
    }
    if (clean.size() > 1 && samePoint(clean.front(), clean.back())) {
        clean.pop_back();
    }
    normalizeCCW(clean);
    return compressConvexPolygon(clean);
}

static bool tryMergeConvexPolygons(const vector<Vec>& a, const vector<Vec>& b, vector<Vec>& merged) {
    int na = static_cast<int>(a.size());
    int nb = static_cast<int>(b.size());
    for (int i = 0; i < na; ++i) {
        int ni = (i + 1) % na;
        for (int j = 0; j < nb; ++j) {
            int nj = (j + 1) % nb;
            if (!samePoint(a[i], b[nj]) || !samePoint(a[ni], b[j])) {
                continue;
            }

            vector<Vec> poly;
            poly.reserve(na + nb - 2);
            int cur = ni;
            poly.push_back(a[cur]);
            while (cur != i) {
                cur = (cur + 1) % na;
                poly.push_back(a[cur]);
            }
            cur = (j + 2) % nb;
            while (cur != j) {
                poly.push_back(b[cur]);
                cur = (cur + 1) % nb;
            }

            poly = cleanupPolygon(poly);
            if (poly.size() >= 3 && isConvexPolygon(poly) && fabs(signedArea2(poly)) > GEOM_EPS) {
                merged.swap(poly);
                return true;
            }
        }
    }
    return false;
}

static vector<vector<Vec> > mergeConvexParts(vector<vector<Vec> > parts) {
    bool changed = true;
    while (changed) {
        changed = false;
        for (size_t i = 0; i < parts.size() && !changed; ++i) {
            for (size_t j = i + 1; j < parts.size(); ++j) {
                vector<Vec> merged;
                if (tryMergeConvexPolygons(parts[i], parts[j], merged)) {
                    parts[i].swap(merged);
                    parts.erase(parts.begin() + j);
                    changed = true;
                    break;
                }
            }
        }
    }
    return parts;
}

static vector<Vec> minkowskiSum(const vector<Vec>& pRaw, const vector<Vec>& qRaw) {
    vector<Vec> p = reorderPolygon(pRaw);
    vector<Vec> q = reorderPolygon(qRaw);
    const int n = static_cast<int>(p.size());
    const int m = static_cast<int>(q.size());

    p.push_back(p[0]);
    p.push_back(p[1]);
    q.push_back(q[0]);
    q.push_back(q[1]);

    vector<Vec> res;
    res.reserve(n + m + 1);
    Vec cur = p[0] + q[0];
    res.push_back(cur);

    int i = 0;
    int j = 0;
    while (i < n || j < m) {
        if (i == n) {
            cur = cur + (q[j + 1] - q[j]);
            ++j;
        } else if (j == m) {
            cur = cur + (p[i + 1] - p[i]);
            ++i;
        } else {
            Vec ep = p[i + 1] - p[i];
            Vec eq = q[j + 1] - q[j];
            double limit = GEOM_EPS * sqrt(norm2(ep) * norm2(eq));
            double cr = cross(ep, eq);
            if (cr > limit) {
                cur = cur + ep;
                ++i;
            } else if (cr < -limit) {
                cur = cur + eq;
                ++j;
            } else {
                cur = cur + ep + eq;
                ++i;
                ++j;
            }
        }
        res.push_back(cur);
    }

    if (res.size() > 1 && samePoint(res.front(), res.back())) {
        res.pop_back();
    }
    return compressConvexPolygon(res);
}

static ConvexRegion buildConvexRegion(vector<Vec> poly) {
    normalizeCCW(poly);
    poly = compressConvexPolygon(poly);

    ConvexRegion region;
    region.pt = poly;
    int n = static_cast<int>(poly.size());
    region.nx.resize(n);
    region.ny.resize(n);
    region.c.resize(n);
    region.edges.resize(n);

    region.minx = region.maxx = poly[0].x;
    region.miny = region.maxy = poly[0].y;
    for (int i = 0; i < n; ++i) {
        const Vec& a = poly[i];
        const Vec& b = poly[(i + 1) % n];
        Vec edge = b - a;
        double len = length(edge);
        double ox = edge.y / len;
        double oy = -edge.x / len;
        region.nx[i] = ox;
        region.ny[i] = oy;
        region.c[i] = ox * a.x + oy * a.y;
        region.edges[i] = Segment{a, b};
        if (a.x < region.minx) region.minx = a.x;
        if (a.x > region.maxx) region.maxx = a.x;
        if (a.y < region.miny) region.miny = a.y;
        if (a.y > region.maxy) region.maxy = a.y;
    }
    return region;
}

static inline bool insideBox(const ConvexRegion& region, double x, double y) {
    return x >= region.minx - QUERY_EPS && x <= region.maxx + QUERY_EPS &&
           y >= region.miny - QUERY_EPS && y <= region.maxy + QUERY_EPS;
}

static inline bool pointInConvexStrict(const ConvexRegion& region, double x, double y) {
    if (!insideBox(region, x, y)) {
        return false;
    }
    for (size_t i = 0; i < region.pt.size(); ++i) {
        double slack = region.c[i] - (region.nx[i] * x + region.ny[i] * y);
        if (slack <= QUERY_EPS) {
            return false;
        }
    }
    return true;
}

static inline bool pointInConvexInclusive(const ConvexRegion& region, double x, double y) {
    if (!insideBox(region, x, y)) {
        return false;
    }
    for (size_t i = 0; i < region.pt.size(); ++i) {
        double slack = region.c[i] - (region.nx[i] * x + region.ny[i] * y);
        if (slack < -QUERY_EPS) {
            return false;
        }
    }
    return true;
}

static inline Vec closestPointOnSegment(const Vec& p, const Segment& seg) {
    Vec ab = seg.b - seg.a;
    double len2 = norm2(ab);
    if (len2 <= GEOM_EPS) {
        return seg.a;
    }
    double t = dot(p - seg.a, ab) / len2;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    return seg.a + ab * t;
}

static inline bool pointInsideAnyRegionStrict(const vector<ConvexRegion>& regions, double x, double y) {
    for (size_t i = 0; i < regions.size(); ++i) {
        if (pointInConvexStrict(regions[i], x, y)) {
            return true;
        }
    }
    return false;
}

static inline double clamp01(double t) {
    if (t < 0.0) return 0.0;
    if (t > 1.0) return 1.0;
    return t;
}

static inline void appendSplitParam(vector<double>& params, double t) {
    params.push_back(clamp01(t));
}

static inline void collectEdgeSplitParams(const Segment& base, const Segment& other, vector<double>& params) {
    Vec u = base.b - base.a;
    Vec v = other.b - other.a;
    double uu = norm2(u);
    double vv = norm2(v);
    if (uu <= GEOM_EPS || vv <= GEOM_EPS) {
        return;
    }

    Vec w = other.a - base.a;
    double den = cross(u, v);
    double limit = GEOM_EPS * (sqrt(uu * vv) + 1.0);
    if (fabs(den) > limit) {
        double alpha = cross(w, v) / den;
        double beta = cross(w, u) / den;
        if (alpha >= -GEOM_EPS && alpha <= 1.0 + GEOM_EPS &&
            beta >= -GEOM_EPS && beta <= 1.0 + GEOM_EPS) {
            appendSplitParam(params, alpha);
        }
        return;
    }

    double col1 = fabs(cross(other.a - base.a, u));
    double col2 = fabs(cross(other.b - base.a, u));
    double colLimit = GEOM_EPS * (uu + 1.0);
    if (col1 > colLimit || col2 > colLimit) {
        return;
    }

    appendSplitParam(params, dot(other.a - base.a, u) / uu);
    appendSplitParam(params, dot(other.b - base.a, u) / uu);
}

vector<Vec> polyA;
vector<Vec> polyB;
bool useConvexFastPath;

vector<double> fastNormalX;
vector<double> fastNormalY;
vector<double> fastOffsetC;
double fastMinX, fastMaxX, fastMinY, fastMaxY;

vector<ConvexRegion> generalRegions;
vector<BoundarySegment> generalBoundarySegments;
vector<vector<int> > generalRegionBuckets;
bool generalUseBucketsForQuery;
int generalGridW, generalGridH;
double generalMinX, generalMaxX, generalMinY, generalMaxY;
double generalCellW, generalCellH;

static inline int clampIndex(int v, int limit) {
    if (v < 0) return 0;
    if (v >= limit) return limit - 1;
    return v;
}

static inline int generalCellX(double x) {
    if (generalGridW <= 1 || generalCellW <= GEOM_EPS) {
        return 0;
    }
    return clampIndex(static_cast<int>((x - generalMinX) / generalCellW), generalGridW);
}

static inline int generalCellY(double y) {
    if (generalGridH <= 1 || generalCellH <= GEOM_EPS) {
        return 0;
    }
    return clampIndex(static_cast<int>((y - generalMinY) / generalCellH), generalGridH);
}

static void buildGeneralRegionBuckets() {
    generalRegionBuckets.clear();
    generalUseBucketsForQuery = false;
    if (generalRegions.empty()) {
        generalGridW = generalGridH = 1;
        generalCellW = generalCellH = 1.0;
        return;
    }

    double width = max(generalMaxX - generalMinX, 1e-6);
    double height = max(generalMaxY - generalMinY, 1e-6);
    int target = static_cast<int>(sqrt(static_cast<double>(generalRegions.size()) * 4.0));
    if (target < 4) target = 4;
    if (target > 32) target = 32;

    double scale = max(width, height);
    generalGridW = max(1, static_cast<int>(round(target * width / scale)));
    generalGridH = max(1, static_cast<int>(round(target * height / scale)));
    generalCellW = width / generalGridW;
    generalCellH = height / generalGridH;
    generalRegionBuckets.assign(generalGridW * generalGridH, vector<int>());

    long long totalBucketRefs = 0;
    for (size_t i = 0; i < generalRegions.size(); ++i) {
        int x0 = generalCellX(generalRegions[i].minx - QUERY_EPS);
        int x1 = generalCellX(generalRegions[i].maxx + QUERY_EPS);
        int y0 = generalCellY(generalRegions[i].miny - QUERY_EPS);
        int y1 = generalCellY(generalRegions[i].maxy + QUERY_EPS);
        for (int y = y0; y <= y1; ++y) {
            for (int x = x0; x <= x1; ++x) {
                generalRegionBuckets[y * generalGridW + x].push_back(static_cast<int>(i));
                ++totalBucketRefs;
            }
        }
    }

    double avgBucket = static_cast<double>(totalBucketRefs) / generalRegionBuckets.size();
    generalUseBucketsForQuery = avgBucket * 3.0 < static_cast<double>(generalRegions.size());
}

static inline bool pointInsideBucketedGeneralRegionStrict(double x, double y) {
    if (x < generalMinX - QUERY_EPS || x > generalMaxX + QUERY_EPS ||
        y < generalMinY - QUERY_EPS || y > generalMaxY + QUERY_EPS) {
        return false;
    }
    if (generalRegionBuckets.empty()) {
        return pointInsideAnyRegionStrict(generalRegions, x, y);
    }
    const vector<int>& bucket = generalRegionBuckets[generalCellY(y) * generalGridW + generalCellX(x)];
    for (size_t i = 0; i < bucket.size(); ++i) {
        if (pointInConvexStrict(generalRegions[bucket[i]], x, y)) {
            return true;
        }
    }
    return false;
}

static inline double boundaryBBoxDist2(const BoundarySegment& seg, double x, double y) {
    double dx = 0.0;
    double dy = 0.0;
    if (x < seg.minx) dx = seg.minx - x;
    else if (x > seg.maxx) dx = x - seg.maxx;
    if (y < seg.miny) dy = seg.miny - y;
    else if (y > seg.maxy) dy = y - seg.maxy;
    return dx * dx + dy * dy;
}

static inline BoundarySegment makeBoundarySegment(const Vec& a, const Vec& b) {
    BoundarySegment seg;
    seg.seg = Segment{a, b};
    seg.minx = min(a.x, b.x);
    seg.maxx = max(a.x, b.x);
    seg.miny = min(a.y, b.y);
    seg.maxy = max(a.y, b.y);
    return seg;
}

static inline bool sameBoundarySegment(const BoundarySegment& a, const BoundarySegment& b) {
    return (samePoint(a.seg.a, b.seg.a) && samePoint(a.seg.b, b.seg.b)) ||
           (samePoint(a.seg.a, b.seg.b) && samePoint(a.seg.b, b.seg.a));
}

static bool tryMergeBoundarySegments(const BoundarySegment& a, const BoundarySegment& b, BoundarySegment& merged) {
    if (sameBoundarySegment(a, b)) {
        merged = a;
        return true;
    }

    Vec shared, pa, pb;
    bool found = false;
    if (samePoint(a.seg.a, b.seg.a)) {
        shared = a.seg.a; pa = a.seg.b; pb = b.seg.b; found = true;
    } else if (samePoint(a.seg.a, b.seg.b)) {
        shared = a.seg.a; pa = a.seg.b; pb = b.seg.a; found = true;
    } else if (samePoint(a.seg.b, b.seg.a)) {
        shared = a.seg.b; pa = a.seg.a; pb = b.seg.b; found = true;
    } else if (samePoint(a.seg.b, b.seg.b)) {
        shared = a.seg.b; pa = a.seg.a; pb = b.seg.a; found = true;
    }
    if (!found) {
        return false;
    }

    Vec va = pa - shared;
    Vec vb = pb - shared;
    double la2 = norm2(va);
    double lb2 = norm2(vb);
    if (la2 <= GEOM_EPS || lb2 <= GEOM_EPS) {
        merged = la2 >= lb2 ? a : b;
        return true;
    }

    double limit = GEOM_EPS * sqrt(la2 * lb2);
    if (fabs(cross(va, vb)) > limit) {
        return false;
    }

    if (dot(va, vb) >= 0.0) {
        merged = makeBoundarySegment(shared, la2 >= lb2 ? pa : pb);
    } else {
        merged = makeBoundarySegment(pa, pb);
    }
    return true;
}

static void mergeGeneralBoundarySegments() {
    bool changed = true;
    while (changed) {
        changed = false;
        for (size_t i = 0; i < generalBoundarySegments.size() && !changed; ++i) {
            for (size_t j = i + 1; j < generalBoundarySegments.size(); ++j) {
                BoundarySegment merged;
                if (!tryMergeBoundarySegments(generalBoundarySegments[i], generalBoundarySegments[j], merged)) {
                    continue;
                }
                generalBoundarySegments[i] = merged;
                generalBoundarySegments.erase(generalBoundarySegments.begin() + j);
                changed = true;
                break;
            }
        }
    }
}

static void buildGeneralBoundarySegments() {
    generalBoundarySegments.clear();

    for (size_t i = 0; i < generalRegions.size(); ++i) {
        const ConvexRegion& region = generalRegions[i];
        for (size_t e = 0; e < region.edges.size(); ++e) {
            const Segment& base = region.edges[e];
            double segMinX = min(base.a.x, base.b.x);
            double segMaxX = max(base.a.x, base.b.x);
            double segMinY = min(base.a.y, base.b.y);
            double segMaxY = max(base.a.y, base.b.y);

            vector<double> params;
            params.reserve(16);
            params.push_back(0.0);
            params.push_back(1.0);

            for (size_t j = 0; j < generalRegions.size(); ++j) {
                if (i == j) {
                    continue;
                }
                const ConvexRegion& otherRegion = generalRegions[j];
                if (segMaxX < otherRegion.minx - QUERY_EPS || segMinX > otherRegion.maxx + QUERY_EPS ||
                    segMaxY < otherRegion.miny - QUERY_EPS || segMinY > otherRegion.maxy + QUERY_EPS) {
                    continue;
                }
                for (size_t oe = 0; oe < otherRegion.edges.size(); ++oe) {
                    collectEdgeSplitParams(base, otherRegion.edges[oe], params);
                }
            }

            sort(params.begin(), params.end());
            vector<double> uniq;
            uniq.reserve(params.size());
            for (size_t k = 0; k < params.size(); ++k) {
                if (uniq.empty() || fabs(params[k] - uniq.back()) > 1e-10) {
                    uniq.push_back(params[k]);
                }
            }

            Vec dir = base.b - base.a;
            Vec out(region.nx[e], region.ny[e]);
            for (size_t k = 0; k + 1 < uniq.size(); ++k) {
                double l = uniq[k];
                double r = uniq[k + 1];
                if (r - l <= 1e-10) {
                    continue;
                }

                Vec mid = base.a + dir * (0.5 * (l + r));
                Vec probe = mid + out * OUT_STEP;
                if (pointInsideBucketedGeneralRegionStrict(probe.x, probe.y)) {
                    continue;
                }

                Vec s = base.a + dir * l;
                Vec t = base.a + dir * r;
                if (norm2(t - s) <= GEOM_EPS) {
                    continue;
                }

                generalBoundarySegments.push_back(makeBoundarySegment(s, t));
            }
        }
    }
}

static void preprocessConvexFastPath() {
    vector<Vec> negB(polyB.size());
    for (size_t i = 0; i < polyB.size(); ++i) {
        negB[i] = Vec(-polyB[i].x, -polyB[i].y);
    }

    vector<Vec> nfp = minkowskiSum(polyA, negB);
    int n = static_cast<int>(nfp.size());

    fastNormalX.resize(n);
    fastNormalY.resize(n);
    fastOffsetC.resize(n);

    fastMinX = fastMaxX = nfp[0].x;
    fastMinY = fastMaxY = nfp[0].y;

    for (int i = 0; i < n; ++i) {
        const Vec& a = nfp[i];
        const Vec& b = nfp[(i + 1) % n];
        Vec edge = b - a;
        double len = length(edge);
        double ox = edge.y / len;
        double oy = -edge.x / len;
        fastNormalX[i] = ox;
        fastNormalY[i] = oy;
        fastOffsetC[i] = ox * a.x + oy * a.y;

        if (a.x < fastMinX) fastMinX = a.x;
        if (a.x > fastMaxX) fastMaxX = a.x;
        if (a.y < fastMinY) fastMinY = a.y;
        if (a.y > fastMaxY) fastMaxY = a.y;
    }
}

static vector<vector<Vec> > buildConvexParts(const vector<Vec>& polyInput) {
    vector<Vec> poly = polyInput;
    normalizeCCW(poly);
    if (isConvexPolygon(poly)) {
        return vector<vector<Vec> >(1, poly);
    }
    return mergeConvexParts(triangulateSimplePolygon(poly));
}

static void preprocessGeneralPath() {
    vector<vector<Vec> > partsA = buildConvexParts(polyA);
    vector<vector<Vec> > partsB = buildConvexParts(polyB);

    generalRegions.clear();
    for (size_t ia = 0; ia < partsA.size(); ++ia) {
        for (size_t ib = 0; ib < partsB.size(); ++ib) {
            vector<Vec> negPartB(partsB[ib].size());
            for (size_t i = 0; i < partsB[ib].size(); ++i) {
                negPartB[i] = Vec(-partsB[ib][i].x, -partsB[ib][i].y);
            }
            vector<Vec> sumPoly = minkowskiSum(partsA[ia], negPartB);
            if (sumPoly.size() < 3 || fabs(signedArea2(sumPoly)) <= GEOM_EPS) {
                continue;
            }
            generalRegions.push_back(buildConvexRegion(sumPoly));
        }
    }

    vector<char> removed(generalRegions.size(), 0);
    for (size_t i = 0; i < generalRegions.size(); ++i) {
        if (removed[i]) continue;
        for (size_t j = 0; j < generalRegions.size(); ++j) {
            if (i == j || removed[j]) continue;
            if (generalRegions[i].minx < generalRegions[j].minx - QUERY_EPS ||
                generalRegions[i].maxx > generalRegions[j].maxx + QUERY_EPS ||
                generalRegions[i].miny < generalRegions[j].miny - QUERY_EPS ||
                generalRegions[i].maxy > generalRegions[j].maxy + QUERY_EPS) {
                continue;
            }
            bool inside = true;
            for (size_t k = 0; k < generalRegions[i].pt.size(); ++k) {
                if (!pointInConvexInclusive(generalRegions[j], generalRegions[i].pt[k].x, generalRegions[i].pt[k].y)) {
                    inside = false;
                    break;
                }
            }
            if (inside) {
                removed[i] = 1;
                break;
            }
        }
    }

    vector<ConvexRegion> filtered;
    filtered.reserve(generalRegions.size());
    generalMinX = generalMinY = INF;
    generalMaxX = generalMaxY = -INF;
    for (size_t i = 0; i < generalRegions.size(); ++i) {
        if (removed[i]) continue;
        filtered.push_back(generalRegions[i]);
        if (generalRegions[i].minx < generalMinX) generalMinX = generalRegions[i].minx;
        if (generalRegions[i].maxx > generalMaxX) generalMaxX = generalRegions[i].maxx;
        if (generalRegions[i].miny < generalMinY) generalMinY = generalRegions[i].miny;
        if (generalRegions[i].maxy > generalMaxY) generalMaxY = generalRegions[i].maxy;
    }
    generalRegions.swap(filtered);

    if (generalRegions.empty()) {
        generalMinX = generalMaxX = generalMinY = generalMaxY = 0.0;
    }

    buildGeneralRegionBuckets();
    buildGeneralBoundarySegments();
    mergeGeneralBoundarySegments();
}

static inline bool outsideGeneralBBox(double x, double y) {
    return x < generalMinX - QUERY_EPS || x > generalMaxX + QUERY_EPS ||
           y < generalMinY - QUERY_EPS || y > generalMaxY + QUERY_EPS;
}

static inline void solveConvexFast(double tx, double ty, double& outX, double& outY) {
    if (tx < fastMinX - QUERY_EPS || tx > fastMaxX + QUERY_EPS ||
        ty < fastMinY - QUERY_EPS || ty > fastMaxY + QUERY_EPS) {
        outX = 0.0;
        outY = 0.0;
        return;
    }

    double best = INF;
    double bestX = 0.0, bestY = 0.0;
    for (size_t i = 0; i < fastNormalX.size(); ++i) {
        double slack = fastOffsetC[i] - (fastNormalX[i] * tx + fastNormalY[i] * ty);
        if (slack <= QUERY_EPS) {
            outX = 0.0;
            outY = 0.0;
            return;
        }
        if (slack < best) {
            best = slack;
            bestX = fastNormalX[i] * slack;
            bestY = fastNormalY[i] * slack;
        }
    }

    if (fabs(bestX) < PRINT_EPS) bestX = 0.0;
    if (fabs(bestY) < PRINT_EPS) bestY = 0.0;
    outX = bestX;
    outY = bestY;
}

static inline bool pointInsideAnyGeneralRegionStrict(double x, double y) {
    if (outsideGeneralBBox(x, y)) {
        return false;
    }
    if (!generalUseBucketsForQuery) {
        return pointInsideAnyRegionStrict(generalRegions, x, y);
    }
    return pointInsideBucketedGeneralRegionStrict(x, y);
}

static inline void solveGeneral(double tx, double ty, double& outX, double& outY) {
    if (generalRegions.empty() || generalBoundarySegments.empty() || outsideGeneralBBox(tx, ty)) {
        outX = 0.0;
        outY = 0.0;
        return;
    }

    if (!pointInsideAnyGeneralRegionStrict(tx, ty)) {
        outX = 0.0;
        outY = 0.0;
        return;
    }

    Vec p(tx, ty);
    double bestDist2 = INF;
    Vec bestVec(0.0, 0.0);

    for (size_t i = 0; i < generalBoundarySegments.size(); ++i) {
        if (boundaryBBoxDist2(generalBoundarySegments[i], tx, ty) >= bestDist2 + GEOM_EPS) {
            continue;
        }
        Vec q = closestPointOnSegment(p, generalBoundarySegments[i].seg);
        Vec diff = q - p;
        double dist2 = norm2(diff);
        if (dist2 + GEOM_EPS >= bestDist2) {
            continue;
        }
        bestDist2 = dist2;
        bestVec = diff;
    }

    if (bestDist2 >= INF / 2.0) {
        outX = 0.0;
        outY = 0.0;
        return;
    }

    if (fabs(bestVec.x) < PRINT_EPS) bestVec.x = 0.0;
    if (fabs(bestVec.y) < PRINT_EPS) bestVec.y = 0.0;
    outX = bestVec.x;
    outY = bestVec.y;
}

static inline void applyOutputMargin(double& x, double& y) {
    double len = sqrt(x * x + y * y);
    if (len <= GEOM_EPS) {
        return;
    }
    x += x / len * OUTPUT_MARGIN;
    y += y / len * OUTPUT_MARGIN;
}

static inline void solveQuery(double tx, double ty, double& outX, double& outY) {
    if (useConvexFastPath) {
        solveConvexFast(tx, ty, outX, outY);
    } else {
        solveGeneral(tx, ty, outX, outY);
    }
    applyOutputMargin(outX, outY);
    if (fabs(outX) < PRINT_EPS) outX = 0.0;
    if (fabs(outY) < PRINT_EPS) outY = 0.0;
}

int main() {
    int n1, n2;
    if (scanf("%d%d", &n1, &n2) != 2) {
        return 0;
    }

    polyA.resize(n1);
    polyB.resize(n2);
    for (int i = 0; i < n1; ++i) {
        if (scanf("%lf%lf", &polyA[i].x, &polyA[i].y) != 2) {
            return 0;
        }
    }
    for (int i = 0; i < n2; ++i) {
        if (scanf("%lf%lf", &polyB[i].x, &polyB[i].y) != 2) {
            return 0;
        }
    }

    normalizeCCW(polyA);
    normalizeCCW(polyB);

    char ok[32];
    if (scanf("%31s", ok) != 1 || strcmp(ok, "OK") != 0) {
        return 0;
    }

    useConvexFastPath = isConvexPolygon(polyA) && isConvexPolygon(polyB);
    if (useConvexFastPath) {
        preprocessConvexFastPath();
    } else {
        preprocessGeneralPath();
    }

    puts("OK");
    fflush(stdout);

    int m;
    if (scanf("%d", &m) != 1) {
        return 0;
    }

    vector<Vec> queries(m);
    for (int i = 0; i < m; ++i) {
        if (scanf("%lf%lf", &queries[i].x, &queries[i].y) != 2) {
            return 0;
        }
    }

    if (scanf("%31s", ok) != 1 || strcmp(ok, "OK") != 0) {
        return 0;
    }

    for (int i = 0; i < m; ++i) {
        double ansX, ansY;
        solveQuery(queries[i].x, queries[i].y, ansX, ansY);
        printf("%.5f %.5f\n", ansX, ansY);
        fflush(stdout);
    }
    puts("OK");
    fflush(stdout);
    return 0;
}
