#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <vector>
#include <memory.h>
#include <string>
using namespace std;
const double EPS = 1e-8;
const double PI = 3.141592653;

struct Point {
    double x, y;
    Point(double x_ = 0.0, double y_ = 0.0) :x(x_), y(y_) {}
    bool operator==(const Point& pt) const { return x == pt.x && y == pt.y; }
    bool operator!=(const Point& pt) const { return x != pt.x || y != pt.y; }
    bool operator<(const Point& pt) const { return x != pt.x ? x < pt.x : y < pt.y; }
    Point operator+(const Point& pt) const { return Point(x + pt.x, y + pt.y); }
    Point operator-(const Point& pt) const { return Point(x - pt.x, y - pt.y); }
    Point operator*(double pt) const { return Point(x * pt, y * pt); }
    double norm() const { return hypot(x, y); }
    Point normalize() const { return Point(x / norm(), y / norm()); }
    double dot(const Point& pt) const { return x * pt.x + y * pt.y; }
    double cross(const Point& pt) const { return x * pt.y - pt.x * y; }
};

class Vector2 {
public:
    Vector2() {}
    int ccw(Point p, Point a, Point b);
    double ccw2(Point p, Point a, Point b);
    double pointToPoint(Point a, Point b);
    double intervalAngle(Point a, Point b);
    double area(const vector<Point>& poly);
    bool lineIntersection(Point a, Point b, Point c, Point d, Point& x);
    bool paralleSegments(Point a, Point b, Point c, Point d, Point& p);
    bool inBoundingRectangle(Point p, Point a, Point b);
    pair<bool, Point> segmentIntersection(Point a, Point b, Point c, Point d);
    int segmentIntersects(Point a, Point b, Point c, Point d);
    vector<Point> convexHull(vector<Point>& v);
    vector<Point> pointInConvexPoly(vector<Point>& hull, vector<Point>& pt);
    double areaNonOverlap(vector<Point>& poly1, vector<Point>& poly2);
};

int Vector2::ccw(Point p, Point a, Point b) {
    double ret = (a - p).cross(b - p);
    if (fabs(ret) < EPS) return 0;
    else if (ret < 0) return -1;
    else return 1;
}

double Vector2::ccw2(Point p, Point a, Point b) {
    double ret = (a - p).cross(b - p);
    return ret;
}

double Vector2::pointToPoint(Point a, Point b) {
    Point d = b - a;
    return sqrt(d.x * d.x + d.y * d.y);
}

double Vector2::intervalAngle(Point a, Point b) {
    double t = a.norm() * b.norm();
    return b.dot(a) / t;
}

double Vector2::area(const vector<Point>& p) {
    double ret = 0.0;
    if (p.size() < 3) return 0.0;
    for (int i = 0; i < p.size(); i++) {
        int j = (i + 1) % p.size();
        ret += p[i].cross(p[j]);
    }
    return fabs(ret) / 2.0;
}

bool Vector2::lineIntersection(Point a, Point b, Point c, Point d, Point& x) {
    double det = (b - a).cross(d - c);
    if (fabs(det) < EPS) return false;
    x = a + (b - a) * ((c - a).cross(d - c) / det);
    return true;
}

bool Vector2::paralleSegments(Point a, Point b, Point c, Point d, Point& p) {
    if (b < a) swap(a, b);
    if (d < c) swap(c, d);
    if (ccw(a, b, c) != 0 || b < c || d < a) return false;
    if (a < c) p = c;
    else p = a;
    return true;
}

bool Vector2::inBoundingRectangle(Point p, Point a, Point b) {
    if (b < a) swap(a, b);
    return p == a || p == b || (a < p && p < b);
}

pair<bool, Point> Vector2::segmentIntersection(Point a, Point b, Point c, Point d) {
    pair<bool, Point> res;
    if (!lineIntersection(a, b, c, d, res.second)) return res;
    res.first = inBoundingRectangle(res.second, a, b) && inBoundingRectangle(res.second, c, d);
    return res;
}

int Vector2::segmentIntersects(Point a, Point b, Point c, Point d) {
    if (ccw(a, b, c) == 0 && ccw(a, b, d) == 0) {
        if (b < a) swap(a, b);
        if (d < c) swap(c, d);
        if (b < c || d < a)  return 0;
        else if (b == c || d == a)
            return 1;
        else return 2;
    }
    else {
        int ab = ccw(a, b, c) * ccw(a, b, d);
        int cd = ccw(c, d, a) * ccw(c, d, b);
        return ab <= 0 && cd <= 0;
    }
}

vector<Point> Vector2::convexHull(vector<Point>& v) {
    if (v.size() < 3) return v;
    swap(v[0], *min_element(v.begin(), v.end()));
    sort(++v.begin(), v.end(), [&](Point a, Point b) -> bool {
        int ret = ccw(v[0], a, b);
        return ret > 0 || (ret == 0 && a < b);
        });

    vector<Point> hull;
    hull.emplace_back(v[0]);
    hull.emplace_back(v[1]);
    for (int i = 2; i < v.size(); i++) {
        while (hull.size() > 1 && ccw(v[i], hull[hull.size() - 2], hull[hull.size() - 1]) <= 0) {
            hull.pop_back();
        }
        hull.emplace_back(v[i]);
    }
    return hull;
}

vector<Point> Vector2::pointInConvexPoly(vector<Point>& hull, vector<Point>& pt) {
    int s = hull.size() - 1;
    int mid = max_element(hull.begin(), hull.end()) - hull.begin();

    vector<Point> lowerHull, upperHull;
    upperHull.emplace_back(hull[0]);

    vector<Point> res;
    for (int i = 0; i <= mid; i++) lowerHull.emplace_back(hull[i]);
    for (int i = s; i >= mid; i--) upperHull.emplace_back(hull[i]);
    for (int i = 0; i < pt.size(); i++) {
        int a = lower_bound(lowerHull.begin(), lowerHull.end(), pt[i]) - lowerHull.begin();
        int b = lower_bound(upperHull.begin(), upperHull.end(), pt[i]) - upperHull.begin();

        if (a == 0) {
            if (pt[i] != lowerHull[0]) continue;
        }
        else if (a == lowerHull.size()) continue;
        else {
            if (ccw(pt[i], lowerHull[a - 1], lowerHull[a]) < 0) continue;
        }
        if (b == 0) {
            if (pt[i] != upperHull[0]) continue;
        }
        else if (b == upperHull.size()) continue;
        else {
            if (ccw(pt[i], upperHull[b], upperHull[b - 1]) < 0) continue;
        }

        res.emplace_back(pt[i]);
    }
    return res;
}

double Vector2::areaNonOverlap(vector<Point>& poly1, vector<Point>& poly2) {
    vector<Point> hull1, hull2;
    vector<Point> ans, res, tmp;
    hull1 = convexHull(poly1);
    hull2 = convexHull(poly2);

    res = pointInConvexPoly(hull1, poly2);
    tmp = pointInConvexPoly(hull2, poly1);
    res.insert(res.end(), tmp.begin(), tmp.end());

    for (int i = 0; i < hull1.size(); i++) {
        for (int j = 0; j < hull2.size(); j++) {
            pair<bool, Point> x = segmentIntersection(hull1[i], hull1[(i + 1) % hull1.size()], hull2[j], hull2[(j + 1) % hull2.size()]);
            if (x.first) {
                res.emplace_back(x.second);
            }
        }
    }

    ans = convexHull(res);
    double area1 = area(hull1);
    double area2 = area(hull2);
    double overlapArea = area(ans);

    if(fabs(area1)<EPS && fabs(area2)<EPS) return 1.0;
    else if(fabs(area1)<EPS || fabs(area2)<EPS) return 0.0;
    else return overlapArea*2.0/(area1+area2);
}

class Landmark {
public:
    vector<Point> all;
    vector<Point> rJaw, lJaw;
    vector<Point> rNose, lNose; 
    vector<Point> rEye, lEye;
    vector<Point> rMouth, lMouth;
    vector<Point> midPoint;
    vector<Point> rFace, lFace;

    void setMarks() {
        int i = 0;
        for (i = 0; i <= 8; i++) rJaw.emplace_back(all[i]);
        for (i = 16; i >= 8; i--) lJaw.emplace_back(all[i]);
        for (i = 27; i <= 30; i++) rNose.emplace_back(all[i]);
        for (i = 33; i >= 31; i--) rNose.emplace_back(all[i]);
        for (i = 27; i <= 30; i++) lNose.emplace_back(all[i]);
        for (i = 33; i <= 35; i++) lNose.emplace_back(all[i]);
        for (i = 36; i <= 41; i++) rEye.emplace_back(all[i]);
        for (i = 45; i >= 42; i--) lEye.emplace_back(all[i]);
        for (i = 47; i >= 46; i--) lEye.emplace_back(all[i]);
        for (i = 51; i >= 48; i--) rMouth.emplace_back(all[i]);
        for (i = 59; i >= 57; i--) rMouth.emplace_back(all[i]);
        for (i = 51; i <= 57; i++) lMouth.emplace_back(all[i]);

        rFace.insert(rFace.end(), rEye.begin(), rEye.end());
        lFace.insert(lFace.end(), lEye.begin(), lEye.end());
        rFace.insert(rFace.end(), rNose.begin(), rNose.end());
        lFace.insert(lFace.end(), lNose.begin(), lNose.end());
        rFace.insert(rFace.end(), rMouth.begin(), rMouth.end());
        lFace.insert(lFace.end(), lMouth.begin(), lMouth.end());
        rFace.insert(rFace.end(), rJaw.begin(), rJaw.end());
        lFace.insert(lFace.end(), lJaw.begin(), lJaw.end());
        
        Point tmp;
        for(int i=0; i<rFace.size(); i++){
            tmp = rFace[i] + lFace[i]; 
            tmp.x/=2.0; tmp.y/=2.0;
            midPoint.push_back(tmp);
        }
    }
};

class Figure {
public:
    vector<double> diff;
    double mean=0.0;
    double median=0.0;
    double stdDev=0.0;
    double variance=0.0;
    double skewness=0.0;
    double kurtosis=0.0;

    Figure(vector<double>& diff) { this->diff = diff; }
    void calValues() {
        sort(diff.begin(), diff.end());
        calMeanAndMid();
        calStdDev();
        calSkewAndKur();
    }
    void calMeanAndMid() {
        double sum = 0.0;
        for (int i = 0; i < diff.size(); i++) sum += (double)diff[i];
        mean = sum / (double)diff.size();
        median = diff[diff.size() / 2];
    }
    void calStdDev() {
        double sum = 0.0;
        for (int i = 0; i < diff.size(); i++)
            sum += (diff[i] - mean) * (diff[i] - mean);
        variance = sum / (double)(diff.size() - 1);
        stdDev = sqrt(variance);
    }
    void calSkewAndKur() {
        skewness = 0.0;
        kurtosis = 0.0;
        double sum1 = 0.0, sum2 = 0.0;
        for (int i = 0; i < diff.size(); i++) {
            sum2 = diff[i] - mean;
            sum1 = sum2 * sum2 * sum2;
            sum2 = sum1 * sum2;
            skewness += sum1;
            kurtosis += sum2;
        }
        double std3 = stdDev * stdDev * stdDev;
        if (fabs(std3) < EPS || fabs(std3 * stdDev) < EPS) {
            skewness = 0.0;
            kurtosis = 0.0;
            return;
        }
        skewness = skewness / ((double)diff.size() * std3);
        kurtosis = kurtosis / ((double)diff.size() * std3 * stdDev);
        kurtosis -= 3.0;
    }
};

class Feature {
public:
    Landmark* lm;
    Vector2* v2d;
    Figure* fg;
    vector<double> diff;
    Feature(Landmark* lm, Vector2* v2d) {
        this->lm = lm;
        this->v2d = v2d;
    }
    virtual void calDiff() = 0;
    void calFig() {
        fg = new Figure(this->diff);
        fg->calValues();
    }

};

class Line : public Feature {
public:
    Line(Landmark* lm, Vector2* v2d) : Feature(lm, v2d) {}

    virtual void calDiff() {
        double angle=0.0;
        for (int i = 0; i < lm->midPoint.size(); i++) {
            for (int j = i + 4; j < lm->midPoint.size(); j+=2) {
                for (int a = 0; a < lm->lFace.size(); a++) {
                    if (lm->lFace[a] == lm->rFace[a]) continue;
                    if (lm->midPoint[j] == lm->midPoint[i]) continue;
                    angle = v2d->intervalAngle((lm->lFace[a]) - (lm->rFace[a]), (lm->midPoint[j]) - (lm->midPoint[i]));
                    if(fabs(angle)<EPS) {this->diff.push_back(0.0); continue;}
                    if(angle < 0) angle = fabs(angle);
                    this->diff.push_back(fabs(angle));
                }
            }
        }

    }

};

class Length : public Feature {
public:
    Length(Landmark* lm, Vector2* v2d) : Feature(lm, v2d) {}
    virtual void calDiff() {
        double rLen=0.0, lLen=0.0;
        for (int a = 0; a < lm->lFace.size(); a++) {
            for (int b = a + 1; b < lm->lFace.size(); b++) {
                if (lm->rFace[a] == lm->lFace[a] && lm->rFace[b] == lm->lFace[b]) continue;
                rLen = v2d->pointToPoint(lm->rFace[a], lm->rFace[b]);
                lLen = v2d->pointToPoint(lm->lFace[a], lm->lFace[b]);
                if(fabs(rLen)<EPS && fabs(lLen)<EPS) {this->diff.push_back(1.0); continue;}
                if(rLen > lLen) swap(rLen, lLen);
                this->diff.push_back(rLen/lLen);
            }
        }

    }

};

class Angle : public Feature {
public:
    Angle(Landmark* lm, Vector2* v2d) : Feature(lm, v2d) {}
    virtual void calDiff() {
        double rAngle=0.0, lAngle=0.0;
        for (int a = 0; a < lm->lFace.size(); a++) {
            for (int b = a + 1; b < lm->lFace.size(); b++) {
                for (int c = b + 1; c < lm->lFace.size(); c++) {
                    if (lm->rFace[a] == lm->lFace[a] && lm->rFace[b] == lm->lFace[b] && lm->rFace[c] == lm->lFace[c]) continue;
                    if (lm->rFace[a] == lm->rFace[b] || lm->rFace[c] == lm->rFace[b]) continue;
                    if (lm->lFace[a] == lm->lFace[b] || lm->lFace[c] == lm->lFace[b]) continue;
                    rAngle = v2d->intervalAngle(lm->rFace[a] - lm->rFace[b], lm->rFace[c] - lm->rFace[b]);
                    lAngle = v2d->intervalAngle(lm->lFace[a] - lm->lFace[b], lm->lFace[c] - lm->lFace[b]);
                    if(fabs(rAngle)<EPS && fabs(lAngle)<EPS) {this->diff.push_back(1.0); continue;}
                    if(fabs(rAngle)<EPS) rAngle = 0.0;
                    if(fabs(lAngle)<EPS) lAngle = 0.0;
                    if(rAngle<0) rAngle = fabs(rAngle) + 1.0;
                    if(lAngle<0) lAngle = fabs(lAngle) + 1.0;
                    if(rAngle > lAngle) swap(rAngle, lAngle);
                    this->diff.push_back(fabs(rAngle/lAngle));
                }
            }
        }

    }

};

class Area : public Feature {
public:
    Area(Landmark* lm, Vector2* v2d) : Feature(lm, v2d) {}
    virtual void calDiff() {
        vector<Point> rPoly, lPoly;
        double rArea=0.0, lArea=0.0;
        for (int a = 0; a < lm->rFace.size(); a++) {
            for (int b = a + 1; b < lm->rFace.size(); b++) {
                for (int c = b + 1; c < lm->rFace.size(); c++) {
                    for (int d = c + 1; d < lm->rFace.size(); d++) {
                        rPoly.clear();
                        rPoly = { lm->rFace[a], lm->rFace[b], lm->rFace[c], lm->rFace[d] };
                        rPoly = v2d->convexHull(rPoly);
                        lPoly.clear();
                        lPoly = { lm->lFace[a], lm->lFace[b], lm->lFace[c], lm->lFace[d] };
                        lPoly = v2d->convexHull(lPoly);
                        rArea = v2d->area(rPoly);
                        lArea = v2d->area(lPoly);
                        if(fabs(rArea)<EPS && fabs(lArea)<EPS) {diff.push_back(1.0); continue;}
                        if(rArea > lArea) swap(rArea, lArea);
                        diff.push_back(rArea/lArea);
                    }
                }
            }
        }
    }

};

class LocalSymmetry : public Feature {
public:
    vector<double> eyeDiff;
    vector<double> noseDiff;
    vector<double> mouthDiff;
    vector<double> jawDiff;
    vector<double> faceDiff;
    double eyeMean=0.0, noseMean=0.0, mouthMean=0.0, jawMean=0.0, faceMean=0.0;

    LocalSymmetry(Landmark* lm, Vector2* v2d) : Feature(lm, v2d) {}
    virtual void calDiff() {
        vector<Point> cpyPoints;
        for (int a = 0; a < lm->midPoint.size(); a++) {
            for (int b = a + 1; b < lm->midPoint.size(); b++) {
                lm->rFace.clear();

                cpyPoints = shiftSymmetrically(lm->midPoint[a], lm->midPoint[b], lm->rEye);
                eyeDiff.push_back(v2d->areaNonOverlap(cpyPoints, lm->lEye));
                lm->rFace.insert(lm->rFace.end(), cpyPoints.begin(), cpyPoints.end());
                cpyPoints.clear();

                cpyPoints = shiftSymmetrically(lm->midPoint[a], lm->midPoint[b], lm->rNose);
                noseDiff.push_back(v2d->areaNonOverlap(cpyPoints, lm->lNose));
                lm->rFace.insert(lm->rFace.end(), cpyPoints.begin(), cpyPoints.end());
                cpyPoints.clear();

                cpyPoints = shiftSymmetrically(lm->midPoint[a], lm->midPoint[b], lm->rMouth);
                mouthDiff.push_back(v2d->areaNonOverlap(cpyPoints, lm->lMouth));
                lm->rFace.insert(lm->rFace.end(), cpyPoints.begin(), cpyPoints.end());
                cpyPoints.clear();

                cpyPoints = shiftSymmetrically(lm->midPoint[a], lm->midPoint[b], lm->rJaw);
                jawDiff.push_back(v2d->areaNonOverlap(cpyPoints, lm->lJaw));
                lm->rFace.insert(lm->rFace.end(), cpyPoints.begin(), cpyPoints.end());
                cpyPoints.clear();

                faceDiff.push_back(v2d->areaNonOverlap(lm->rFace, lm->lFace));
                lm->rFace.clear();

            }
        }
        this->diff.insert(this->diff.end(), eyeDiff.begin(), eyeDiff.end());
        this->diff.insert(this->diff.end(), noseDiff.begin(), noseDiff.end());
        this->diff.insert(this->diff.end(), mouthDiff.begin(), mouthDiff.end());
        this->diff.insert(this->diff.end(), jawDiff.begin(), jawDiff.end());
        this->diff.insert(this->diff.end(), faceDiff.begin(), faceDiff.end());
        
        for(int i=0; i<faceDiff.size(); i++){
            eyeMean += eyeDiff[i];
            noseMean += noseDiff[i];
            mouthMean += mouthDiff[i];
            jawMean += jawDiff[i];
            faceMean += faceDiff[i];
        }
        double sz = (double)faceDiff.size();
        eyeMean /= sz;
        noseMean /= sz;
        mouthMean /= sz;
        jawMean /= sz;
        faceMean /= sz;
    }

    vector<Point> shiftSymmetrically(Point pt1, Point pt2, vector<Point> &part) {
        vector<Point> res;
        Point mL, tmp;
        mL = pt1-pt2;
        if(mL.x==0) {
            for(int i=0; i<part.size(); i++){
                tmp.x = pt2.x - (part[i].x - pt2.x);
                tmp.y = part[i].y;
                res.push_back(tmp);
            }
            return res;
        }
        double a = mL.y / mL.x;
        if(fabs(a)>10000.0) a = a>0 ? 10000.0 : -10000.0;
        double b = pt2.y + (-a * pt2.x);
        double x, y;
        for (int i = 0; i < part.size(); i++) {
            x = part[i].x;
            y = part[i].y;
            double nx = ((1.0 - a * a) * x + 2.0 * a * y - 2.0 * a * b) / (1.0 + a * a);
            double ny = (2.0 * a * x - (1.0 - a * a) * y + 2.0 * b) / (1.0 + a * a);
            nx = round(nx * 10000.0) / 10000.0;
            ny = round(ny * 10000.0) / 10000.0;
            tmp.x = nx;
            tmp.y = ny;
            res.push_back(tmp);
        }
        return res;
    }

};

int main(int argc, char** args) {
    ios_base::sync_with_stdio(0);
    cin.tie(0);
    
    // Input data from args

    Landmark *lm = new Landmark();
    for(int i=1; i<=68*2; i+=2){
        lm->all.push_back(Point(atoi(args[i]), atoi(args[i+1])));
    }
    lm->setMarks();

    double dd =atof(args[68*2+1]);
    double eos =atof(args[68*2+2]);

    Vector2* v2d = new Vector2;

    Feature* ft[4];
    ft[0] = new Line(lm, v2d);
    ft[1] = new Length(lm, v2d);
    ft[2] = new Angle(lm, v2d);
    ft[3] = new Area(lm, v2d);
    ft[4] = new LocalSymmetry(lm, v2d);
    
    for(int i=0; i<5; i++) {
        ft[i]->calDiff();
        ft[i]->calFig();
    }

    string result="";
    for(int i=0; i<5; i++){
        result += to_string(ft[i]->fg->mean); result += '\n';
        result += to_string(ft[i]->fg->median); result += '\n';
        result += to_string(ft[i]->fg->variance); result += '\n';
        result += to_string(ft[i]->fg->stdDev); result += '\n';
        result += to_string(ft[i]->fg->skewness); result += '\n';
        result += to_string(ft[i]->fg->kurtosis); result += '\n';
    }
    result += to_string(dd); result += '\n';
    result += to_string(eos); result += '\n';

    // making Score start //
    LocalSymmetry* LSObj = dynamic_cast<LocalSymmetry*>(ft[4]);
    double resScore = 0.0, eye=0.0, nose=0.0, mouth=0.0, jaw=0.0, face=0.0;

    eye = LSObj->eyeMean + 0.132;
    nose = LSObj->noseMean + 0.124;
    mouth = LSObj->mouthMean + 0.059;
    jaw = LSObj->jawMean + 0.078;
    face = LSObj->faceMean + 0.030;
    eye = eye > 1.0 ? 100.0 : eye*100;
    nose = nose > 1.0 ? 100.0 : nose*100;
    mouth = mouth > 1.0 ? 100.0 : mouth*100;
    jaw = jaw > 1.0 ? 100.0 : jaw*100;
    face = face > 1.0 ? 100.0 : face*100;

    resScore += eye;
    resScore += nose;
    resScore += mouth;
    resScore += jaw;
    resScore += face;
    resScore /= 5.0;
    // making Score end //
    
    result += to_string(eye); result += '\n';
    result += to_string(nose); result += '\n';
    result += to_string(mouth); result += '\n';
    result += to_string(jaw); result += '\n';
    result += to_string(face); result += '\n';
    result += to_string(resScore); result += '\n';

    ofstream out("C:\\Users\\ChanuiJeon\\Desktop\\landmarkResult\\result_0.txt");
    out<<result;
    out.close();
    
    cout<<result;
}

