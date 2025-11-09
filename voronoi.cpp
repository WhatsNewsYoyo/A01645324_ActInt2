#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/draw_polygon_2.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <CGAL/Polygon_2.h>
using namespace std;

#define EPS 1e-9

struct Pt { double x, y; };
struct HalfPlane { double a,b,c; };

vector<Pt> clip_polygon_halfplane(const vector<Pt>& poly, const HalfPlane &hp){
    vector<Pt> out;
    int m = poly.size();
    if(m==0) return out;
    for(int i=0;i<m;i++){
        Pt P = poly[i];
        Pt Q = poly[(i+1)%m];
        double fP = hp.a*P.x + hp.b*P.y - hp.c;
        double fQ = hp.a*Q.x + hp.b*Q.y - hp.c;
        bool inP = fP <= EPS;
        bool inQ = fQ <= EPS;
        if(inP && inQ){
            out.push_back(Q);
        } else if(inP && !inQ){
            double dx = Q.x - P.x, dy = Q.y - P.y;
            double denom = hp.a*dx + hp.b*dy;
            if(fabs(denom) > EPS){
                double t = (hp.c - hp.a*P.x - hp.b*P.y) / denom;
                out.push_back({P.x + t*dx, P.y + t*dy});
            }
        } else if(!inP && inQ){
            double dx = Q.x - P.x, dy = Q.y - P.y;
            double denom = hp.a*dx + hp.b*dy;
            if(fabs(denom) > EPS){
                double t = (hp.c - hp.a*P.x - hp.b*P.y) / denom;
                out.push_back({P.x + t*dx, P.y + t*dy});
            }
            out.push_back(Q);
        }
    }
    return out;
}

vector<vector<Pt>> compute_voronoi(const vector<Pt>& sites){
    int n = sites.size();
    double minx=1e9,maxx=-1e9,miny=1e9,maxy=-1e9;
    for(auto&s:sites){
        minx=min(minx,s.x); maxx=max(maxx,s.x);
        miny=min(miny,s.y); maxy=max(maxy,s.y);
    }
    double margin=max(maxx-minx,maxy-miny)+100;
    double bx0=minx-margin, by0=miny-margin;
    double bx1=maxx+margin, by1=maxy+margin;
    vector<Pt> bbox={{bx0,by0},{bx1,by0},{bx1,by1},{bx0,by1}};
    vector<vector<Pt>> cells(n);
    for(int i=0;i<n;i++){
        vector<Pt> poly=bbox;
        for(int j=0;j<n;j++){
            if(i==j) continue;
            double a=2*(sites[j].x - sites[i].x);
            double b=2*(sites[j].y - sites[i].y);
            double c=sites[j].x*sites[j].x + sites[j].y*sites[j].y - sites[i].x*sites[i].x - sites[i].y*sites[i].y;
            poly=clip_polygon_halfplane(poly,{a,b,c});
            if(poly.empty()) break;
        }
        cells[i]=poly;
    }
    return cells;
}

int main(){
    using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Point_2 = Kernel::Point_2;
    using Polygon_2 = CGAL::Polygon_2<Kernel>;

    vector<Pt> sites = {{0,0},{100,50},{50,100},{150,80}};
    auto cells = compute_voronoi(sites);

    for(int i=0;i<sites.size();i++){
        vector<Point_2> pts;
        for(auto&p:cells[i]) pts.push_back(Point_2(p.x,p.y));
        Polygon_2 poly(pts.begin(), pts.end());
        cout << "Celda " << i << ": " << poly.size() << " vértices\n";
        CGAL::draw(poly); // Dibuja cada celda
    }

    return 0;
}
