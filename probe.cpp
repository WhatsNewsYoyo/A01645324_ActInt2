// network_tools.cpp
// Compile: g++ -O2 -std=c++17 network_tools.cpp -o network_tools
// Run: ./network_tools   (reads input.txt from current folder)

#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>
#include <string>
#include <unordered_map>
#include <sstream>    // Para stringstream, stoi, stod, getline
#include <fstream>    // Para ifstream
#include <iomanip>    // Para setprecision, fixed

using namespace std;

using namespace std;

const string INPUT_FILENAME = "input.txt";
const double EPS = 1e-9;

struct Edge {
    int u, v;
    double w;
    bool operator<(Edge const& o) const { return w < o.w; }
};

struct DSU {
    int n;
    vector<int> p, r;
    DSU(int n=0): n(n), p(n), r(n,0) { for(int i=0;i<n;i++) p[i]=i; }
    int find(int a){ return p[a]==a ? a : p[a]=find(p[a]); }
    bool unite(int a,int b){
        a=find(a); b=find(b);
        if(a==b) return false;
        if(r[a]<r[b]) swap(a,b);
        p[b]=a;
        if(r[a]==r[b]) r[a]++;
        return true;
    }
};

// parse coordinate like "(200,500)" or "200,500"
pair<double,double> parse_pair(const string &s0){
    string s = s0;
    // remove parentheses and spaces
    string t;
    for(char c: s) if(c!='(' && c!=')' && c!=' ') t.push_back(c);
    size_t comma = t.find(',');
    if(comma==string::npos) throw runtime_error("Bad coordinate: "+s0);
    double x = stod(t.substr(0,comma));
    double y = stod(t.substr(comma+1));
    return {x,y};
}

// -------------------- MST (Kruskal) --------------------
vector<pair<int,int>> compute_mst(const vector<vector<double>>& dist, double &totalLen){
    int n = dist.size();
    vector<Edge> edges;
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){
        if(dist[i][j] >= 0) edges.push_back({i,j,dist[i][j]});
    }
    sort(edges.begin(), edges.end());
    DSU dsu(n);
    vector<pair<int,int>> mst;
    totalLen = 0;
    for(auto &e: edges){
        if(dsu.unite(e.u,e.v)){
            mst.emplace_back(e.u,e.v);
            totalLen += e.w;
        }
    }
    return mst;
}

// -------------------- TSP Held-Karp (exact) --------------------
pair<vector<int>, double> tsp_held_karp(const vector<vector<double>>& dist){
    int n = dist.size();
    if(n==0) return {{},0.0};
    if(n==1) return {{0,0},0.0};
    int FULL = 1<<n;
    const double INF = 1e18;
    // dp[mask][i] stored as unordered_map<int,double> per mask (only entries with end node i)
    vector<unordered_map<int,double>> dp(FULL);
    vector<unordered_map<int,int>> parent(FULL);
    dp[1<<0][0] = 0.0;
    for(int mask=0; mask<FULL; ++mask){
        if(!(mask & 1)) continue; // must include start
        for(auto &p : dp[mask]){
            int u = p.first;
            double cost = p.second;
            for(int v=0; v<n; ++v){
                if(mask & (1<<v)) continue;
                int nm = mask | (1<<v);
                double nc = cost + dist[u][v];
                auto it = dp[nm].find(v);
                if(it==dp[nm].end() || nc < it->second - EPS){
                    dp[nm][v] = nc;
                    parent[nm][v] = u;
                }
            }
        }
    }
    int full = FULL-1;
    double best = INF; int last = -1;
    for(auto &p: dp[full]){
        int u = p.first;
        if(u==0) continue;
        double c = p.second + dist[u][0];
        if(c < best - EPS){ best = c; last = u; }
    }
    if(last==-1) return {{}, INF};
    // reconstruct (in reverse)
    vector<int> path;
    int mask = full;
    int cur = last;
    while(mask != (1<<0)){
        path.push_back(cur);
        int prev = parent[mask][cur];
        mask ^= (1<<cur);
        cur = prev;
    }
    path.push_back(0);
    reverse(path.begin(), path.end());
    // add final return to 0
    path.push_back(0);
    return {path, best};
}

// -------------------- Max Flow (Edmonds-Karp) --------------------
double edmonds_karp(vector<vector<double>> cap, int s, int t){
    int n = cap.size();
    double maxflow = 0.0;
    while(true){
        vector<int> parent(n, -1);
        vector<double> flow(n, 0.0);
        queue<int> q;
        q.push(s);
        parent[s] = -2; flow[s] = 1e18;
        bool found = false;
        while(!q.empty() && !found){
            int u = q.front(); q.pop();
            for(int v=0; v<n; ++v){
                if(parent[v]==-1 && cap[u][v] > EPS){
                    parent[v] = u;
                    flow[v] = min(flow[u], cap[u][v]);
                    if(v==t){
                        double add = flow[t];
                        maxflow += add;
                        // update residuals
                        int cur = t;
                        while(cur != s){
                            int prev = parent[cur];
                            cap[prev][cur] -= add;
                            if(cap[prev][cur] < EPS) cap[prev][cur] = 0.0;
                            cap[cur][prev] += add;
                            cur = prev;
                        }
                        found = true;
                        break;
                    }
                    q.push(v);
                }
            }
        }
        if(!found) break;
    }
    return maxflow;
}

// -------------------- Voronoi via half-plane intersection --------------------
struct Pt { double x,y; };
struct HalfPlane { double a,b,c; }; // a*x + b*y <= c

// Sutherland-Hodgman style clipping of polygon by half-plane a*x + b*y <= c
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
            // leaving: add intersection
            double dx = Q.x - P.x, dy = Q.y - P.y;
            double denom = hp.a*dx + hp.b*dy;
            if(fabs(denom) > EPS){
                double t = (hp.c - hp.a*P.x - hp.b*P.y) / denom;
                Pt I{P.x + t*dx, P.y + t*dy};
                out.push_back(I);
            }
        } else if(!inP && inQ){
            // entering: add intersection then Q
            double dx = Q.x - P.x, dy = Q.y - P.y;
            double denom = hp.a*dx + hp.b*dy;
            if(fabs(denom) > EPS){
                double t = (hp.c - hp.a*P.x - hp.b*P.y) / denom;
                Pt I{P.x + t*dx, P.y + t*dy};
                out.push_back(I);
            }
            out.push_back(Q);
        } else {
            // both outside: nothing
        }
    }
    // remove near duplicates
    vector<Pt> clean;
    for(auto &p: out){
        if(clean.empty() || fabs(clean.back().x - p.x) > 1e-8 || fabs(clean.back().y - p.y) > 1e-8)
            clean.push_back(p);
    }
    if(clean.size()>1){
        Pt &f = clean.front(), &l = clean.back();
        if(fabs(f.x - l.x) < 1e-9 && fabs(f.y - l.y) < 1e-9) clean.pop_back();
    }
    return clean;
}

vector<vector<Pt>> compute_voronoi(const vector<Pt>& sites){
    int n = sites.size();
    vector<double> xs(n), ys(n);
    for(int i=0;i<n;i++){ xs[i]=sites[i].x; ys[i]=sites[i].y; }
    double minx = *min_element(xs.begin(), xs.end());
    double maxx = *max_element(xs.begin(), xs.end());
    double miny = *min_element(ys.begin(), ys.end());
    double maxy = *max_element(ys.begin(), ys.end());
    double margin = max(maxx-minx, maxy-miny) + 1000.0;
    double bx0 = minx - margin, by0 = miny - margin;
    double bx1 = maxx + margin, by1 = maxy + margin;
    // initial bounding rectangle CCW
    vector<Pt> bbox = { {bx0, by0}, {bx1, by0}, {bx1, by1}, {bx0, by1} };
    vector<vector<Pt>> cells(n);
    for(int i=0;i<n;i++){
        vector<Pt> poly = bbox;
        for(int j=0;j<n;j++){
            if(i==j) continue;
            // half-plane for points closer (or equal) to i than j:
            // 2*(xj-xi)*x + 2*(yj-yi)*y <= xj^2 + yj^2 - xi^2 - yi^2
            double a = 2.0*(sites[j].x - sites[i].x);
            double b = 2.0*(sites[j].y - sites[i].y);
            double c = sites[j].x*sites[j].x + sites[j].y*sites[j].y - sites[i].x*sites[i].x - sites[i].y*sites[i].y;
            HalfPlane hp{a,b,c};
            poly = clip_polygon_halfplane(poly, hp);
            if(poly.empty()) break;
        }
        cells[i] = poly;
    }
    return cells;
}


// -------------------- Input read --------------------
int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    ifstream fin(INPUT_FILENAME);
    if(!fin){
        cerr << "Cannot open input file '" << INPUT_FILENAME << "'. Place your input in that file and run the program.\n";
        return 1;
    }

    // Read all non-empty lines
    vector<string> lines;
    string line;
    while(getline(fin,line)){
        // trim
        string t;
        for(char c: line) if(c!='\r' && c!='\n') t.push_back(c);
        // keep line even if blank? we will ignore empty later
        lines.push_back(t);
    }
    fin.close();
    // Keep only non-empty trimmed lines
    vector<string> L;
    for(auto &s: lines){
        string t;
        for(char c: s) if(!isspace((unsigned char)c)) t.push_back(c);
        if(!t.empty()) L.push_back(s); // keep original spacing trimmed line
    }
    if(L.empty()){
        cerr << "Input file seems empty.\n";
        return 1;
    }

    // parse N
    int idx = 0;
    // skip completely blank lines at top
    while(idx < (int)L.size()){
        string s = L[idx];
        // try parse as integer by extracting digits
        string trimmed;
        for(char c: s) if(!isspace((unsigned char)c)) trimmed.push_back(c);
        // ensure trimmed is numeric (possibly with + -)
        bool ok = !trimmed.empty() && ( (trimmed[0]=='-' || trimmed[0]=='+' || isdigit(trimmed[0])) );
        if(!ok){ idx++; continue; }
        try{
            int N = stoi(trimmed);
            // advance and parse matrices
            idx++;
            // read distance matrix: next N non-empty lines with N numbers each
            vector<vector<double>> dist(N, vector<double>(N, 0.0));
            int read = 0;
            while(read < N && idx < (int)L.size()){
                string sline = L[idx++];
                // parse numbers in sline
                stringstream ss(sline);
                double v;
                vector<double> row;
                while(ss >> v) row.push_back(v);
                if((int)row.size() == N){
                    dist[read] = row;
                    read++;
                } else {
                    // maybe multiple spaces and strange formatting; try manual parse
                    row.clear();
                    string token;
                    string tmp;
                    for(char c: sline) if(c!='\t') tmp.push_back(c); else tmp.push_back(' ');
                    stringstream ss2(tmp);
                    while(ss2 >> token){
                        try{ row.push_back(stod(token)); } catch(...) {}
                    }
                    if((int)row.size()==N){
                        dist[read] = row;
                        read++;
                    } else {
                        // skip badly formed line
                    }
                }
            }
            // capacity matrix
            vector<vector<double>> cap(N, vector<double>(N,0.0));
            read = 0;
            while(read < N && idx < (int)L.size()){
                string sline = L[idx++];
                stringstream ss(sline);
                double v;
                vector<double> row;
                while(ss >> v) row.push_back(v);
                if((int)row.size() == N){
                    cap[read] = row;
                    read++;
                } else {
                    // try manual
                    row.clear();
                    string token;
                    string tmp;
                    for(char c: sline) if(c!='\t') tmp.push_back(c); else tmp.push_back(' ');
                    stringstream ss2(tmp);
                    while(ss2 >> token){
                        try{ row.push_back(stod(token)); } catch(...) {}
                    }
                    if((int)row.size()==N){
                        cap[read] = row;
                        read++;
                    } else {
                        // skip
                    }
                }
            }
            // coordinates: read N lines containing (x,y)
            vector<pair<double,double>> coords;
            read = 0;
            while(read < N && idx < (int)L.size()){
                string sline = L[idx++];
                // find first '(' or digits; we will extract until newline
                // remove spaces at ends
                string s2 = sline;
                // remove leading/trailing spaces
                auto lpos = s2.find_first_not_of(" \t");
                auto rpos = s2.find_last_not_of(" \t");
                if(lpos==string::npos) continue;
                s2 = s2.substr(lpos, rpos-lpos+1);
                try {
                    auto p = parse_pair(s2);
                    coords.push_back(p);
                    read++;
                } catch(...){
                    // try to remove spaces and parse again
                    string tmp;
                    for(char c: s2) if(!isspace((unsigned char)c)) tmp.push_back(c);
                    try {
                        auto p = parse_pair(tmp);
                        coords.push_back(p);
                        read++;
                    } catch(...){
                        // skip
                    }
                }
            }

            // At this point we have N, dist, cap, coords
            // Convert coords to Pt
            vector<Pt> sites;
            for(auto &pr: coords) sites.push_back({pr.first, pr.second});

            // 1) MST
            double mstLen = 0.0;
            auto mst = compute_mst(dist, mstLen);

            // 2) TSP exact (Held-Karp)
            auto tsp_res = tsp_held_karp(dist);
            vector<int> tsp_path = tsp_res.first;
            double tsp_cost = tsp_res.second;

            // 3) Max flow from 0 to N-1
            double maxflow = edmonds_karp(cap, 0, N-1);

            // 4) Voronoi polygons
            auto vor = compute_voronoi(sites);

            // Outputs
            auto idx_to_letter = [](int i)->string{
                string s;
                if(i<26) s.push_back(char('A'+i));
                else {
                    int x = i;
                    while(x>=0){
                        s.push_back(char('A' + (x%26)));
                        x = x/26 - 1;
                    }
                    reverse(s.begin(), s.end());
                }
                return s;
            };

            cout.setf(std::ios::fixed); cout<<setprecision(2);
            cout << "1) Optimal wiring (MST) as list of arcs (A,B):\n";
            for(auto &e: mst){
                cout << "(" << idx_to_letter(e.first) << "," << idx_to_letter(e.second) << ")\n";
            }
            cout << "   Total MST length: " << mstLen << " km\n\n";

            cout << "2) Shortest route for mail (TSP exact), start and end at A:\n";
            if(tsp_path.empty()){
                cout << "   No TSP solution found or N too small.\n\n";
            } else {
                cout << "   Route: ";
                for(size_t i=0;i<tsp_path.size();++i){
                    cout << idx_to_letter(tsp_path[i]);
                    if(i+1<tsp_path.size()) cout << " -> ";
                }
                cout << "\n";
                cout << "   Total distance: " << tsp_cost << " km\n\n";
            }

            cout << "3) Maximum information flow from initial node (A) to final node (" << idx_to_letter(N-1) << "):\n";
            cout << "   Max flow value: " << maxflow << "\n\n";

            cout << "4) Voronoi polygons (each polygon is a list of (x,y) points):\n";
            int nSites = (int)vor.size();
            cout.setf(std::ios::fixed);
            for (int i = 0; i < nSites; ++i) {
                cout << "   Site " << idx_to_letter(i)
                     << " at (" << sites[i].x << "," << sites[i].y << ") -> polygon points ("
                     << vor[i].size() << " vertices):\n";
                for (auto &p : vor[i]) {
                    cout << "      (" << setprecision(6) << p.x << ", " << p.y << ")\n";
                }
                cout << "\n";
            }
            cout << setprecision(2); // rest of output uses 2 decimals if quieres
            
            // finished main processing
            return 0;

        } catch(exception &e){
            // couldn't parse N in this line, try next line
            idx++;
            continue;
        }
    }

    cerr << "Failed to parse input. Make sure file format matches:\n"
         << "N\n"
         << "N lines with N numbers (distance matrix)\n"
         << "N lines with N numbers (capacity matrix)\n"
         << "N lines with coordinates like (200,500)\n";
    return 1;
}
