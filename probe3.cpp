#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <limits>
#include <algorithm>
#include <iomanip>
#include <queue>
#include <fstream>

using namespace std;

const double INF = numeric_limits<double>::max();
const double EPS = 1e-9;

struct Point {
    double x, y;
};

// Helper to convert index to neighborhood name (0 -> A, 1 -> B, etc.)
char indexToChar(int i) {
    return 'A' + i;
}

// ==========================================
// 1. Prim's Algorithm for MST (Optimal Wiring)
// ==========================================
void solveMST(int N, const vector<vector<double>>& distMatrix) {
    cout << "1. Optimal wiring (MST):" << endl;
    vector<double> key(N, INF);
    vector<int> parent(N, -1);
    vector<bool> inMST(N, false);

    key[0] = 0;

    double totalWeight = 0;
    for (int count = 0; count < N; count++) {
        double minKey = INF;
        int u = -1;

        for (int v = 0; v < N; v++) {
            if (!inMST[v] && key[v] < minKey) {
                minKey = key[v];
                u = v;
            }
        }

        if (u == -1) break; 

        inMST[u] = true;
        if (parent[u] != -1) {
            cout << "(" << indexToChar(parent[u]) << ", " << indexToChar(u) << ")" << endl;
            totalWeight += distMatrix[parent[u]][u];
        }

        for (int v = 0; v < N; v++) {
            if (distMatrix[u][v] > 0 && !inMST[v] && distMatrix[u][v] < key[v]) {
                parent[v] = u;
                key[v] = distMatrix[u][v];
            }
        }
    }
    cout << "Total wiring distance: " << totalWeight << " km" << endl;
    cout << "--------------------------------------------------------" << endl;
}

// ==========================================
// 2. TSP Nearest Neighbor (Mail Route)
// ==========================================
void solveTSP(int N, const vector<vector<double>>& distMatrix) {
    cout << "2. Shortest route (TSP - Nearest Neighbor):" << endl;

    if (N == 0) return;

    vector<bool> visited(N, false);
    vector<int> path;
    int current = 0; // Start at 'A'
    double totalDist = 0;

    path.push_back(current);
    visited[current] = true;

    // Find nearest unvisited neighbor N-1 times
    for (int i = 0; i < N - 1; ++i) {
        int nextNode = -1;
        double minDist = INF;

        for (int j = 0; j < N; ++j) {
            if (!visited[j] && distMatrix[current][j] < minDist) {
                minDist = distMatrix[current][j];
                nextNode = j;
            }
        }

        if (nextNode != -1) {
            path.push_back(nextNode);
            visited[nextNode] = true;
            totalDist += minDist;
            current = nextNode;
        }
    }

    // Return to start
    totalDist += distMatrix[current][path[0]];
    path.push_back(path[0]);

    for (size_t i = 0; i < path.size(); i++) {
        cout << indexToChar(path[i]) << (i == path.size() - 1 ? "" : " -> ");
    }
    cout << "\nTotal route distance: " << totalDist << " km" << endl;
    cout << "--------------------------------------------------------" << endl;
}

// ==========================================
// 3. Edmonds-Karp for Max Flow
// ==========================================
bool bfs(int N, const vector<vector<double>>& rGraph, int s, int t, vector<int>& parent) {
    fill(parent.begin(), parent.end(), -1);
    vector<bool> visited(N, false);
    queue<int> q;

    q.push(s);
    visited[s] = true;
    parent[s] = -1;

    while (!q.empty()) {
        int u = q.front();
        q.pop();

        for (int v = 0; v < N; v++) {
            if (!visited[v] && rGraph[u][v] > 0) {
                if (v == t) {
                    parent[v] = u;
                    return true;
                }
                q.push(v);
                parent[v] = u;
                visited[v] = true;
            }
        }
    }
    return false;
}

void solveMaxFlow(int N, vector<vector<double>> capacity) {
    cout << "3. Maximum information flow (A to " << indexToChar(N-1) << "):" << endl;
    int s = 0;
    int t = N - 1;
    vector<vector<double>> rGraph = capacity; // Residual graph
    vector<int> parent(N);
    double max_flow = 0;

    while (bfs(N, rGraph, s, t, parent)) {
        double path_flow = INF;
        for (int v = t; v != s; v = parent[v]) {
            int u = parent[v];
            path_flow = min(path_flow, rGraph[u][v]);
        }

        for (int v = t; v != s; v = parent[v]) {
            int u = parent[v];
            rGraph[u][v] -= path_flow;
            rGraph[v][u] += path_flow;
        }
        max_flow += path_flow;
    }

    cout << "Max Flow Value: " << max_flow << endl;
    cout << "--------------------------------------------------------" << endl;
}

// ==========================================
// 4. Voronoi Regions (Half-plane intersection)
// ==========================================

// Geometry helpers
struct Line {
    double a, b, c; // ax + by + c = 0
};

Line getBisector(Point p1, Point p2) {
    // Midpoint
    double mx = (p1.x + p2.x) / 2.0;
    double my = (p1.y + p2.y) / 2.0;
    // Vector p1->p2 is normal to bisector: (p2.x - p1.x, p2.y - p1.y)
    double a = p2.x - p1.x;
    double b = p2.y - p1.y;
    // c = -ax - by. Plug in midpoint.
    double c = -a * mx - b * my;
    return {a, b, c};
}

// Check which side of line a point is on. >0 means "inside" standard halfplane if oriented correctly.
// For Voronoi, we want the side containing the site.
double distToLine(Line l, Point p) {
    return l.a * p.x + l.b * p.y + l.c;
}

Point intersection(Line l1, Line l2) {
    double det = l1.a * l2.b - l2.a * l1.b;
    if (abs(det) < EPS) return {NAN, NAN}; // Parallel
    return { (l1.b * l2.c - l2.b * l1.c) / det, (l2.a * l1.c - l1.a * l2.c) / det };
}

// Clip polygon by line. Keep the side where 'site' is.
vector<Point> clipPoly(const vector<Point>& poly, Line l, Point site) {
    vector<Point> newPoly;
    // Determine which side of the line the site is on.
    // We want to keep points on the SAME side as the site.
    double siteSide = distToLine(l, site);

    for (size_t i = 0; i < poly.size(); ++i) {
        Point cur = poly[i];
        Point prev = poly[(i - 1 + poly.size()) % poly.size()];

        double curSide = distToLine(l, cur);
        double prevSide = distToLine(l, prev);

        // If current point is on the correct side (same sign as siteSide, or 0)
        // We use EPS to handle float inaccuracies.
        bool curIn = (siteSide >= -EPS && curSide >= -EPS) || (siteSide <= EPS && curSide <= EPS);
        bool prevIn = (siteSide >= -EPS && prevSide >= -EPS) || (siteSide <= EPS && prevSide <= EPS);

        if (curIn) {
            if (!prevIn) {
                // Entering true side, add intersection
                 Line segLine = { -(cur.y - prev.y), cur.x - prev.x, -(-(cur.y - prev.y)*cur.x + (cur.x - prev.x)*cur.y) };
                 newPoly.push_back(intersection(l, segLine));
            }
            newPoly.push_back(cur);
        } else if (prevIn) {
            // Leaving true side, add intersection
             Line segLine = { -(cur.y - prev.y), cur.x - prev.x, -(-(cur.y - prev.y)*cur.x + (cur.x - prev.x)*cur.y) };
             newPoly.push_back(intersection(l, segLine));
        }
    }
    return newPoly;
}

void solveVoronoi(int N, vector<Point> sites) {
    cout << "4. Voronoi Regions (nearest exchange areas):" << endl;

    // Define a bounding box for the city.
    // Based on example input (coords around 100-500), let's use a generous 0-1000 box.
    double MIN_C = 0, MAX_C = 1000;
    vector<Point> boundingBox = { {MIN_C, MIN_C}, {MAX_C, MIN_C}, {MAX_C, MAX_C}, {MIN_C, MAX_C} };

    cout << fixed << setprecision(1);

    for (int i = 0; i < N; ++i) {
        vector<Point> region = boundingBox;
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;
            Line bisector = getBisector(sites[i], sites[j]);
            region = clipPoly(region, bisector, sites[i]);
        }

        cout << "Region for " << indexToChar(i) << " (" << sites[i].x << "," << sites[i].y << "): [";
        for (size_t k = 0; k < region.size(); ++k) {
            cout << "(" << region[k].x << "," << region[k].y << ")";
            if (k < region.size() - 1) cout << ", ";
        }
        cout << "]" << endl;
    }
    cout << "--------------------------------------------------------" << endl;
}

// ==========================================
// Main Driver
// ==========================================
int main() {
    ifstream inFile("input.txt");

    if (!inFile) {
        cerr << "Error: Could not open input.txt" << endl;
        return 1;
    }

    int N;
    if (!(inFile >> N)) {
        cerr << "Error: Could not read N from input file." << endl;
        return 1;
    }

    vector<vector<double>> distMatrix(N, vector<double>(N));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            inFile >> distMatrix[i][j];
        }
    }

    vector<vector<double>> flowMatrix(N, vector<double>(N));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            inFile >> flowMatrix[i][j];
        }
    }

    vector<Point> sites(N);
    string line;
    // consume newline after the last matrix element
    getline(inFile, line); 

    for (int i = 0; i < N; i++) {
        getline(inFile, line);
        // Parse (x,y)
        size_t openP = line.find('(');
        size_t comma = line.find(',');
        size_t closeP = line.find(')');
        if (openP != string::npos && comma != string::npos && closeP != string::npos) {
            sites[i].x = stod(line.substr(openP + 1, comma - openP - 1));
            sites[i].y = stod(line.substr(comma + 1, closeP - comma - 1));
        }
    }

    inFile.close();

    cout << "--- RESULTS ---" << endl;
    solveMST(N, distMatrix);
    solveTSP(N, distMatrix);
    solveMaxFlow(N, flowMatrix);
    solveVoronoi(N, sites);

    return 0;
}