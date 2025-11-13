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
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/draw_voronoi_diagram_2.h>
using namespace std;

// Auxiliary constants
const double INF = numeric_limits<double>::max();   // Infinite is defined as the maximum value a double can store
const double EPS = 1e-9;                            // Epsilon constant to stabilize the comparison of two small numbers

struct Point {
    double x, y;
};

// Helper to convert index to neighborhood name (0 -> A, 1 -> B, etc)
char indexToChar(int i) {
    return 'A' + i;
}

// ==========================================
// Prim's Algorithm for MST (Optimal Wiring)
// ==========================================
/*
    This function finds a set of connections of edges that join all the vertices, which represent the cities, with the
    lowest possible cost. N represent the number of cities, and distMatrix the distances between city i and city j.
    
    Complexity: Because of the nested for loops, the complexity of this algorithm is O(N^2), where N represents the number
    of nodes in the graph.
*/
void primMST(int N, const vector<vector<double>>& distMatrix) {
    cout << "1. Optimal wiring (MST):" << endl;
    vector<double> key(N, INF);         // Will store the shortest known distance to connect node i
    vector<int> parent(N, -1);          // Stores the predecessor of a node i so we can reconstruct the path
    vector<bool> inMST(N, false);       // Helps us identify already included nodes

    key[0] = 0;     // We start from the node 0

    double totalWeight = 0;     // Sum of the weight of all edges in the MST

    // We iterate over all the nodes in the graph
    for (int count = 0; count < N; count++) {
        double minKey = INF;
        int u = -1;

        // Search for a node that hasn't been included yet with the lowest connection cost. We will add it to the MST
        for (int v = 0; v < N; v++) {
            // We try to find a node that has not been included in the MST
            if (!inMST[v] && key[v] < minKey) {
                minKey = key[v];
                u = v;
            }
        }

        if (u == -1) break; // There are no more possible connections

        inMST[u] = true;    // The node is marked as included

        // If the node has a parent, we print the edge that join both nodes and add the distance to the total
        if (parent[u] != -1) {
            cout << "(" << indexToChar(parent[u]) << ", " << indexToChar(u) << ")" << endl;     // Convert numbers to letters
            totalWeight += distMatrix[parent[u]][u];
        }

        // For each neighbor of u
        for (int v = 0; v < N; v++) {
            // We try to see if connecting each node v through u is cheaper, so we change the connection
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
// TSP Nearest Neighbor (Mail Route)
// ==========================================
/*
    This function tries to solve the TSP problem by always visiting the nearest node. N represents the number
    of cities (nodes) and distMatrix the distances between city i and city j.

    Complexity:
*/
void nearestNeighbor(int N, const vector<vector<double>>& distMatrix) {
    cout << "2. Shortest route (Nearest Neighbor):" << endl;

    if (N == 0) return;

    vector<bool> visited(N, false);     // Helps us identify already visited cities
    vector<int> path;                   // Stores the order in which the nodes were visited
    int current = 0;                    // The starting node, that represents 'A'
    double totalDist = 0;               // Acumulated distance of the path

    // Start the trip from A
    path.push_back(current);
    visited[current] = true;

    // Find nearest unvisited neighbor iterating over all the nodes of the graph
    for (int i = 0; i < N - 1; ++i) {
        int nextNode = -1;
        double minDist = INF;

        // For all the nodes, we look for their neighbors
        for (int j = 0; j < N; ++j) {
            // Search between all the unvisited cities
            if (!visited[j] && distMatrix[current][j] < minDist) {
                minDist = distMatrix[current][j];       // The nearest one
                nextNode = j;                           // Store the index
            }
        }

        // Add the nearest city to the path
        if (nextNode != -1) {
            path.push_back(nextNode);
            visited[nextNode] = true;
            totalDist += minDist;       // Update the distance traveled
            current = nextNode;         // Move to that city
        }
    }

    // Return to start once all cities had been visited
    totalDist += distMatrix[current][path[0]];
    path.push_back(path[0]);

    for (size_t i = 0; i < path.size(); i++) {
        cout << indexToChar(path[i]) << (i == path.size() - 1 ? "" : " -> ");
    }
    cout << "\nTotal route distance: " << totalDist << " km" << endl;
    cout << "--------------------------------------------------------" << endl;
}

// ==========================================
// Edmonds-Karp for Max Flow
// ==========================================
/*
    This function finds the maximum flow of a network, represented by a graph.

    Complexity of the whole algorithm: O(V*E^2), because each one of the E edges can cause
    an iteration of the BFS, which costs approximately O(E) (in the worst case, there are more edges
    than vertices), O(V) times.
*/

// Auxiliary function to find augmenting paths. The Edmonds-Karp implementation uses breadth-first search
// Complexity: O(V + E)
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

// Edmonds-Karp
void edmondsKarp(int N, vector<vector<double>> capacity) {
    cout << "3. Maximum information flow (A to " << indexToChar(N-1) << "):" << endl;
    int s = 0;          // Source node
    int t = N - 1;      // Sink node
    vector<vector<double>> rGraph = capacity;       // Will be used as the residual graph
    vector<int> parent(N);      // Node's predecessors
    double max_flow = 0;        // Accumulator for the max flow

    // While there exists an augmenting path
    while (bfs(N, rGraph, s, t, parent)) {
        double path_flow = INF;

        // We take the edge with the minimum capacity
        for (int v = t; v != s; v = parent[v]) {
            int u = parent[v];
            path_flow = min(path_flow, rGraph[u][v]);
        }

        // Update the residual graph
        for (int v = t; v != s; v = parent[v]) {
            int u = parent[v];
            rGraph[u][v] -= path_flow;      // Capacity
            rGraph[v][u] += path_flow;      // Remaining capacity
        }
        max_flow += path_flow;      // Update the maximum flow
    }

    cout << "Max Flow Value: " << max_flow << endl;
    cout << "--------------------------------------------------------" << endl;
}

// ==========================================
// Voronoi Regions (Half-plane intersection)
// ==========================================
/*
    These following functions are an approximation to the Voronoi diagram by a polygon cut approach.
    Complexity: The auxiliary functions have constant O(1) time. The Voronoi function has quadratic
    complexity because of the nested for loops + the cost of cutting the polygon, which is linear
    but has to be done N times, transforming it in quadratic. Considering we have N sites, the total
    complexity of this procedure is O(n^3).
*/

// Geometry helpers
struct Line {
    double a, b, c;     // Represent ax + by + c = 0
};

Line getBisector(Point p1, Point p2) {
    // Midpoint between points
    double mx = (p1.x + p2.x) / 2.0;
    double my = (p1.y + p2.y) / 2.0;

    // Vector p1->p2 is normal to bisector (p2.x - p1.x, p2.y - p1.y)
    double a = p2.x - p1.x;
    double b = p2.y - p1.y;

    double c = -a * mx - b * my;    // Adjusts c value

    return {a, b, c};       // Return the constants for the line equation
}

// Checks which side of line (bisector) a point is on
double distToLine(Line l, Point p) {
    return l.a * p.x + l.b * p.y + l.c;
}

// Calculates the intersection point between two lines
Point intersection(Line l1, Line l2) {
    double det = l1.a * l2.b - l2.a * l1.b;

    if (abs(det) < EPS)
        return {NAN, NAN};

    return { (l1.b * l2.c - l2.b * l1.c) / det, (l2.a * l1.c - l1.a * l2.c) / det };
}

// Clip polygon by line, keep the side where a point is
vector<Point> clipPoly(const vector<Point>& poly, Line l, Point site) {
    vector<Point> newPoly;
    // Determine which side of the line the site is on. We want to keep points on the same side as the point
    double siteSide = distToLine(l, site);

    // Traverse the sides of the polygon looking on each segment
    for (size_t i = 0; i < poly.size(); ++i) {
        Point cur = poly[i];
        Point prev = poly[(i - 1 + poly.size()) % poly.size()];

        double curSide = distToLine(l, cur);
        double prevSide = distToLine(l, prev);

        // If current point is on the correct side (same sign as siteSide, or 0). We use EPS to handle float inaccuracies.
        bool curIn = (siteSide >= -EPS && curSide >= -EPS) || (siteSide <= EPS && curSide <= EPS);
        bool prevIn = (siteSide >= -EPS && prevSide >= -EPS) || (siteSide <= EPS && prevSide <= EPS);

        if (curIn) {
            if (!prevIn) {
                // Entering true side, add intersection
                 Line segLine = { -(cur.y - prev.y), cur.x - prev.x, -(-(cur.y - prev.y)*cur.x + (cur.x - prev.x)*cur.y) };
                 newPoly.push_back(intersection(l, segLine));
            }
            newPoly.push_back(cur);
        }
        else if (prevIn) {
            // Leaving true side, add intersection
             Line segLine = { -(cur.y - prev.y), cur.x - prev.x, -(-(cur.y - prev.y)*cur.x + (cur.x - prev.x)*cur.y) };
             newPoly.push_back(intersection(l, segLine));
        }
    }
    return newPoly;     // Returns the line that passes through the polygon to calculate the intersection with the cut line
}

// CGAL types
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> DT;
typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT> AT;
typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT> AP;
typedef CGAL::Voronoi_diagram_2<DT, AT, AP> VD;
typedef K::Point_2 CGAL_Point;

// This function prints in text the Voronoi regions for each site, calculated by the bounding box method.
// The graphical representation is done with CGAL, using the Delaunay triangulation function.
void voronoi(int N, const vector<Point>& sites) {
    cout << "4. Voronoi Regions (nearest exchange areas):" << endl;

    // Text output
    double MIN_C = 0, MAX_C = 1000;     // Defines the "box" that limits the total space
    vector<Point> boundingBox = { {MIN_C, MIN_C}, {MAX_C, MIN_C}, {MAX_C, MAX_C}, {MIN_C, MAX_C} };     // Initial polygon

    cout << fixed << setprecision(1);       // Print numbers with one decimal

    // Calculation of the Voronoi regions
    for (int i = 0; i < N; ++i) {
        vector<Point> region = boundingBox;     // Region is the whole box

        // For each site j, calculate the bisector and cut the polygon
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;
            Line bisector = getBisector(sites[i], sites[j]);
            region = clipPoly(region, bisector, sites[i]);
        }

        // Print results mapping 0->A, 1->B, etc.
        cout << "Region for " << indexToChar(i) << " (" << sites[i].x << "," << sites[i].y << "): [";
        for (size_t k = 0; k < region.size(); ++k) {
            cout << "(" << region[k].x << "," << region[k].y << ")";
            if (k < region.size() - 1) cout << ", ";
        }
        cout << "]" << endl;
    }

    // Graphic representation
    cout << "\n[CGAL] Preparing interactive Voronoi window..." << endl;
    
    // Convert Point to CGAL_Point
    vector<CGAL_Point> cgalPoints;
    for (const auto& p : sites)
        cgalPoints.push_back(CGAL_Point(p.x, p.y));
    
    // Create the Delaunay triangulation and inserts the points
    DT dt;
    dt.insert(cgalPoints.begin(), cgalPoints.end());
    
    // Create de Voronoir diagram based on the triangulation
    VD vd(dt);

    cout << "Opening CGAL Voronoi viewer..." << endl;
    CGAL::draw(vd);
    
    cout << "--------------------------------------------------------" << endl;
}


int main() {
    // Input reading
    ifstream inFile("../../input.txt");

    if (!inFile) {
        cerr << "Error: Could not open input.txt" << endl;
        return 1;
    }

    // Parsing the first line of the file (number of cities)
    int N;
    if (!(inFile >> N)) {
        cerr << "Error: Could not read N from input file." << endl;
        return 1;
    }

    // Matrix that stores the distance between the different N cities
    vector<vector<double>> distMatrix(N, vector<double>(N));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            inFile >> distMatrix[i][j];
        }
    }

    // Matrix that stores the flow of a network
    vector<vector<double>> flowMatrix(N, vector<double>(N));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            inFile >> flowMatrix[i][j];
        }
    }

    vector<Point> sites(N);
    string line;
    getline(inFile, line); 

    for (int i = 0; i < N; i++) {
        while (getline(inFile, line)) {
            if (line.find('(') != string::npos) {
                break;
            }
        }       
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
    primMST(N, distMatrix);
    nearestNeighbor(N, distMatrix);
    edmondsKarp(N, flowMatrix);
    voronoi(N, sites);

    return 0;
}