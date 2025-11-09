#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <cctype>

// CGAL Includes
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

// Definiciones de tipos para CGAL
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Delaunay_triangulation_2<K> Delaunay;

// Estructura de tus puntos originales
struct Pt { double x,y; };

// Nombre de archivo de entrada (debes asegurarte que es correcto)
#define INPUT_FILENAME "input.txt"

// Función auxiliar para parsear un par (x, y)
std::pair<double, double> parse_pair(const std::string& s) {
    // Busca números flotantes. Elimina paréntesis y sustituye comas por espacios
    std::string clean_s;
    for (char c : s) {
        if (c == '(' || c == ')') continue;
        if (c == ',') clean_s.push_back(' '); // Tratar ',' como espacio para stringstream
        else clean_s.push_back(c);
    }

    std::stringstream ss(clean_s);
    double x, y;

    // Intenta leer dos dobles.
    if (ss >> x >> y) {
        return {x, y};
    }
    
    throw std::runtime_error("Failed to parse pair in string: " + s);
}


// Función principal para calcular Voronoi con CGAL
void compute_voronoi_cgal(const std::vector<Pt>& sites) {
    if (sites.empty()) return;

    // 1. Convertir estructuras Pt a Point de CGAL
    std::vector<Point> cgal_sites;
    for (const auto& p : sites) {
        cgal_sites.emplace_back(p.x, p.y);
    }

    // 2. Construir la Triangulación de Delaunay
    Delaunay dt;
    dt.insert(cgal_sites.begin(), cgal_sites.end());

    std::cout << "\n4) Voronoi polygons (CGAL - dual of Delaunay):\n";
    std::cout.setf(std::ios::fixed);

    // Utilizamos un iterador para los sitios, ya que CGAL no garantiza el orden de inserción.
    int i = 0;
    // Iterar sobre los vértices de Delaunay (sitios de entrada)
    for (Delaunay::Finite_vertices_iterator vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit, ++i) {
        Point p_site = vit->point();
        
        std::cout << "   Site " << i // Usamos el índice de lectura como identificador
                  << " at (" << std::setprecision(2) << p_site.x() << "," << p_site.y() << ")\n";
        
        // Circulador para recorrer las caras incidentes a este vértice (triángulos de Delaunay)
        Delaunay::Face_circulator fc = dt.incident_faces(vit);
        Delaunay::Face_circulator done = fc;

        std::vector<Point> polygon_vertices;

        if (fc != 0) {
            do {
                if (!dt.is_infinite(fc)) {
                    // El circuncentro (dual de la cara) es un vértice de Voronoi
                    Point voronoi_vertex = dt.dual(fc);
                    polygon_vertices.push_back(voronoi_vertex);
                } 
                // Los polígonos no acotados (caras infinitas) se omiten para simplicidad
                ++fc;
            } while (fc != done);
        }

        std::cout << "     Polygon points (" << polygon_vertices.size() << " finite vertices):\n";
        for (const auto& v : polygon_vertices) {
             std::cout << "      (" << std::setprecision(6) << v.x() << ", " << v.y() << ")\n";
        }
        std::cout << "\n";
    }
}

// -------------------- Lógica de Lectura del Archivo (Integrada) --------------------
int main(){
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    std::ifstream fin(INPUT_FILENAME);
    if(!fin){
        std::cerr << "Cannot open input file '" << INPUT_FILENAME << "'. Place your input in that file and run the program.\n";
        return 1;
    }

    // Leer todas las líneas (incluyendo vacías)
    std::vector<std::string> lines;
    std::string line;
    while(getline(fin,line)){
        std::string t;
        for(char c: line) if(c!='\r' && c!='\n') t.push_back(c);
        lines.push_back(t);
    }
    fin.close();

    // L: vector con solo líneas NO vacías o que contienen contenido no-espacio.
    std::vector<std::string> L;
    for(auto &s: lines){
        bool has_content = false;
        for(char c: s) if(!std::isspace((unsigned char)c)) has_content = true;
        if(has_content) L.push_back(s); 
    }
    if(L.empty()){
        std::cerr << "Input file seems empty.\n";
        return 1;
    }
    
    int N = 0;
    int idx = 0;
    
    // 1. Parsear N (la primera línea con contenido numérico)
    while(idx < (int)L.size()){
        std::string s = L[idx];
        std::string trimmed;
        for(char c: s) if(!std::isspace((unsigned char)c)) trimmed.push_back(c);
        
        bool ok = !trimmed.empty() && ( (trimmed[0]=='-' || trimmed[0]=='+' || std::isdigit(trimmed[0])) );
        if(!ok){ idx++; continue; } 

        try{
            N = std::stoi(trimmed);
            idx++; // Avanza al inicio de la primera matriz
            break;
        } catch(const std::exception&){
            idx++;
            continue;
        }
    }
    
    if (N <= 0) {
        std::cerr << "Failed to parse N or N is invalid.\n";
        return 1;
    }

    // 2. Saltar la matriz de distancia (N líneas)
    int read_dist = 0;
    while(read_dist < N && idx < (int)L.size()){
        // Simulamos la lectura de la línea de la matriz
        read_dist++;
        idx++;
    }
    if (read_dist != N) {
        std::cerr << "Incomplete distance matrix.\n";
        return 1;
    }

    // 3. Saltar la matriz de capacidad (N líneas)
    int read_cap = 0;
    while(read_cap < N && idx < (int)L.size()){
        // Simulamos la lectura de la línea de la matriz
        read_cap++;
        idx++;
    }
    if (read_cap != N) {
        std::cerr << "Incomplete capacity matrix.\n";
        return 1;
    }

    // 4. Leer N coordenadas
    std::vector<Pt> sites;
    int read_coords = 0;
    
    while(read_coords < N && idx < (int)L.size()){
        std::string sline = L[idx++];
        
        // Recortar espacios en blanco
        auto lpos = sline.find_first_not_of(" \t");
        auto rpos = sline.find_last_not_of(" \t");
        if(lpos==std::string::npos) continue;
        std::string s2 = sline.substr(lpos, rpos-lpos+1);
        
        try {
            auto p = parse_pair(s2);
            sites.push_back({p.first, p.second});
            read_coords++;
        } catch(...){
            // Intento de fallback: quitar todos los espacios
            std::string tmp;
            for(char c: s2) if(!std::isspace((unsigned char)c)) tmp.push_back(c);
            try {
                auto p = parse_pair(tmp);
                sites.push_back({p.first, p.second});
                read_coords++;
            } catch(...){
                std::cerr << "Warning: Could not parse coordinate line: " << sline << "\n";
            }
        }
    }
    
    if ((int)sites.size() != N) {
        std::cerr << "Error: Expected " << N << " coordinates but found " << sites.size() << ".\n";
        return 1;
    }

    // Llama al algoritmo de Voronoi de CGAL
    compute_voronoi_cgal(sites);

    return 0;
}