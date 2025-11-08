// Diego Núñez García A01644371

#include <bits/stdc++.h>

using namespace std;

struct Neighborhood {
    int id;     // Identificador
    double x, y; // Coordenadas
};

struct Network {
    int N; // Número de vecindarios
    vector<vector<double>> distanceMatrix; // Matriz de distancias
    vector<vector<double>> capacityMatrix; // Matriz de capacidades
    vector<Neighborhood> neighborhoods;    // Lista de vecindarios
};

Network readNetworkFile(const string& filename) {
    ifstream file(filename); // Abre el archivo
    Network network;
    string line;

    if (!file.is_open()) {
        cerr << "Error: no se pudo abrir el archivo :'v \n";
        return network;
    }

    file >> network.N; // Lee el número de vecindarios

    // lee los valores de la matriz de distancias del archivo
    for (int i = 0; i < network.N; i++) { 

        vector<double> row(network.N); // fila vacía con N columnas

        for (int j = 0; j < network.N; j++) { 
            file >> row[j]; 
        }
        network.distanceMatrix.push_back(row); // guarda la fila en la matriz
    }

    // lee los valores de la matriz de capacidades del archivo
    for (int i = 0; i < network.N; i++) { 

        vector<double> row(network.N); // fila vacía con N columnas

        for (int j = 0; j < network.N; j++) { 
            file >> row[j]; 
        }
        network.capacityMatrix.push_back(row); // guarda la fila en la matriz
    }

    // coordenadas 
    for (int i = 0; i < network.N; i++) {

        string coord;
        file >> coord; 

        // esto pues pa solo ver el numero y los parentesis
        coord.erase(remove(coord.begin(), coord.end(), '('), coord.end());
        coord.erase(remove(coord.begin(), coord.end(), ')'), coord.end());
        replace(coord.begin(), coord.end(), ',', ' ');

        stringstream ss(coord); // leer los 2 numeros dentro del parentesis

        Neighborhood n;
        n.id = i + 1;
        ss >> n.x >> n.y;
        network.neighborhoods.push_back(n);
    }

    return network;
}

int main() {
    string filename = "input.txt";
    Network network = readNetworkFile(filename);

    if (network.N == 0) {
        cerr << "Error :'v\n";
        return 1;
    }

    cout << "Número de vecindarios: " << network.N << "\n\n";

    cout << "Matriz de distancias:\n";
    for (auto& row : network.distanceMatrix) {
        for (auto val : row) cout << val << " ";
        cout << "\n";
    }

cout << "\nMatriz de capacidades:\n";
    for (auto& row : network.capacityMatrix) {
        for (auto val : row) cout << val << " ";
        cout << "\n";
    }   

cout << "\nCoordenadas:\n";
    for (auto& n : network.neighborhoods){
        cout << n.id << " (" << n.x << ", " << n.y << ")\n";
    }
        
    return 0;
}
