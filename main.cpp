#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>
#include <chrono>
#include <algorithm>
#include <iomanip>

// Estructura para almacenar las estadísticas globales (Tarea B)
struct Estadisticas {
    long long total_A = 0;
    long long total_C = 0;
    long long total_G = 0;
    long long total_T = 0;
    double porcentaje_gc = 0.0;
    double longitud_media = 0.0;
};

// ============================================================================
// TAREA A: Lectura desde archivo FASTA o generación aleatoria reproducible
// ============================================================================
std::vector<std::string> cargar_o_generar_secuencias(int N, int L, const std::string& ruta_archivo) {
    std::vector<std::string> secuencias;

    if (!ruta_archivo.empty()) {
        std::ifstream file(ruta_archivo);
        if (file.is_open()) {
            std::string line, current_seq = "";
            while (std::getline(file, line)) {
                if (line.empty()) continue;
                if (line[0] == '>') { // Cabecera FASTA
                    if (!current_seq.empty()) {
                        secuencias.push_back(current_seq);
                        current_seq = "";
                    }
                } else {
                    current_seq += line;
                }
            }
            if (!current_seq.empty()) secuencias.push_back(current_seq);
            file.close();
            std::cout << "[Tarea A] Cargadas " << secuencias.size() << " secuencias desde " << ruta_archivo << "\n";
            return secuencias;
        } else {
            std::cerr << "[Advertencia] No se pudo abrir el archivo, generando datos aleatorios...\n";
        }
    }

    // Generación sintética si no hay archivo FASTA
    secuencias.resize(N);
    const char bases[] = {'A', 'C', 'G', 'T'};
    srand(42); // Semilla fija para asegurar mediciones reproducibles en la práctica
    for (int i = 0; i < N; ++i) {
        secuencias[i].resize(L);
        for (int j = 0; j < L; ++j) {
            secuencias[i][j] = bases[rand() % 4];
        }
    }
    std::cout << "[Tarea A] Generadas " << N << " secuencias sinteticas de longitud " << L << "\n";
    return secuencias;
}

// ============================================================================
// TAREA B: Análisis de estadísticas globales (Secuencial)
// ============================================================================
Estadisticas calcular_estadisticas_globales(const std::vector<std::string>& secuencias) {
    Estadisticas est;
    long long total_bases = 0;

    for (const auto& seq : secuencias) {
        total_bases += seq.size();
        for (char c : seq) {
            switch (c) {
                case 'A': est.total_A++; break;
                case 'C': est.total_C++; break;
                case 'G': est.total_G++; break;
                case 'T': est.total_T++; break;
                default: break;
            }
        }
    }

    if (total_bases > 0) {
        est.porcentaje_gc = (double)(est.total_G + est.total_C) / total_bases * 100.0;
    }
    if (!secuencias.empty()) {
        est.longitud_media = (double)total_bases / secuencias.size();
    }
    return est;
}

// ============================================================================
// Rutina de alineamiento par a par (Needleman-Wunsch / Edit Distance)
// ============================================================================
float calcular_distancia_par(const std::string& s1, const std::string& s2) {
    int n = s1.size();
    int m = s2.size();

    // Matriz de programación dinámica aplanada en 1D para favorecer el acceso a la caché L1/L2
    std::vector<int> dp((n + 1) * (m + 1), 0);

    for (int i = 0; i <= n; ++i) dp[i * (m + 1)] = i;
    for (int j = 0; j <= m; ++j) dp[j] = j;

    for (int i = 1; i <= n; ++i) {
        int idx_actual = i * (m + 1);
        int idx_previo = (i - 1) * (m + 1);
        for (int j = 1; j <= m; ++j) {
            int coste = (s1[i - 1] == s2[j - 1]) ? 0 : 1;
            int ins_del = std::min(dp[idx_previo + j] + 1, dp[idx_actual + (j - 1)] + 1);
            dp[idx_actual + j] = std::min(ins_del, dp[idx_previo + (j - 1)] + coste);
        }
    }
    return static_cast<float>(dp[n * (m + 1) + m]) / static_cast<float>(n);
}

// ============================================================================
// TAREA C: Rellenado secuencial de la Matriz de Distancias (N x N)
// ============================================================================
void calcular_matriz_distancias(const std::vector<std::string>& secuencias, std::vector<float>& matriz, int N) {
    for (int i = 0; i < N; ++i) {
        matriz[i * N + i] = 0.0f; // Diagonal a cero
        for (int j = i + 1; j < N; ++j) {
            float dist = calcular_distancia_par(secuencias[i], secuencias[j]);
            matriz[i * N + j] = dist;
            matriz[j * N + i] = dist; // Simetría
        }
    }
}

// ============================================================================
// TAREA D: Agrupamiento jerárquico (Paso inicial del Árbol Filogenético)
// ============================================================================
void encontrar_par_mas_cercano(const std::vector<float>& matriz, int N, int& mejor_i, int& mejor_j, float& min_dist) {
    min_dist = 1e9f;
    mejor_i = -1;
    mejor_j = -1;

    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            float dist = matriz[i * N + j];
            if (dist < min_dist) {
                min_dist = dist;
                mejor_i = i;
                mejor_j = j;
            }
        }
    }
}

// ============================================================================
// TAREA E: Programa principal (Medición y salida)
// ============================================================================
int main(int argc, char* argv[]) {
    // Valores por defecto
    int N = 200;                // Número de secuencias
    int L = 1000;               // Longitud de las secuencias
    std::string ruta_fasta = ""; // Archivo FASTA opcional

    // Argumentos CLI para permitir ejecuciones sin intervención
    if (argc > 1) N = std::atoi(argv[1]);
    if (argc > 2) L = std::atoi(argv[2]);
    if (argc > 3) ruta_fasta = argv[3];

    std::cout << "========================================================\n";
    std::cout << "  Alineamiento de Secuencias y Arbol Filogenetico (Secuencial)\n";
    std::cout << "========================================================\n";
    std::cout << "Parametros de entrada: N = " << N << ", L = " << L << "\n";

    // Medición del tiempo de ejecución
    auto inicio = std::chrono::high_resolution_clock::now();

    // Tarea A: Cargar / Generar datos
    std::vector<std::string> secuencias = cargar_o_generar_secuencias(N, L, ruta_fasta);
    int num_secuencias = secuencias.size();

    // Tarea B: Estadísticas globales
    Estadisticas est = calcular_estadisticas_globales(secuencias);

    // Tarea C: Matriz de distancias N x N
    std::vector<float> matriz_distancias(num_secuencias * num_secuencias, 0.0f);
    calcular_matriz_distancias(secuencias, matriz_distancias, num_secuencias);

    // Tarea D: Búsqueda del par más cercano en la matriz
    int mejor_i = -1, mejor_j = -1;
    float min_dist = 0.0f;
    encontrar_par_mas_cercano(matriz_distancias, num_secuencias, mejor_i, mejor_j, min_dist);

    auto fin = std::chrono::high_resolution_clock::now();
    double tiempo_total = std::chrono::duration<double>(fin - inicio).count();

    // Salida final
    std::cout << "\n--- RESULTADOS DE LA EJECUCION ---\n";
    std::cout << "[Tarea B] Longitud media: " << est.longitud_media << " bases\n";
    std::cout << "[Tarea B] Contenido GC: " << std::fixed << std::setprecision(2) << est.porcentaje_gc << "%\n";
    std::cout << "[Tarea D] Primeras secuencias en agruparse en el arbol: (" 
              << mejor_i << ", " << mejor_j << ") con distancia = " << min_dist << "\n";
    std::cout << "[Tarea E] Tiempo total de ejecucion secuencial: " 
              << std::setprecision(4) << tiempo_total << " segundos.\n";
    std::cout << "========================================================\n";

    return 0;
}