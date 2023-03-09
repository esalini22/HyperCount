#ifndef HYPERCOUNT_H
#define HYPERCOUNT_H
#include <vector>
#include <unordered_map>
#include <string>
#include <immintrin.h> 
typedef unsigned long long int ullint;
using namespace std;
class HyperCount{
    private:
        unsigned char p,b;
        unsigned int bits_v2; //bits para hacer & en operacion bitwise
        long double N; //corresponde a |M|: numero de buckets - se usa long double para evitar hacer casting al calcular cardinalidad
        float a_m;
        int ciclos_red;
        string kmer_length;
        string sketch_size;
        vector<ullint> bit_mask; //vector de máscaras de bits, usadas para insertar registros
        //se usa vector de ullint, cada celda (64 bits) a su vez tiene hasta 12 buckets
        vector<ullint> sketch; //n sketches, n: numero de genomas
        vector<int> aux; //vector de indicadores auxiliares
        vector<int> freq; //vector de contadores de frecuencias
        vector<ullint> demo; //vector de muestras demograficas
        string name; //vector con los nombres de los genomas
        unsigned int min_val; //minimo valor en v2 para no sobrepasar registro 15
        //unordered_map<ullint,int> counter; //cuenta repeticiones de demo en sketch
    public:
        HyperCount(unsigned char n1,unsigned char n2,unsigned char k);
        ~HyperCount();
        inline float hsum_sse3(__m128 v);
        inline float hsum_avx(__m256 v);
        void addSketch(string genome); //añade un sketch vacio al mapa de sketches, para funciones hll y dist
        void insert(ullint kmer); //inserta en sketch, para funciones hll y dist
        void loadSketch(char* infilename); //carga un archivo de sketch en memoria, para hll y dist
        void saveSketch(); //guarda el sketch individual en un archivo, para funcion sketch
        //float estCard(); //cardinalidad individual
        vector<ullint> merge(HyperCount *hll); //hace la union y devuelve un nuevo sketch
        const vector<ullint> &getSketch();
        string getName();
        //void saveOutput(char* filename); //guarda las cardinalidades y distancias en un txt
        void print();
};
#endif
