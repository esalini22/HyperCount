#include <stdio.h>
#include <string.h>
#include <string>
#include <fstream>
#include <unistd.h>
#include <cmath>
#include <omp.h>
#include <thread>
#include <vector>
#include <algorithm>
#include <vector>
#include <immintrin.h> 
#include "HyperCount.h"

#include <unordered_map>
#include <map>

#define lim 4294967296 //2^32
#define lim_int 2147483647

using namespace std;
typedef unsigned long long int ullint;

unordered_map<ullint,int> counter;
multimap<int,ullint,greater<int>> truecont;
int elements;

void leer(char *genome,HyperCount *hll){
	char c; //para lectura
	ullint kmer=0,comp=0;
	string linea;
	ifstream indata(genome);
	if(!indata){
		printf("No se pudo abrir el archivo %s\n",genome);
		exit(1);
	}

	hll->addSketch(genome); //añade el sketch inicializado con 0s

	getline(indata,linea); //para primer kmer

	//al leerse cada caracter, se insert la base al final del kmer no canonico
	//y luego se insert su complemento al inicio del complemento del reverso del kmer
	//el inicio estara dado por el largo del kmer
	hll->insert(atoi(linea.c_str()));
	pair<unordered_map<ullint,int>::iterator,bool> ret=counter.insert(pair<ullint,int>(atoi(linea.c_str()),1));
	elements++;
	if(!ret.second) counter[atoi(linea.c_str())]++;
	while(!indata.eof()){
		getline(indata,linea);
		hll->insert(atoi(linea.c_str()));
		pair<unordered_map<ullint,int>::iterator,bool> ret=counter.insert(pair<ullint,int>(atoi(linea.c_str()),1));
		if(!ret.second) counter[atoi(linea.c_str())]++;
		elements++;
	}
	indata.close();
}

//obtiene archivos de un txt
//formato de archivos:
//genoma1
//genoma2
//etc. (1 genoma por linea)
vector<string> readFromFile(char* paths){
	ifstream indata(paths);
	if(!indata){
		printf("No se pudo abrir el archivo %s\n",paths);
		exit(1);
	}
	vector<string> genomes;
	string filename;
	while(!indata.eof()){
		getline(indata,filename);
		if(filename!="") genomes.push_back(filename);
	}
	indata.close();
	return genomes;
}

//obtiene archivos de la linea de argumentos
vector<string> getPaths(char** argv, int argc){
	vector<string> genomes;
	for(int i=1;i<argc;++i){
		if(!strcmp(argv[i],"-k") || !strcmp(argv[i],"-p") || !strcmp(argv[i],"-t") || !strcmp(argv[i],"-o") || !strcmp(argv[i],"-d") || !strcmp(argv[i],"-r")) ++i;
		else if(strcmp(argv[i],"-s")) genomes.push_back(argv[i]);
	}
	return genomes;
}

vector<string> readCompressedFromFile(char* paths){
	ifstream indata(paths);
	if(!indata){
		printf("No se pudo abrir el archivo %s\n",paths);
		exit(1);
	}
	vector<string> genomes;
	string filename;
	while(!indata.eof()){
		getline(indata,filename);
		if(filename!="") genomes.push_back(filename);
	}
	indata.close();
	return genomes;
}

//obtiene archivos de la linea de argumentos
vector<string> getCompressed(char** argv, int argc){
	vector<string> genomes;
	for(int i=1;i<argc;++i){
		if(!strcmp(argv[i],"-p") || !strcmp(argv[i],"-t") || !strcmp(argv[i],"-o") || !strcmp(argv[i],"-f") || !strcmp(argv[i],"-r")) ++i;
		else if(!strcmp(argv[i],"-d")) genomes.push_back(argv[i+1]);
	}
	return genomes;
}

//formato: ./hll -opcion valor genomas
//o bien ./hll genomas -opcion valor
//no detecta caso en que se introduza opcion o valor invalido
int main(int argc, char *argv[]){
	if(argc<2) {
		printf("No hay suficientes argumentos\n");
		exit(1);
	}
	unsigned char p=12;
	
	char** option;
	char** end=argv+argc;
	option=std::find((char**)argv,end,(const std::string&)"-p");
	if(option!=end){
		char val=atoi(*(option+1));
		if(val<16 && val>8) p=val;
	}
	
	vector<string> genomes,compressed;
	option=std::find((char**)argv,end,(const std::string&)"-f");
	if(option!=end) genomes=readFromFile((char*)(*(option+1)));
	else genomes=getPaths(argv,argc);

	option=std::find((char**)argv,end,(const std::string&)"-r");
	if(option!=end) compressed=readCompressedFromFile((char*)(*(option+1)));
	else compressed=getCompressed(argv,argc);

	int tam=genomes.size(),tam2=compressed.size();
	printf("tam: %d, tam2: %d\n",tam,tam2);
	printf("p: %d\n",p);
	
	int numThreads=min(tam+tam2,(int)std::thread::hardware_concurrency());

	option=std::find((char**)argv,end,(const std::string&)"-t");
	if(option!=end) numThreads=atoi((*(option+1)));

	printf("threads: %d\n",numThreads);

	vector<HyperCount*> v_hll;
	for(int i=0;i<tam+tam2;++i){
		HyperCount *hll;
		hll = new HyperCount(p,32-p);
		v_hll.push_back(hll);
	}

	elements=0;

	//lee paralelamente cada archivo
	omp_set_num_threads(numThreads);
	#pragma omp parallel
	{
		#pragma omp single
		{
			for(int i=0;i<tam;i++){
				#pragma omp task
				leer((char*)genomes[i].c_str(),v_hll[i]);
			}
			for(int i=0;i<tam2;i++){
				#pragma omp task
				v_hll[i+tam]->loadSketch((char*)compressed[i].c_str());
			}
		}
	}
	
	option=std::find((char**)argv,end,(const std::string&)"-s");
	if(option!=end){
		for(int i=0;i<tam+tam2;++i)
			v_hll[i]->saveSketch();
	}


	//for(int i=0;i<tam+tam2;++i)
		//v_hll[i]->print();


	/*printf("==============================================================\n");
	for(unordered_map<ullint,int>::iterator it=counter.begin();it!=counter.end();++it){
		truecont.insert(pair<int,ullint>(it->second,it->first));
	}
	for(multimap<int,ullint>::iterator it=truecont.begin();it!=truecont.end();++it){
		printf("%d:	%llu\n",it->first,it->second);
	}*/

	for(int i=0;i<tam+tam2;++i)
		printf("Entropy estimate: %f\n",v_hll[i]->entropy());

	float true_entropy=0;
	for(unordered_map<ullint,int>::iterator it=counter.begin();it!=counter.end();++it){
		true_entropy+=((float)it->second/(float)elements)*log2((float)it->second/(float)elements);
	}
	printf("True entropy: %f\n",-(float)true_entropy/(float)log2(counter.size()));

	for(int i=0;i<tam+tam2;++i)
		delete v_hll[i];

	return 0;
}
