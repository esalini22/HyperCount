#include "HyperCount.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <omp.h>
#include <zlib.h>
#include <algorithm>
#include "wyhash32.h"

#include <map>

using namespace std;
#define seed 5 //seed para hash
#define lim 4294967296 //2^32

inline float HyperCount::hsum_sse3(__m128 v) {
    __m128 shuf = _mm_movehdup_ps(v);        // broadcast elements 3,1 to 2,0
    __m128 maxs = _mm_add_ps(v, shuf);
    shuf        = _mm_movehl_ps(shuf, maxs); // high half -> low half
    maxs        = _mm_add_ss(maxs, shuf);
    return        _mm_cvtss_f32(maxs);
}

inline float HyperCount::hsum_avx(__m256 v) {
    __m128 lo = _mm256_castps256_ps128(v);   // low 128
    __m128 hi = _mm256_extractf128_ps(v, 1); // high 128
           lo = _mm_add_ps(lo, hi);          // max the low 128
    return hsum_sse3(lo);                    // and inline the sse3 version
}

HyperCount::HyperCount(unsigned char n1, unsigned char n2){
	srand(time(NULL));

	//omp_set_num_threads(numThreads);
	p=n1;
	b=n2;
	bits_v2=(1<<b)-1; //2^b -1
	N=1<<p; //2^V1
	for(unsigned char i=0;i<16;++i)
		bit_mask.emplace_back(~((ullint)0xF<<(4*i)));
	sketch_size=to_string(p);
	min_val=1<<(b-15);
	cont=0;
}
HyperCount::~HyperCount(){}

void HyperCount::addSketch(string genome){
	//inicializa sketch en 0s, y añade su nombre a un vector
	for(unsigned int j=0;j<N/16;++j) //cada celda tendra 12 buckets del sketch array
		sketch.emplace_back(0);
	for(unsigned int j=0;j<N;++j) //probar explotando con loop unrolling
		aux.emplace_back(0);
	for(unsigned int j=0;j<N;++j)
		freq.emplace_back(0);
	for(unsigned int j=0;j<N;++j)
		demo.emplace_back(0);
	//if((long)N%16)
		//sketch.emplace_back(0);
	name=genome;
	printf("%s\n",name.c_str());
}

void HyperCount::insert(ullint kmer){
	cont++;
	//calculamos el hash usando wyhash de 32 bits (hashing muy rapido)
	const ullint *key = &kmer;
	uint32_t hash=wyhash32(key,8,seed);

	unsigned int v1=hash>>b; //se desplaza el hash b bits a la derecha
	unsigned int v2=hash&bits_v2; // AND 2^b -1 
	unsigned char w;
	if(v2<min_val) w=15;
	else w=__builtin_clz(v2)+1-p; //cuenta cantidad de bits 0 mas significativo
	//lo cual sirve para encontrar posicion de 1 mas a la izquierda
	//se resta p, ya que al trabajar con int, se tiene 32 bits,
	//y siempre habran al menos p ceros a la izquierda
	unsigned int indice=v1/16; //indice de celda
	unsigned char bucket=v1%16; //12 buckets por celda
	ullint temp=(sketch[indice]>>(4*bucket))&0xF;

	if(temp < w || (temp==w && aux[v1]<v2)){ //
		aux[v1]=v2;
		freq[v1]=1;
		//if(demo[v1]!=0) counter[demo[v1]]--;
		//if(counter[demo[v1]]==0) counter.erase(demo[v1]);
		demo[v1]=kmer;
		//pair<unordered_map<ullint,int>::iterator,bool> ret=counter.insert(pair<ullint,int>(kmer,1));
		//if(ret.second==false) counter[kmer]++;
	}
	else if(temp==w && aux[v1]==v2){
		freq[v1]++;
		if(demo[v1]!=kmer){ //nos quedamos con el mas repetido
			/*if(counter[kmer]>counter[demo[v1]]){
				counter[demo[v1]]--;
				if(counter[demo[v1]]==0) counter.erase(demo[v1]);
				demo[v1]=kmer;
				counter[kmer]++;
			}*/
			int cont_demo=0,cont_kmer=0;
			for(unsigned int j=0;j<N;++j){
				if(demo[j]==demo[v1]) cont_demo++;
				if(demo[j]==kmer) cont_kmer++;
			}
			if(cont_kmer>cont_demo) demo[v1]=kmer;
		}
	}

	else if(temp>w && aux[v1]>v2){ //linea añadida por mi, mejora significativamente la precision de la entropia de la muestra
		freq[v1]++;
	}

	if(w > temp) sketch[indice]=(sketch[indice]&(bit_mask[bucket]))|((ullint)v2<<(4*bucket));
}

void HyperCount::loadSketch(char *infilename){ //¿que pasa si tamaño del sketch es distinto al de argumento?
	string genome=infilename; //se debe determinar nombre del genoma a partir del archivo .hll
	printf("infilename: %s\n",infilename);
	genome.erase(genome.find(".w."),string::npos);
	printf("genome: %s\n",genome.c_str());

	//debemos descomprimir el archivo y ponerlo en un vector
	vector<ullint> v;
	gzFile infile = gzopen(infilename, "rb");
	//if(!infile) return -1;
	char buffer[128];
	int num_read = 0;
	char temp[19];
	char* end;
	int cont=0;
	while((num_read = gzread(infile,buffer, sizeof(buffer))) > 0){
		for(int i=0;i<num_read;++i){
			temp[cont%19]=buffer[i];
			++cont;
			if(cont%19==0){
				cont=0;
				v.emplace_back(strtoull(temp,&end,10));
			}
		}
	}
	if(cont>0) v.emplace_back(strtoull(temp,&end,10));
	v.emplace_back(0);
	gzclose(infile);

	sketch=v;
	name=genome;
	printf("done\n");
}

void HyperCount::saveSketch(){ //usar zlib para comprimir
	printf("save sketch\n");
	string temp_name=name;
	temp_name+=".spacing."+sketch_size+".hll";
	char* outfilename=(char*)temp_name.c_str();
	printf("%s\n",outfilename);
	gzFile outfile = gzopen(outfilename, "wb");
	//if(s.empty() || !outfile) return -1;
	char inbuffer[128];
	char cont=0;
	for(vector<ullint>::iterator it=sketch.begin();it!=sketch.end();++it){
		string temp=to_string(*it)+'\n';
		int len=temp.length();
		for(char i=0;i<len;++i){
			inbuffer[cont]=temp[i];
			++cont;
			if(cont==127){
				gzwrite(outfile, inbuffer, cont);
				cont=0;
			}
		}
	}
	if(cont>0) gzwrite(outfile, inbuffer, cont);
	gzclose(outfile);
}

vector<ullint> HyperCount::merge(HyperCount *hll){ //hace la union y devuelve un nuevo sketch
	vector<ullint> sketch2=hll->getSketch();
	vector<ullint> ret=max(sketch,sketch2);
	ret.resize(max(sketch.size(),sketch2.size()));
	vector<ullint>::iterator it2=ret==sketch?sketch2.begin():sketch.begin();
	
	/*for(unsigned int j=0;j<N;++j){
		if(){
		}
		else if(){
		}
		else{
		}
	}*/

	return ret;
}

const vector<ullint> &HyperCount::getSketch(){
	return sketch;
}
string HyperCount::getName(){
	return name;
}

void HyperCount::print(){
	unordered_map<ullint,int> sumcont;
	for(unsigned int j=0;j<N;++j){ //sumamos las frecuencias
		pair<unordered_map<ullint,int>::iterator,bool> ret = sumcont.insert(pair<ullint,int>(demo[j],freq[j]));
		if(ret.second==false) sumcont[demo[j]]+=freq[j];
	}
	multimap<int,ullint,greater<int>> estcont; 
	//ordenamos por frecuencia
	for(unordered_map<ullint,int>::iterator it=sumcont.begin();it!=sumcont.end();++it){
		if(it->first!=0) estcont.insert(pair<int,ullint>(it->second,it->first));
	}
	for(multimap<int,ullint>::iterator it=estcont.begin();it!=estcont.end();++it){
		printf("%d:	%llu\n",it->first,it->second);
	}

}

float HyperCount::entropy(){
	unordered_map<ullint,int> table;
	int size=0;
	for(unsigned int j=0;j<N;++j){
		if(freq[j]==0) continue;
		pair<unordered_map<ullint,int>::iterator,bool> ret=table.insert(pair<ullint,int>(demo[j],freq[j]));
		if(!ret.second) table[demo[j]]+=freq[j];
		size+=freq[j];
	}
	float ent=0;
	for(unordered_map<ullint,int>::iterator it=table.begin();it!=table.end();++it){
		ent+=((float)it->second/(float)size)*log2((float)it->second/(float)size);
	}
	return -(float)ent*(float)size*9.3/(float)(log2(table.size())*cont); //entropia normalizada
}
