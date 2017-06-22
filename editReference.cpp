#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <queue>
using namespace std;
// clang++ -std=c++11 -stdlib=libc++ editReference.cpp -o editReference


struct Mutation{
  int chr, position;
  char base;
};


const double AF_THRESHOLD=0.9; //allele freq 
const int HEADER_CHARS=6; //number of bases on a header line in hg19 (+1 to account for \n)
const int SEQUENCE_CHARS=51; //number of bases on a line in hg19 (+1 to account for \n)
const int CHR_SIZES[] = {16900, 1321120100, 675007400, 201982800, 236608100, 184533500, 213748200, 162508000, 149368900, 185015300, 138245400, 178191400, 136528900, 117473200, 109496500, 104582000, 92161800, 88624000, 79643200, 100409800, 64286000, 88933800, 52330600, 158375900, 60561000};

fstream fRef;
ifstream fNormal;


bool aboveThreshold(int depth, int totalPairs){
  if((double)totalPairs/(double)depth > AF_THRESHOLD)
    return true;
  return false;
}


queue<Mutation> normals; //assumes freq file is sorted
int normalPairs[8];
void storeNormals(){
  string filler;
  fNormal >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler;

  while(!fNormal.eof()){
    string chr;
    int pos, depth;
    fNormal >> chr >> pos >> depth >> filler >> filler >> filler >> normalPairs[0] >> normalPairs[1] >> normalPairs[2] >> normalPairs[3] >> normalPairs[4] >> normalPairs[5] >> normalPairs[6] >> normalPairs[7];
    if(chr == "") continue; 

    int chrNum=-1;
    chr = chr.substr(3);
    if(chr == "X")
      chrNum = 23;
    else if(chr == "Y")
      chrNum = 24;
    else 
      chrNum = stoi(chr);

    if(aboveThreshold(depth, normalPairs[0]+normalPairs[1])) 
      normals.push({chrNum, pos, 'A'});
    else if(aboveThreshold(depth, normalPairs[2]+normalPairs[3])) 
      normals.push({chrNum, pos, 'C'});
    else if(aboveThreshold(depth, normalPairs[4]+normalPairs[5]))
      normals.push({chrNum, pos, 'T'});
    else if(aboveThreshold(depth, normalPairs[6]+normalPairs[7])) 
      normals.push({chrNum, pos, 'G'});
  }
}


void editRef(){
  int chrCtr=-1, sequencePos=-1; //chrCtr starts at -1 because chrM is the first chromosome in hg19
  char c;
  while(!normals.empty()){
    Mutation current = normals.front();

    // cout << current.chr << "-" << current.position << " - " << current.base << "  (" << chrCtr << " - " << sequencePos << ')' <<  '\n';
    // if(current.chr ==1 && current.position > 1289991)
    //   return;
    if(current.chr == chrCtr && current.position-sequencePos<SEQUENCE_CHARS-1){ //mutation is on current line or nearing the end of the chromosome (where the SEQUENCE_CHARS is not constant)
      for(int i=0; i<SEQUENCE_CHARS-1; i++){
        if(current.position == i+sequencePos){
          fRef.put(current.base);
          normals.pop();
          current = normals.front();
        }
        else
          fRef.seekg(1, ios::cur);
      }
      fRef.seekg(1, ios::cur);
      sequencePos+=SEQUENCE_CHARS-1;
    }
    else if(fRef.peek() == '>'){ //on header line
      fRef.seekg(HEADER_CHARS, ios::cur);
      chrCtr++;
      sequencePos=1;
    }
    else{
      fRef.seekg(SEQUENCE_CHARS, ios::cur); //skips line
      sequencePos+=SEQUENCE_CHARS-1;
    }

  }
}


int main(int argc, char *argv[]){ //<ref_genome> <normal_freq>
  std::ios::sync_with_stdio(false); cin.tie(NULL);
  fRef.open(argv[1]);
  fNormal.open(argv[2]);

  clock_t t;
  t = clock();

  // cout << "----Storing normal mutations..." << '\n';
  // storeNormals();
  normals.push({1, 15, '.'});
  normals.push({2, 45, '.'});
  normals.push({2, 115, '.'});
  normals.push({3, 100, '.'});
  cout << "----Substituting normal mutations into reference genome..." << '\n';
  editRef();

  cout << "----Done" << '\n';
  cout << "------Executed in " << ((float)(clock()-t))/CLOCKS_PER_SEC << " seconds." << '\n';

  fRef.close();
  fNormal.close();
  return 0;
}