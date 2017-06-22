#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>
using namespace std;
// clang++ -std=c++11 -stdlib=libc++ editReference.cpp -o editReference


const double AF_THRESHOLD=0.9;


fstream fRef;
ifstream fNormal;


bool aboveThreshold(int depth, int totalPairs){
  if((double)totalPairs/(double)depth > AF_THRESHOLD)
    return true;
  return false;
}


unordered_map<string, char> normals;
int normalPairs[8];
void storeNormals(){
  string filler;
  fNormal >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler;

  while(!fNormal.eof()){
    string chr;
    int pos, depth;
    fNormal >> chr >> pos >> depth >> filler >> filler >> filler >> normalPairs[0] >> normalPairs[1] >> normalPairs[2] >> normalPairs[3] >> normalPairs[4] >> normalPairs[5] >> normalPairs[6] >> normalPairs[7];
    if(chr == "") continue; 

    chr += '-'; chr += to_string(pos); 
    if(aboveThreshold(depth, normalPairs[0]+normalPairs[1])) 
      normals[chr] = 'A';
    if(aboveThreshold(depth, normalPairs[2]+normalPairs[3])) 
      normals[chr] = 'C';
    if(aboveThreshold(depth, normalPairs[4]+normalPairs[5]))
      normals[chr] = 'T';
    if(aboveThreshold(depth, normalPairs[6]+normalPairs[7])) 
      normals[chr] = 'G';
  }
}


void editRef(){
  char c;
  while(fRef.get(c)){
    // fRef.seekg(-1, ios::cur);
    // fRef.put(c);
  }
}

int main(int argc, char *argv[]){ //<ref_genome> <normal_freq>
  std::ios::sync_with_stdio(false); cin.tie(NULL);
  fRef.open(argv[1]);
  fNormal.open(argv[2]);

  clock_t t;
  t = clock();

  cout << "----Storing normal mutations..." << '\n';
  storeNormals();
  cout << "----Substituting normal mutations into ref genome..." << '\n';
  editRef();

  cout << "----Done" << '\n';
  cout << "------Executed in " << ((float)(clock()-t))/CLOCKS_PER_SEC << " seconds." << '\n';

  fRef.close();
  fNormal.close();
  return 0;
}