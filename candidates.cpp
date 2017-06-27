#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>
#include <vector>
using namespace std;
// clang++ -std=c++11 -stdlib=libc++ candidates.cpp -o candidates
// ./candidates DLBCL021_Tumor.remapped.left.txt DLBCL021_Tumor.remapped.right.txt 20 DLBCL021_Tumor.candidates.txt


ifstream fLeft, fRight;
ofstream fOut;


string alphabetize(string s1, string s2){
  if(s1 < s2 || s1.size() < s2.size())
    return s1 + '/' + s2;
  return s2 + '/' + s1;
}


void finishLine(ifstream &fin){
  char c;
  fin.get(c);
  while(c != '\n')
    fin.get(c);
}

  
unordered_map<string, int> fusionCtr; //{discordant_reads, counter} - how many times a discordant pair combination (eg chr1/chr2) happens
unordered_map<string, vector<string> > molecules; //{header, chromosomes} - keeps track of the different chromosomes a read-pair has seen so far
void generateCandidates(int minMapQ, ifstream &fin){
  string header, filler, chr;
  int mQ;
  while(fin >> header >> filler >> chr >> filler >> mQ){
    finishLine(fin);
    if(mQ < minMapQ || chr.size() > 6 || chr == "chrM") continue;  

    unordered_map<string, vector<string> >::iterator it = molecules.find(header);
    if(it != molecules.end()){ //check if this fragment is the first of 4 possible (2 for each read pair)
      for(string s : it->second){
        if(s == chr) continue; //2 chromosomes mapped are the same, no fusion

        pair<unordered_map<string, int>::iterator, bool> ret = fusionCtr.emplace(alphabetize(chr, s), 1); //found discordant reads
        if(ret.second == false) 
          (ret.first->second)++; //keep track of how many times certain discordant read pairs show up
      }
      it->second.push_back(chr);
    }
    else{
      vector<string> vec = {chr};
      molecules.emplace(header, vec);
    }

  }
}


void printCandidates(){
  for(auto p : fusionCtr)
    cout << p.second << ' ' << p.first << '\n';
}


int main(int argc, char* argv[]){ // <left> <right> <min_mapping_quality> <output> 
  std::ios::sync_with_stdio(false); cin.tie(NULL);
  
  fLeft.open(argv[1]);
  fRight.open(argv[2]);
  fOut.open(argv[4]);

  generateCandidates(atoi(argv[3]), fLeft);
  generateCandidates(atoi(argv[3]), fRight);
  printCandidates();

  fLeft.close();
  fRight.close();
  fOut.close();
  return 0;
}