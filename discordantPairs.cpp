#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>
#include <vector>
using namespace std;
// clang++ -std=c++11 -stdlib=libc++ discordantPairs.cpp -o discordantPairs
// ./discordantPairs DLBCL021_Tumor.remapped.left.txt DLBCL021_Tumor.remapped.right.txt 40 DLBCL021_Tumor.discordantpairs.txt


struct Read{
  string chr, geneName;
  int pos, geneStart;
};


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

  
unordered_map<string, int> fusionCtr; //eg. {MXRA8-PIP5K1B, 12} - how many times a discordant pair combination (eg chr1/chr2) happens
unordered_map<string, vector<Read> > molecules; //{D87PMJN1:331:C528CACXX:7:2205:14374:10502:TGCCATCG, {read1, read2}} - keeps track of the different chromosomes a read-pair has seen so far
void generateCandidates(int minMapQ, ifstream &fin){
  string header, filler, chr;
  int mQ;
  while(fin >> header >> filler >> chr >> filler >> mQ){
    finishLine(fin);
    if(mQ < minMapQ || chr.size() > 6 || chr == "chrM") continue; 

    if(header[header.size()-2] == '/'){
      header.pop_back(); header.pop_back(); 
    } 

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


int main(int argc, char* argv[]){ // <left> <right> <min_mapping_quality> <output> 
  std::ios::sync_with_stdio(false); cin.tie(NULL);
  
  fLeft.open(argv[1]);
  fRight.open(argv[2]);
  fOut.open(argv[4]);

  fOut << "GENE1\tGENE2 FUSION\tCHR1\tCHR2\tGENE1_START\tGENE2_START\tDIST READ1_POS\tREAD2_POS\n";
  generateCandidates(atoi(argv[3]), fLeft);
  generateCandidates(atoi(argv[3]), fRight);
  printCandidates();

  fLeft.close();
  fRight.close();
  fOut.close();
  return 0;
}