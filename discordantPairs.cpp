#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>
#include <vector>
#include <algorithm>
#include <assert.h>
using namespace std;
// sh local_discordant_c_r.sh

struct Read{
  string chr, geneName;
  int pos, geneStart;
};
struct Exon{
  string name;
  int start, end;
};


ifstream fLeft, fRight, fExons;
ofstream fDetails, fDepth;


void finishLine(ifstream &fin){
  char c;
  fin.get(c);
  while(c != '\n')
    fin.get(c);
}

bool exonCmpr(const Exon &a, const Exon &b){
  return a.start < b.start;
}


int chrToInt(string chr){
  if(chr.size() > 7 || chr=="chrM") return -1;
  if(chr=="chrX") return 23;
  if(chr=="chrY") return 24;
  return stoi(chr.substr(3));
}


vector<vector<Exon> > exons(25, vector<Exon>(1)); //[chr][exon]
void storeExons(){
  string chr, name;
  int start, end;
  while(fExons >> chr >> start >> end >> name){
    int chrNum = chrToInt(chr);
    if(chrNum == -1)
      continue;
    
    exons[chrNum].push_back({name, start, end});
  }
  for(int i=1; i<=24; i++)
    sort(exons[i].begin(), exons[i].end(), exonCmpr);
}


int getDist(int pos, int start, int end){
  if(pos >= start && pos <= end)
    return 0;
  if(pos > end)
    return pos-end;
  return start-pos;
}


pair<string, int> findGene(int chr, int pos){
  vector<Exon> &chromosome = exons[chr];
  auto it = lower_bound(chromosome.begin(), chromosome.end(), pos,
                        [](const Exon &exon, const int value){  return exon.start < value; });

  int minDist = INT_MAX;
  pair<string, int> closest;
  for(int i=-5; i<=5; i++)
    if(it+i >= chromosome.begin() && it+i <= chromosome.end()){
      int dist = getDist(pos, (it+i)->start, (it+i)->end);
      if(dist < minDist){
        closest = {(it+i)->name, (it+i)->start};
        minDist = dist;
      }
    }
  
  assert(minDist != INT_MAX && "Error: No genes found in chromosome.");
  return closest;
}


string alphabetize(string s1, string s2){
  if(s1 < s2)
    return s1 + '-' + s2;
  return s2 + '-' + s1;
}

  
unordered_map<string, int> fusionCtr; //eg. {MXRA8-PIP5K1B, 12} - how many times a discordant pair combination (eg chr1/chr2) happens
unordered_map<string, vector<Read> > molecules; //{D87PMJN1:331:C528CACXX:7:2205:14374:10502:TGCCATCG, {read1, read2}} - keeps track of the different chromosomes a read-pair has seen so far
void generateCandidates(int minMapQ, ifstream &fin){
  string header, filler, chr;
  int mQ, pos;
  while(fin >> header >> filler >> chr >> pos >> mQ){
    finishLine(fin);
    if(mQ < minMapQ || chr.size() > 6 || chr == "chrM") continue; 

    if(header[header.size()-2] == '/'){
      header.pop_back(); header.pop_back(); 
    } 
    pair<string, int> curGene = findGene(chrToInt(chr), pos);

    unordered_map<string, vector<Read> >::iterator it = molecules.find(header);
    if(it != molecules.end()){ 
      for(Read r : it->second){
        if(r.chr == chr && r.geneName == curGene.first)
           continue;

        if(r.geneName < curGene.first) //alphabetized output
          fDetails << r.geneName << '\t' << curGene.first << '\t' << r.geneName << '-' << curGene.first << '\t' << r.chr << '\t' << chr << '\t' << r.geneStart << '\t' << curGene.second << '\t' << r.pos << '\t' << pos << '\n'; 
        else 
          fDetails << curGene.first << '\t' << r.geneName << '\t' << curGene.first << '-' << r.geneName << '\t' << chr << '\t' << r.chr << '\t' << curGene.second << '\t' << r.geneStart << '\t' << pos << '\t' << r.pos << '\n'; 

        pair<unordered_map<string, int>::iterator, bool> ret = fusionCtr.emplace(alphabetize(curGene.first, r.geneName), 1); //found discordant reads
        if(ret.second == false) 
          (ret.first->second)++; //keep track of how many times certain discordant read pairs show up
      }
      it->second.push_back({chr, curGene.first, pos, curGene.second});
    }
    else{
      vector<Read> vec = {{chr, curGene.first, pos, curGene.second}};
      molecules.emplace(header, vec);
    }

  }
}


void printDepths(){
  for(auto p : fusionCtr)
    fDepth << p.first << '\t' << p.second << '\n';
}


int main(int argc, char* argv[]){ // <left> <right> <allexons> <min_mapping_quality> <details> <depth> 
  std::ios::sync_with_stdio(false); cin.tie(NULL);
  
  fLeft.open(argv[1]);
  fRight.open(argv[2]);
  fExons.open(argv[3]);
  fDetails.open(argv[5]);
  fDepth.open(argv[6]);


  fDetails << "GENE1\tGENE2\tFUSION\tCHR1\tCHR2\tGENE1_START\tGENE2_START\tREAD1_POS\tREAD2_POS\n";
  storeExons();
  generateCandidates(atoi(argv[4]), fLeft);
  generateCandidates(atoi(argv[4]), fRight);

  fDepth << "FUSION\tDEPTH\n";
  printDepths();

  fLeft.close();
  fRight.close();
  fExons.close();
  fDetails.close();
  fDepth.close();
  return 0;
}