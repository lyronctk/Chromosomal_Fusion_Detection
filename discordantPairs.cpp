#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;


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
void storeExons(int fragmentLength){
  string chr, name;
  int start, end;
  while(fExons >> chr >> start >> end >> name){
    int chrNum = chrToInt(chr);
    if(chrNum == -1)
      continue;
    
    exons[chrNum].push_back({name, start-fragmentLength+1, end}); // subtract fragmentLength from start because the fragment does not need to completely be within the exon
  }
  for(int i=1; i<=24; i++)
    sort(exons[i].begin(), exons[i].end(), exonCmpr);
}


pair<string, int> findGene(int chr, int pos, int fragmentLength){
  vector<Exon> &chromosome = exons[chr];
  auto it = lower_bound(chromosome.begin(), chromosome.end(), pos,
                        [](const Exon &exon, const int value){  return exon.start < value; });

  for(int i=-5; i<=5; i++)
    if(pos >= (it+i)->start && pos <= (it+i)->end && it >= chromosome.begin() && it <= chromosome.end())
      return {(it+i)->name, (it+i)->start+fragmentLength-1};
  return {"NON_EXONIC", 0};
}


string alphabetize(string s1, string s2){
  if(s1 < s2)
    return s1 + '-' + s2;
  return s2 + '-' + s1;
}

  
unordered_map<string, int> fusionCtr; //eg. {MXRA8-PIP5K1B, 12} - how many times a discordant pair combination (eg chr1/chr2) happens
unordered_map<string, vector<Read> > molecules; //{D87PMJN1:331:C528CACXX:7:2205:14374:10502:TGCCATCG, {read1, read2}} - keeps track of the different chromosomes a read-pair has seen so far
void generateCandidates(int minMapQ, int fragmentLength, ifstream &fin){
  string header, filler, chr;
  int mQ, pos;
  while(fin >> header >> filler >> chr >> pos >> mQ){
    finishLine(fin);
    if(mQ < minMapQ || chr.size() > 6 || chr == "chrM") continue; 

    if(header[header.size()-2] == '/'){
      header.pop_back(); header.pop_back(); 
    } 
    pair<string, int> curGene = findGene(chrToInt(chr), pos, fragmentLength);

    unordered_map<string, vector<Read> >::iterator it = molecules.find(header);
    if(it != molecules.end()){ //check if this fragment is the first of 4 possible (2 for each read pair)
      bool isDup = false;
      for(Read r : it->second){
        if(r.chr == chr){ //2 chromosomes mapped are the same, no fusion
          if(r.pos == pos)
            isDup = true;
          continue;  
        }

        if(r.geneName < curGene.first) //alphabetized output
          fDetails << r.geneName << '\t' << curGene.first << '\t' << r.geneName << '-' << curGene.first << '\t' << r.chr << '\t' << chr << '\t' << r.geneStart << '\t' << curGene.second << '\t' << r.pos << '\t' << pos << '\n'; 
        else 
          fDetails << curGene.first << '\t' << r.geneName << '\t' << curGene.first << '-' << r.geneName << '\t' << chr << '\t' << r.chr << '\t' << curGene.second << '\t' << r.geneStart << '\t' << pos << '\t' << r.pos << '\n'; 

        pair<unordered_map<string, int>::iterator, bool> ret = fusionCtr.emplace(alphabetize(curGene.first, r.geneName), 1); //found discordant reads
        if(ret.second == false) 
          (ret.first->second)++; //keep track of how many times certain discordant read pairs show up
      }
      if(!isDup)
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
    fDepth << p.second << '\t' << p.first << '\n';
}


int main(int argc, char* argv[]){ // <left> <right> <allexons> <end_length> <min_mapping_quality> <details> <depth> 
  std::ios::sync_with_stdio(false); cin.tie(NULL);
  
  fLeft.open(argv[1]);
  fRight.open(argv[2]);
  fExons.open(argv[3]);
  fDetails.open(argv[6]);
  fDepth.open(argv[7]);


  fDetails << "GENE1\tGENE2\tFUSION\tCHR1\tCHR2\tGENE1_START\tGENE2_START\tREAD1_POS\tREAD2_POS\n";
  storeExons(atoi(argv[4]));
  generateCandidates(atoi(argv[5]), atoi(argv[4]), fLeft);
  generateCandidates(atoi(argv[5]), atoi(argv[4]), fRight);

  fDepth << "FUSION\tDEPTH\n";
  printDepths();

  fLeft.close();
  fRight.close();
  fExons.close();
  fDetails.close();
  fDepth.close();
  return 0;
}