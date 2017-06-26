#include <iostream>
#include <fstream>
#include <string>
#include <unordered_set>
#include <assert.h>
#include <set>
using namespace std;
// clang++ -std=c++11 -stdlib=libc++ findEndMutations.cpp -o findEndMutations
// ./findEndMutations DLBCL021_Tumor.mutatedrows.sorted.txt 25 DLBCL021_Tumor.deduped.sam DLBCL021_Tumor.left.fq DLBCL021_Tumor.right.fq


struct Read{
  string header="";
  int length, startPos;
  set<int> mutations;
};


ofstream fLeft, fRight;
ifstream fRows, fDeduped;


Read reads[2];
bool isReadPaired(string &current, string &before){
  if(current.size() != before.size())
    return false;
  for(int i=0; i<current.size(); i++)
    if(current[i] != before[i] && (before[i] != '1' || current[i] != '2')) 
      return false;
  return true;
}


bool containsEndMutations(set<int> &mutations, int start, int end, int endLength){
  bool left=false, right=false;
  for(int n : mutations){
    if(n < start+endLength)
      left = true;
    else if(n > end-endLength)
      right = true;
  }
  return left && right;
}


unordered_set<string> endMutations;
void storeEndMutations(int endLength){
  int readCtr = 0;
  while(!fRows.eof()){ 
    Read &current = reads[readCtr], &before = reads[(readCtr+1)%2];
    current.mutations.clear();

    string L;
    fRows >> current.header >> current.startPos >> L;
    if(L == "") continue;
    L.pop_back(); current.length = stoi(L);
    
    assert(current.length > endLength*2 && "Error: Read length is too short, input longer reads or change the constant endLength in the code");

    char c;
    int readPos = -1;
    while(true){
      fRows.get(c);
      if(c == '\n')
        break;
      if(c == 'A' || c == 'C' || c == 'T' || c=='G')
        current.mutations.insert(current.startPos+readPos);
      readPos++;
    }

    if(isReadPaired(current.header, before.header)){
      current.mutations.insert(before.mutations.begin(), before.mutations.end());
      if(containsEndMutations(current.mutations, before.startPos, current.startPos+current.length-1, endLength))
        endMutations.emplace(current.header);
      current.mutations.clear();
    }
    else if(!before.mutations.empty() && containsEndMutations(before.mutations, before.startPos, before.startPos+before.length-1, endLength))
      endMutations.emplace(before.header);

    readCtr = (readCtr+1)%2;
  }

  // repeated code for final read in file 
  Read &before = reads[(readCtr+1)%2];
  if(!before.mutations.empty() && containsEndMutations(before.mutations, before.startPos, before.startPos+before.length-1, endLength))
    endMutations.emplace(before.header);
}


void printModifiedHeader(string &header, ofstream &out){
  out.put('@');
  for(int i=0; i<header.size(); i++)
    if(i < 8 || i > 11)
      out.put(header[i]);
  out.put('\n');
}


void findInDeduped(int endLength){
  while(!fDeduped.eof()){
    string header="", filler, sequence, score;
    fDeduped >> header >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> sequence >> score;
    if(header == "") continue;

    if(endMutations.find(header) != endMutations.end()){
      printModifiedHeader(header, fLeft);
      fLeft << sequence.substr(0, endLength) + "\n+\n" + score.substr(0, endLength) + '\n';
      int endPos = sequence.length()-endLength;
      printModifiedHeader(header, fRight);
      fRight << sequence.substr(endPos) + "\n+\n" + score.substr(endPos) + '\n';
    }
  }
}


int main(int argc, char* argv[]){ // <mutated_rows> <end_length> <deduped_bam> <left> <right>
  std::ios::sync_with_stdio(false); cin.tie(NULL);
  
  fRows.open(argv[1]);
  fDeduped.open(argv[3]);
  fLeft.open(argv[4], ios_base::app);
  fRight.open(argv[5], ios_base::app);

  // clock_t t;
  // t = clock();

  storeEndMutations(atoi(argv[2]));
  findInDeduped(atoi(argv[2]));

  // cout << "----Done" << '\n';
  // cout << "------Executed in " << ((float)(clock()-t))/CLOCKS_PER_SEC << " seconds." << '\n';

  fRows.close();
  fDeduped.close();
  fLeft.close();
  fRight.close();
  return 0;
}