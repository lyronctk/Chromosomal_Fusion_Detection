#include <iostream>
#include <fstream>
using namespace std;
// clang++ -std=c++11 -stdlib=libc++ splitUnmapped.cpp -o splitUnmapped
// ./splitUnmapped DLBCL021_Tumor.unmapped.fq 25 DLBCL021_Tumor.left.fq DLBCL021_Tumor.right.fq


ifstream fUnmapped;
ofstream fLeft, fRight;


void split(int endLength){
  while(!fUnmapped.eof()){
    string title="", sequence, plus, score;
    fUnmapped >> title >> sequence >> plus >> score;
    if(title=="") continue;

    fLeft << title + '\n' + sequence.substr(0, endLength) + "\n+\n" + score.substr(0, endLength) + '\n';
    int endPos = sequence.length()-endLength;
    fRight << title + '\n' + sequence.substr(endPos) + "\n+\n" + score.substr(endPos) + '\n';
  }
}


int main(int argc, char* argv[]){ // <unmapped> <end_length> <left> <right> 
  std::ios::sync_with_stdio(false); cin.tie(NULL);
  
  fUnmapped.open(argv[1]);
  fLeft.open(argv[3]);
  fRight.open(argv[4]);

  split(atoi(argv[2]));

  fUnmapped.close();
  fLeft.close();
  fRight.close();
  return 0;
}