#include <iostream>
#include <fstream>
using namespace std;
// clang++ -std=c++11 -stdlib=libc++ splitUnmapped.cpp -o splitUnmapped
// ./splitUnmapped DLBCL021_Tumor.unmapped.fq DLBCL021_Tumor.left.fq DLBCL021_Tumor.right.fq


const int END_LENGTH=20;


ifstream fUnmapped;
ofstream fLeft, fRight;


void split(){
  while(!fUnmapped.eof()){
    string title="", sequence, plus, score;
    fUnmapped >> title >> sequence >> plus >> score;
    if(title=="") continue;

    fLeft << title + '\n' + sequence.substr(0, END_LENGTH) + "\n+\n" + score.substr(0, END_LENGTH) + '\n';
    int endPos = sequence.length()-END_LENGTH;
    fRight << title + '\n' + sequence.substr(endPos) + "\n+\n" + score.substr(endPos) + '\n';
  }
}


int main(int argc, char* argv[]){ // <unmapped> <left> <right>
  std::ios::sync_with_stdio(false); cin.tie(NULL);
  
  fUnmapped.open(argv[1]);
  fLeft.open(argv[2]);
  fRight.open(argv[3]);

  split();

  fUnmapped.close();
  fLeft.close();
  fRight.close();
  return 0;
}