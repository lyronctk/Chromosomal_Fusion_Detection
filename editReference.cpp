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
const int CHR_SIZES[] = {16571, 249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129845, 51304556, 155270560, 59373566};
// const int CHR_SIZES[] = {459, 498, 451, 481, 51};

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
    if(chrCtr!=-1 && CHR_SIZES[chrCtr]-sequencePos<SEQUENCE_CHARS*2) //getting close to end of chr, where SEQUENCE_CHARS is not constant
      while(true){
        char c;
        fRef.get(c);
        if(c=='A' || c=='T' || c=='C' || c=='G' || c=='a' || c=='t' || c=='c' || c=='g' || c=='N'){
          if(chrCtr == current.chr && sequencePos==current.position){
            fRef.seekg(-1, ios::cur);
            fRef.put(current.base);
            normals.pop();
            if(normals.empty()) return;
            current = normals.front();
          }
          sequencePos++;
        }
        else if(c=='>'){
          fRef.putback(c);
          break;
        }
      }

    if(fRef.peek() == '>'){ //on header line
      fRef.seekg(HEADER_CHARS, ios::cur);
      chrCtr++;
      sequencePos=1;
    }
    else if(current.chr == chrCtr && current.position-sequencePos<SEQUENCE_CHARS-1){ //mutation is on current line or nearing the end of the chromosome (where the SEQUENCE_CHARS is not constant)
      for(int i=0; i<SEQUENCE_CHARS-1; i++){
        if(current.chr == chrCtr && current.position == i+sequencePos){
          fRef.seekg(0, ios::cur);
          fRef.put(current.base);
          normals.pop();
          if(normals.empty()) return;
          current = normals.front();
        }
        else
          fRef.seekg(1, ios::cur);
      }
      fRef.seekg(1, ios::cur);
      sequencePos+=SEQUENCE_CHARS-1;
    }
    else{
      int skip=0;
      if(current.chr == chrCtr)
        skip = current.position; 
      else
        skip = CHR_SIZES[chrCtr];
      skip -= sequencePos; skip -= SEQUENCE_CHARS; skip /= SEQUENCE_CHARS; //subtracting SEQUENCE_CHARS to avoid corner cases
      skip = max(skip, 1);
      fRef.seekg(SEQUENCE_CHARS*skip, ios::cur); //skips multiple lines
      sequencePos+=(SEQUENCE_CHARS-1)*skip;
    }
  }
}


int main(int argc, char *argv[]){ //<ref_genome> <normal_freq>
  std::ios::sync_with_stdio(false); cin.tie(NULL);
  fRef.open(argv[1]);
  fNormal.open(argv[2]);

  // clock_t t;
  // t = clock();

  cout << "----Storing normal mutations..." << '\n';
  storeNormals();
  cout << "----Substituting normal mutations into reference genome..." << '\n';
  editRef();

  // cout << "----Done" << '\n';
  // cout << "------Executed in " << ((float)(clock()-t))/CLOCKS_PER_SEC << " seconds." << '\n';

  fRef.close();
  fNormal.close();
  return 0;
}