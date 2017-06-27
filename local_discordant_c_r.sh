clang++ -std=c++11 -stdlib=libc++ discordantPairs.cpp -o discordantPairs
./discordantPairs DLBCL021_Tumor.remapped.left.txt DLBCL021_Tumor.remapped.right.txt RefSeq_Gencodev17_022314.allexons.noRPnoNA.sorted.bed 30 40 DLBCL021_Tumor.discordantpairs.details.unsorted.txt DLBCL021_Tumor.discordantpairs.depth.unsorted.txt
(head -1 DLBCL021_Tumor.discordantpairs.details.unsorted.txt ; tail -n +2 DLBCL021_Tumor.discordantpairs.details.unsorted.txt | sort) > A
(head -1 DLBCL021_Tumor.discordantpairs.depth.unsorted.txt ; tail -n +2 DLBCL021_Tumor.discordantpairs.depth.unsorted.txt | sort -k1,1nr) > B