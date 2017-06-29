#!/bin/sh
clang++ -std=c++11 -stdlib=libc++ discordantPairs.cpp -o discordantPairs
./discordantPairs DLBCL021_Tumor.remapped.left.txt DLBCL021_Tumor.remapped.right.txt RefSeq_Gencodev17_022314.allexons.noRPnoNA.sorted.bed 40 DLBCL021_Tumor.discordantpairs.details.unsorted.txt DLBCL021_Tumor.discordantpairs.depth.unsorted.txt
(head -1 DLBCL021_Tumor.discordantpairs.details.unsorted.txt ; tail -n +2 DLBCL021_Tumor.discordantpairs.details.unsorted.txt | sort) > DLBCL021_Tumor.discordantpairs.details.txt
(head -1 DLBCL021_Tumor.discordantpairs.depth.unsorted.txt ; tail -n +2 DLBCL021_Tumor.discordantpairs.depth.unsorted.txt | sort -k2,2nr) > DLBCL021_Tumor.discordantpairs.depth.txt
rm DLBCL021_Tumor.discordantpairs.details.unsorted.txt DLBCL021_Tumor.discordantpairs.depth.unsorted.txt