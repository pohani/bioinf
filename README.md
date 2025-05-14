# bioinf

get to this folder.



g++ -std=c++17 hirgc_single.cpp -O3 -o hirgc

./hirgc.exe compress ./data/ref.fa ./data/target.fa ./data/out.bin

./hirgc.exe decompress ./data/ref.fa ./data/out.bin ./data/reconstructed.fa
