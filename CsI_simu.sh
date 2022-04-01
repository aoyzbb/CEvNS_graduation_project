rm -rf build
cmake -B build
cmake --build build --target all -- -j $(nproc) 
cp ./CsI_Cs134.card ./build
cp ./CsI_Cs134.mac ./build
cd build
./nuHunter CsI_Cs134.card