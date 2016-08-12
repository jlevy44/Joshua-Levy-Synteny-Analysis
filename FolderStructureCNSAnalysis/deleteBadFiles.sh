#!bin/bash/

for folder in HDS111 NKH011 NKH101 NKH111 NKHDS01100 NKHDS01111 NKHDS10100 NKHDS10111 NKHDS11100 NKHDS11111; do
    rm -R "$folder"
    mkdir "$folder"
done

