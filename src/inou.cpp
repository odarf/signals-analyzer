#include "inou.h"
#include <stdio.h>
#include <fstream>

inou::inou(){}

QVector<double> inou::load(string path){
    QVector<double> outputData;
    ifstream file(path, ios::in | ios::binary);
    file.seekg(0, ios::end);
    int tam = file.tellg();
    tam /= 4;
    file.seekg(0, ios::beg);
    QVector<float> foo(tam);
    for (int i=0; i<tam; i++) {
        file.read(reinterpret_cast< char * >(&foo[i]), sizeof(float));
    }
    file.close();
    for (int i=0; i<foo.length(); i++) {
        outputData.append(static_cast<double>(foo[i]));
    }
    return outputData;
}
