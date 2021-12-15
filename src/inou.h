#ifndef INOU_H
#define INOU_H

#include "mainwindow.h"
#include <stdio.h>

using namespace std;

class inou
{
public:
    inou();
    QVector<double> load(string path);
    QVector<double> loadWave(const char* fileName);
    void exportWave(QVector<double> inputdData, int samples_count, const char* fileToSave, double volume=1);
    QVector<double> createFormant(QVector<double> inputData, const char* fileToSave);

};

#endif // INOU_H
