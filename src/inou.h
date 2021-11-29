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
};

#endif // INOU_H
