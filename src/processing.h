#ifndef PROCESSING_H
#define PROCESSING_H

#include "analysis.h"
#include <QMainWindow>
#include <iostream>

class processing
{
public:
    processing();
    QVector<double> antiShift(QVector<double>);
    QVector<double> antiSpike();
};

#endif // PROCESSING_H
