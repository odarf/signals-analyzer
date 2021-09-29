#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <QMainWindow>

class analysis
{
public:
    analysis();
    float get_statistics();
    QVector<float> minMax();
    float averageValue();
};

#endif // ANALYSIS_H
