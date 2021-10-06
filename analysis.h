#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <QMainWindow>

class analysis
{
private:
    QVector<float> min_max;
public:
    analysis();
    float get_statistics();
    QVector<float> minMax();
    float averageValue();
};

#endif // ANALYSIS_H
