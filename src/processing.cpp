#include "processing.h"

processing::processing(){}

QVector<double> processing::antiShift(QVector<double> input){
    QVector<double> output(input.length());
    analysis t_analysis;
    float mean = t_analysis.averageValue(input);
    //auto disp = analysis.dispersion(input);
    std::cout << mean;
    return output;
}
