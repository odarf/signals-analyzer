#include "processing.h"
#include "mainwindow.h"

processing::processing(){}


QVector<double> processing::offset(QVector<double> data, double coeff){
    for(int i = 0; i<data.length(); i++){ data[i] += coeff; }
    return data;
}

float processing::shift(float n){
    return MainWindow().randomGenerator(n) * 1000.0;
}

QVector<double> processing::spikes(QVector<double> data){
    int randomValue = 0;
    for(int i=0; i<5; i++){
        random_device rd;
        mt19937 rng(rd());
        uniform_int_distribution<int> uid(0,1000);
        randomValue = uid(rng);
        data[randomValue] = shift(2); //спайки
    }
    return data;
}

QVector<double> processing::aim(QVector<double> data){
    QVector<double> output;
    int coeff = 10;
    double temp = 0;
    for(int i = 0; i<data.length(); i++){
        temp = data[i] + coeff;
        output.append(temp);
    }
    return output;
}

QVector<double> processing::antiShift(QVector<double> data){
    QVector<double> aimData = aim(data);
    QVector<double> output;
    analysis analysis;
    float mean = analysis.mean(aimData);
    double temp = 0;
    for(int i = 0; i<data.length(); i++){
        temp = aimData[i] - mean;
        output.append(temp);
    }
    return output;
}

QVector<double> processing::antiSpike(QVector<double> spikedData){
    QThread::sleep(1);
    auto output = spikedData;
    for(int i = 1; i<output.length()-1; i++){
        if (qAbs(output[i] / output[i-1]) > 30){
            output[i] = (output[i-1] + output[i+1]) / 2.0;
        }
        else{ output[i] = output[i]; }
    }
    return output;
}
