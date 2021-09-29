#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "analysis.h"

#include <iostream>
#include <QLogValueAxis>
#include <QLineSeries>
#include <QValueAxis>
#include <QChart>
#include <QChartView>

#include <math.h>
#include <time.h>

#include <random>

using namespace QtCharts;

// Генератор случайных чисел
float random_generator(float coeff){
    float left = -1 * coeff;
    float right = 1 * coeff;
    static std::default_random_engine engine{std::random_device()()};
    static std::uniform_real_distribution<float> distribution{left, right};
    return distribution(engine);
}

MainWindow::MainWindow(QWidget *parent): QMainWindow(parent), ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{
    std::cout << random_generator(3.0) << std::endl;
    ui->widget->clearGraphs();
    double a = -1;
    double b = 1;
    //double h = 0.01;
    int N = 1000;
    //int len = (b - a) / h + 2;
    QVector<double> x_lin_inc(N), y_lin_inc(N);
    QVector<double> x_lin_dec(N), y_lin_dec(N);
    QVector<double> x_exp_inc(N), y_exp_inc(N);
    QVector<double> x_exp_dec(N), y_exp_dec(N);
    QVector<double> x_random(N), y_random(N);

    int i = 0;
    int beta = 100;
    float alpha = -1;
    for(int X=0; X<N-1; X++){
        x_lin_inc[i] = X;
        x_lin_dec[i] = X;
        x_exp_inc[i] = X;
        x_exp_dec[i] = X;
        x_random[i] = X;
        y_lin_inc[i] = -a * X;
        y_lin_dec[i] = a * X + 1000;
        y_exp_inc[i] = beta * std::exp(alpha * X);
        y_exp_dec[i] = beta * std::exp(-alpha * X);
        y_random[i] = random_generator(3.0);

        i++;
    }

    ui->widget->addGraph();
    ui->widget->graph(0)->setData(x_lin_inc, y_lin_inc);

    ui->widget->xAxis->setLabel("x");
    ui->widget->yAxis->setLabel("y");

    ui->widget->xAxis->setRange(0, N+1);
    ui->widget->yAxis->setRange(0, N+1);
/*
    double minY, maxY = y[0];
    for(int i = 0; i<len; i++){
        if (y[i]<minY) minY = y[i];
        if (y[i]>maxY) maxY = y[i];
    }
*/
    //ui->widget->yAxis->setRange(minY, maxY);

    ui->widget->replot();

    ui->widget_2->addGraph();
    ui->widget_2->graph(0)->setData(x_lin_dec, y_lin_dec);

    ui->widget_2->xAxis->setLabel("x");
    ui->widget_2->yAxis->setLabel("y");

    ui->widget_2->xAxis->setRange(0, N+1);

/*    minY = y[0];
    maxY = y[0];
    for(int i = 0; i<len; i++){
        if (y[i]<minY) minY = y[i];
        if (y[i]>maxY) maxY = y[i];
    }
    ui->widget_2->yAxis->setRange(minY, maxY);
*/
    ui->widget_2->yAxis->setRange(0, 1001);
    ui->widget_2->replot();

// --------------------------------------------------

    ui->widget_3->addGraph();
    ui->widget_3->graph(0)->setData(x_exp_inc, y_exp_inc);

    ui->widget_3->xAxis->setLabel("x");
    ui->widget_3->yAxis->setLabel("y");

    ui->widget_3->xAxis->setRange(0, N+1);

/*    minY = y[0];
    maxY = y[0];
    for(int i = 0; i<len; i++){
        if (y[i]<minY) minY = y[i];
        if (y[i]>maxY) maxY = y[i];
    }
    ui->widget_3->yAxis->setRange(minY, maxY);
*/
    ui->widget_3->yAxis->setRange(0, 1001);
    ui->widget_3->replot();

// --------------------------------------------------

    ui->widget_4->addGraph();
    ui->widget_4->graph(0)->setData(x_exp_dec, y_exp_dec);

    ui->widget_4->xAxis->setLabel("x");
    ui->widget_4->yAxis->setLabel("y");

    ui->widget_4->xAxis->setRange(0, N+1);
/*    minY = y[0];
    maxY = y[0];
    for(int i = 0; i<len; i++){
        if (y[i]<minY) minY = y[i];
        if (y[i]>maxY) maxY = y[i];
    }
    ui->widget_4->yAxis->setRange(minY, maxY);
*/
    ui->widget_4->yAxis->setRange(0, 1001);
    ui->widget_4->replot();

// --------------------------------------------------

    ui->widget_5->addGraph();
    ui->widget_5->graph(0)->setData(x_random, y_random);

    ui->widget_5->xAxis->setLabel("x");
    ui->widget_5->yAxis->setLabel("y");

    ui->widget_5->xAxis->setRange(0, N+1);
/*
    minY = y[0];
    maxY = y[0];
    for(int i = 0; i<len; i++){
        if (y[i]<minY) minY = y[i];
        if (y[i]>maxY) maxY = y[i];
    }
    ui->widget_5->yAxis->setRange(minY, maxY);
*/
    ui->widget_5->yAxis->setRange(-5, 5);
    ui->widget_5->replot();
}



/*void MainWindow::on_pushButton_clicked(){
    QChartView *chrtView = new QChartView(this);

}
*/
