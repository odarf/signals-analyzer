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

const int LENGTH = 1000;

MainWindow::MainWindow(QWidget *parent): QMainWindow(parent), ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

// ---------------- Генератор случайных чисел ----------------
float MainWindow::randomGenerator(float coeff){
    float left = -1.0 * coeff;
    float right = 1.0 * coeff;
    std::default_random_engine engine{std::random_device()()};
    std::uniform_real_distribution<float> distribution{left, right};
    return distribution(engine);
}

void MainWindow::on_pushButton_clicked(){
    //int selectedTask =  ui->cbSelectTask->currentText().toInt();

    randomCoeff = (float)ui->sbCoeff->value();
    std::cout << randomCoeff << std::endl;

    QVector<double> yLinInc(LENGTH), yLinDec(LENGTH),
                    yExpInc(LENGTH), yExpDec(LENGTH),
                    yRandom(LENGTH), x(LENGTH);

    double a = -1;
    double b = 1; // ?
    double beta = 1;
    double alpha = -0.01;

    for(int X=0; X<LENGTH-1; X++){
        x[X] = X;

        yLinInc[X] = -a * X;
        yLinDec[X] = a * X + LENGTH;
        yExpInc[X] = (double)(beta * exp(-alpha * X)); // определить
        yExpDec[X] = (double)(beta * exp(alpha * X));  // границы
        yRandom[X] = randomGenerator(randomCoeff);
    }

// ---------------- Рисование ----------------

    ui->widget->addGraph();
    ui->widget->graph(0)->setData(x, yLinInc);
    ui->widget->xAxis->setRange(0, LENGTH+1);
    ui->widget->yAxis->setRange(0, LENGTH+1);
    ui->widget->replot();

    ui->widget_2->addGraph();
    ui->widget_2->graph(0)->setData(x, yLinDec);
    ui->widget_2->xAxis->setRange(0, LENGTH+1);
    ui->widget_2->yAxis->setRange(0, LENGTH+1);
    ui->widget_2->replot();

    ui->widget_3->addGraph();
    ui->widget_3->graph(0)->setData(x, yExpInc);
    ui->widget_3->xAxis->setRange(0, LENGTH+1);
    ui->widget_3->yAxis->setRange(0, LENGTH+1);
    ui->widget_3->replot();

    ui->widget_4->addGraph();
    ui->widget_4->graph(0)->setData(x, yExpDec);
    ui->widget_4->xAxis->setRange(0, LENGTH+1); // определить границы
    ui->widget_4->yAxis->setRange(0, LENGTH+1); // для у и х
    ui->widget_4->replot();

    ui->widget_5->addGraph();
    ui->widget_5->graph(0)->setData(x, yRandom);
    ui->widget_5->xAxis->setRange(0, LENGTH+1);
    ui->widget_5->yAxis->setRange(-randomCoeff-0.5, randomCoeff+0.5);
    ui->widget_5->replot();

    ui->widget->clearGraphs();
    ui->widget_2->clearGraphs();
    ui->widget_3->clearGraphs();
    ui->widget_4->clearGraphs();
    ui->widget_5->clearGraphs();

// --------------------------------------------------

/*
 *  вычисление границ, может пригодиться
 *  minY = y[0];
 *  maxY = y[0];
 *  for(int i = 0; i<len; i++){
 *      if (y[i]<minY) minY = y[i];
 *      if (y[i]>maxY) maxY = y[i];
 *  }
 *  ui->widget_2->yAxis->setRange(minY, maxY);
 */


}
