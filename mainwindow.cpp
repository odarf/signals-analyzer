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
const double DELTA = 50.0;

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

float MainWindow::fuckUp(float n){
    return randomGenerator(n) * 1000.0;
}

void MainWindow::on_pushButton_clicked(){
    //int selectedTask =  ui->cbSelectTask->currentText().toInt();

    randomCoeff = (float)ui->sbCoeff->value();
    std::cout << randomCoeff << std::endl;

    QVector<double> yLinInc(LENGTH), yLinDec(LENGTH),
                    yExpInc(LENGTH), yExpDec(LENGTH),
                    yRandom(LENGTH), ySin(3*LENGTH),
                    x(LENGTH), xSin(LENGTH);

    double a = -1;
    double b = 1; // ?
    double beta = 1;
    double alpha = -0.01;
    int amplitude[] = {10, 100, 15};
    int frequency[] = {4, 37, 173};
    float delta_t = 0.001;

    for(int X=0; X<LENGTH-1; X++){
        x[X] = X;

        yLinInc[X] = -a * X;
        yLinDec[X] = a * X + LENGTH;
        yExpInc[X] = (double)(beta * exp(-alpha * X)); // определить
        yExpDec[X] = (double)(beta * exp(alpha * X));  // границы
        yRandom[X] = randomGenerator(randomCoeff);
        //ySin[X] = a_first * sin(2 * 3.14 * f_first * X * delta_t);
    }
    float fUpCoeff = 0;
    for(int count=0; count<3; count++){
        for(int i=0; i<5; i++){
            fUpCoeff = fuckUp(1);
            for(int X=0; X<LENGTH-1; X++){
                xSin[X] = X;
                ySin[X] += amplitude[0] * sin(2 * 3.14 * frequency[0] * X * delta_t) + DELTA;
            }
        }
    }

    for(int i=0; i<5; i++){
        std::random_device rd;
        std::mt19937 rng(rd());
        std::uniform_int_distribution<int> uni(0,1000);
        int randomValue = uni(rng);
        ySin[randomValue] = fuckUp(2);
    }

    ui->widget_6->addGraph();
    ui->widget_6->graph(0)->setData(xSin, ySin);
    ui->widget_6->xAxis->setRange(0, LENGTH);
    ui->widget_6->yAxis->setRange(-2500, 2500);
    ui->widget_6->replot();
    QThread::sleep(1);

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
    ui->widget_6->clearGraphs();

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
