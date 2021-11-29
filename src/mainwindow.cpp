#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "modeling.h"
#include "analysis.h"
#include "inou.h"
#include "qmath.h"

#include <iostream>
#include <QLogValueAxis>
#include <QLineSeries>
#include <QValueAxis>
#include <QChart>
#include <QChartView>
#include <QtCharts>
#include <QtWidgets>
#include <QtWidgets/QHBoxLayout>
#include <QHBoxLayout>
#include <cstdlib>

#include <math.h>
#include <time.h>
#include <random>

#include <fstream>

using namespace QtCharts;
using namespace std;

const int LENGTH = 1000;
const double DELTA = 50.0;

MainWindow::MainWindow(QWidget *parent): QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);
}

MainWindow::~MainWindow() { delete ui; }

void MainWindow::on_pushButton_clicked(){
    modeling  modeling;
    processing processing;
    analysis analysis;

    double k = -1;

    double beta = 11;
    double alpha = -0.01;
    int amplitude[] = {10, 100, 15};
    int frequency[] = {4, 37, 173};
    float randomCoeff = (float)ui->sbCoeff->value();

    QVector<double> x(LENGTH),             xLarge(3*LENGTH),
                    yLinInc,               yLinDec,                 yExpInc,
                    yExpDec,               yHarmonic,               yPolyharmonic,
                    yCombinated,           yCardio,
                    yMyRandom(LENGTH),     yEmbedRandom1(LENGTH),   yEmbedRandom2(LENGTH);

    yLinInc.append(modeling.linearIncrement(LENGTH, k));
    yLinDec.append(modeling.linearDecrement(LENGTH, k));
    yExpInc.append(modeling.exponentialIncrement(LENGTH));
    yExpDec.append(modeling.exponentialDecrement(LENGTH));
    yHarmonic.append(modeling.harmonic(LENGTH, 10, 4));
    yPolyharmonic.append(modeling.polyharmonic(LENGTH, amplitude, frequency));
    yCardio.append(modeling.cardiogram());

    for(int X=0; X<LENGTH; X++){
        x[X] = X;
        yMyRandom[X] = modeling.randomGenerator(randomCoeff);
        yEmbedRandom1[X] = modeling.embedRandom();
        yEmbedRandom2[X] = modeling.embedRandom();
    }

    for(int i = 0; i<3000; i++){
        xLarge[i] = i;
    }

//------Для первого задания ------
    yCombinated.append(yLinInc);
    yCombinated.append(yLinDec);
    yCombinated.append(yExpInc);
//--------------------------------

    QVector<double> dataFile = inou().load("data.dat");
    QVector<double> xfile(dataFile.length());
    for (int i=0; i<xfile.length(); i++) { xfile[i] = i; }

    ui->graphFromFile->addGraph();
    ui->graphFromFile->graph(0)->setData(xfile, dataFile);
    ui->graphFromFile->xAxis->setRange(0, xfile.length());
    ui->graphFromFile->yAxis->setRange(analysis.minValue(dataFile), analysis.maxValue(dataFile));
    ui->graphFromFile->replot();

    QVector<double> fouramp = analysis.fourierAmplitude(dataFile);
    fouramp = processing.offset(fouramp, 10);

    ui->graphFourierAmplitude->addGraph();
    ui->graphFourierAmplitude->graph(0)->setData(xfile, fouramp);
    ui->graphFourierAmplitude->xAxis->setRange(0, xfile.length());
    ui->graphFourierAmplitude->yAxis->setRange(analysis.minValue(fouramp), analysis.maxValue(fouramp));
    ui->graphFourierAmplitude->replot();

    QVector<double> fourspec = analysis.fourierSpectrum(dataFile, 0.91);

    ui->graphFourierSpectrum->addGraph();
    ui->graphFourierSpectrum->graph(0)->setData(xfile, fourspec);
    ui->graphFourierSpectrum->xAxis->setRange(0, xfile.length());
    ui->graphFourierSpectrum->yAxis->setRange(analysis.minValue(fourspec), analysis.maxValue(fourspec));
    ui->graphFourierSpectrum->replot();


// ---------------- Рисование ----------------
    int currentTask = ui->tabWidget->currentIndex() + 1;
    switch(currentTask){
        case 1:
        // ---------------- Task 1 ----------------
            ui->graphLinInc->addGraph();
            ui->graphLinInc->graph(0)->setData(x, yLinInc);
            ui->graphLinInc->xAxis->setRange(0, LENGTH+1);
            ui->graphLinInc->yAxis->setRange(analysis.minValue(yLinInc), analysis.maxValue(yLinInc)+1);
            ui->graphLinInc->replot();

            ui->graphLinDec->addGraph();
            ui->graphLinDec->graph(0)->setData(x, yLinDec);
            ui->graphLinDec->xAxis->setRange(0, LENGTH+1);
            ui->graphLinDec->yAxis->setRange(analysis.minValue(yLinDec), analysis.maxValue(yLinDec)+1);
            ui->graphLinDec->replot();

            ui->graphExpInc->addGraph();
            ui->graphExpInc->graph(0)->setData(x, yExpInc);
            ui->graphExpInc->xAxis->setRange(0, LENGTH+1);
            ui->graphExpInc->yAxis->setRange(analysis.minValue(yExpInc), analysis.maxValue(yExpInc)+1);
            ui->graphExpInc->replot();

            ui->graphExpDec->addGraph();
            ui->graphExpDec->graph(0)->setData(x, yExpDec);
            ui->graphExpDec->xAxis->setRange(0, LENGTH+1);
            ui->graphExpDec->yAxis->setRange(analysis.minValue(yExpDec), analysis.maxValue(yExpDec));
            ui->graphExpDec->replot();

            ui->graphCombinated->addGraph();
            ui->graphCombinated->graph(0)->setData(xLarge, yCombinated);
            ui->graphCombinated->xAxis->setRange(0, 3*LENGTH);
            ui->graphCombinated->yAxis->setRange(0, 3000);
            ui->graphCombinated->replot();
        // ----------------------------------------
            break;
        case 2:
        // ---------------- Task 2 ----------------
            ui->graphMyRandom->addGraph();
            ui->graphMyRandom->graph(0)->setData(x, yMyRandom);
            ui->graphMyRandom->xAxis->setRange(0, LENGTH+1);
            ui->graphMyRandom->yAxis->setRange(-randomCoeff-0.5, randomCoeff+0.5);
            ui->graphMyRandom->replot();

            ui->graphEmbedRandom->addGraph();
            ui->graphEmbedRandom->graph(0)->setData(x, yEmbedRandom1);
            ui->graphEmbedRandom->xAxis->setRange(0, LENGTH+1);
            ui->graphEmbedRandom->yAxis->setRange(analysis.minValue(yEmbedRandom1), analysis.maxValue(yEmbedRandom1));
            ui->graphEmbedRandom->replot();
        // ----------------------------------------
            break;
        case 3:
        // ---------------- Task 3 ----------------
            ui->graphTaskThree->addGraph();
            ui->graphTaskThree->graph(0)->setData(x, yMyRandom);
            ui->graphTaskThree->xAxis->setRange(0, LENGTH);
            ui->graphTaskThree->yAxis->setRange(analysis.minValue(yMyRandom) - 1, analysis.maxValue(yMyRandom) + 1);
            ui->graphTaskThree->replot();
        //---------- Статистики для Task3-----------
            ui->labelMinimumValue->setText(QString::number(analysis.minValue(yMyRandom)));
            ui->labelMaximumValue->setText(QString::number(analysis.maxValue(yMyRandom)));
            ui->labelMeanValue->setText(QString::number(analysis.mean(yMyRandom)));
            ui->labelDispersionValue->setText(QString::number(analysis.dispersion(yMyRandom)));
            ui->labelStdDevValue->setText(QString::number(analysis.standartDeviation(yMyRandom)));
            ui->labelSqrtDevValue->setText(QString::number(analysis.sqrtDeviation(yMyRandom)));
            ui->labelAssymetryValue->setText(QString::number(analysis.assymetry(yMyRandom)));
            ui->labelAssymetryCoeffValue->setText(QString::number(analysis.assymetryCoeff(yMyRandom)));
            ui->labelExcessValue->setText(QString::number(analysis.excess(yMyRandom)));
            ui->labelCurtosisValue->setText(QString::number(analysis.curtosis(yMyRandom)));
            ui->labelIsStationarity->setText(analysis.isStationary(yMyRandom) ? "Стационарен" : "Не стационарен"); //использовать только для случайных процессов
        // -----------------------------------------
            break;
        case 4:{
        // ---------------- Task 4 ----------------
            ui->graphTaskFourRandom1->addGraph();
            ui->graphTaskFourRandom1->graph(0)->setData(x, yEmbedRandom1);
            ui->graphTaskFourRandom1->xAxis->setRange(0, LENGTH);
            ui->graphTaskFourRandom1->yAxis->setRange(analysis.minValue(yEmbedRandom1) - 1, analysis.maxValue(yEmbedRandom1) + 1);
            ui->graphTaskFourRandom1->replot();

            ui->graphTaskFourRandom2->addGraph();
            ui->graphTaskFourRandom2->graph(0)->setData(x, yEmbedRandom2);
            ui->graphTaskFourRandom2->xAxis->setRange(0, LENGTH);
            ui->graphTaskFourRandom2->yAxis->setRange(analysis.minValue(yEmbedRandom2) - 1, analysis.maxValue(yEmbedRandom2) + 1);
            ui->graphTaskFourRandom2->replot();

            ui->graphCovar->addGraph();
            ui->graphCovar->graph(0)->setData(x, analysis.covariance(yEmbedRandom1, yEmbedRandom2));
            ui->graphCovar->xAxis->setRange(0, LENGTH);
            ui->graphCovar->yAxis->setRange(analysis.minValue(analysis.covariance(yEmbedRandom1, yEmbedRandom2)), analysis.maxValue(analysis.covariance(yEmbedRandom1, yEmbedRandom2)));
            ui->graphCovar->replot();

            ui->graphAutocorr->addGraph();
            ui->graphAutocorr->graph(0)->setData(x, analysis.autocorrelation(yEmbedRandom1));
            ui->graphAutocorr->xAxis->setRange(0, LENGTH);
            ui->graphAutocorr->yAxis->setRange(analysis.minValue(analysis.autocorrelation(yEmbedRandom1)), 0.2);
            ui->graphAutocorr->replot();

            QBarSet *barSet = new QBarSet("");
            QBarSeries *series = new QBarSeries();
            int M = 20;
            QVector<double> hist(20);
            double kHist = analysis.minValue(yEmbedRandom1);
            double shag = (qAbs(analysis.minValue(yEmbedRandom1)) + analysis.maxValue(yEmbedRandom1))/M;

            //0..20, 20..40, 40..60, 60..80, 80..100, ... гистограмма
            for(int i = 0; i<M; i++){
                for(int j = 0; j<yEmbedRandom1.length(); j++){
                    if(yEmbedRandom1[j] >= kHist && yEmbedRandom1[j] <= (kHist+shag)){
                        hist[i] = hist[i] + 1;
                    }
                }
                kHist += shag;
                *barSet << hist[i];
            }

            series->append(barSet);
            QChart *chart = new QChart();
            QValueAxis *axisY = new QValueAxis();
            QValueAxis *axisX = new QValueAxis();
            axisY->setRange(0, round(analysis.maxValue(hist)));
            axisX->setRange(0, analysis.maxValue(yEmbedRandom1));
            axisX->setTickCount(M);
            axisX->setTickInterval(shag);
            chart->addAxis(axisY, Qt::AlignLeft);
            chart->addAxis(axisX, Qt::AlignBottom);
            chart->legend()->setVisible(false);
            series->attachAxis(axisY);
            series->attachAxis(axisX);
            chart->addSeries(series);
            chart->setAnimationOptions(QChart::SeriesAnimations);
            QChartView *chartView = new QChartView(chart);
            chartView->setRenderHint(QPainter::Antialiasing);
            ui->widgetTest->setContentsMargins(0,0,0,0);
            QHBoxLayout *lay = new QHBoxLayout(ui->widgetTest);
            lay->addWidget(chartView);
        // ----------------------------------------
            break;
        }
        case 5:{
        // ---------------- Task 5 ----------------
            ui->graphTaskFiveHarmo->addGraph();
            ui->graphTaskFiveHarmo->graph(0)->setData(x, yHarmonic);
            ui->graphTaskFiveHarmo->xAxis->setRange(0, LENGTH);
            ui->graphTaskFiveHarmo->yAxis->setRange(analysis.minValue(yHarmonic)-1, analysis.maxValue(yHarmonic)+1);
            ui->graphTaskFiveHarmo->replot();

            ui->graphTaskFivePolyharmo->addGraph();
            ui->graphTaskFivePolyharmo->graph(0)->setData(x, yPolyharmonic);
            ui->graphTaskFivePolyharmo->xAxis->setRange(0, LENGTH);
            ui->graphTaskFivePolyharmo->yAxis->setRange(analysis.minValue(yPolyharmonic)-1, analysis.maxValue(yPolyharmonic)+1);
            ui->graphTaskFivePolyharmo->replot();

            ui->graphTaskFiveHarmoAutocorr->addGraph();
            ui->graphTaskFiveHarmoAutocorr->graph(0)->setData(x, analysis.autocorrelation(yHarmonic));
            ui->graphTaskFiveHarmoAutocorr->xAxis->setRange(0, LENGTH);
            ui->graphTaskFiveHarmoAutocorr->yAxis->setRange(analysis.minValue(analysis.autocorrelation(yHarmonic))-1, analysis.maxValue(analysis.autocorrelation(yHarmonic))+1);
            ui->graphTaskFiveHarmoAutocorr->replot();

            ui->graphTaskFivePolyAutocorr->addGraph();
            ui->graphTaskFivePolyAutocorr->graph(0)->setData(x, analysis.autocorrelation(yPolyharmonic));
            ui->graphTaskFivePolyAutocorr->xAxis->setRange(0, LENGTH);
            ui->graphTaskFivePolyAutocorr->yAxis->setRange(analysis.minValue(analysis.autocorrelation(yPolyharmonic))-1, analysis.maxValue(analysis.autocorrelation(yPolyharmonic))+1);
            ui->graphTaskFivePolyAutocorr->replot();

            ui->graphTaskFiveCrosscorr->addGraph();
            ui->graphTaskFiveCrosscorr->graph(0)->setData(x, analysis.covariance(yHarmonic, yPolyharmonic));
            ui->graphTaskFiveCrosscorr->xAxis->setRange(0, LENGTH);
            ui->graphTaskFiveCrosscorr->yAxis->setRange(analysis.minValue(analysis.covariance(yHarmonic, yPolyharmonic))-1, analysis.maxValue(analysis.covariance(yHarmonic, yPolyharmonic))+1);
            ui->graphTaskFiveCrosscorr->replot();

            QBarSet *barSet2 = new QBarSet("");
            QBarSeries *series2 = new QBarSeries();
            int M = 20;
            QVector<double> hist2(20);
            double kHist = analysis.minValue(yHarmonic);
            double shag = (qAbs(analysis.minValue(yHarmonic)) + analysis.maxValue(yHarmonic))/M;

            //0..20, 20..40, 40..60, 60..80, 80..100, ... гистограмма
            for(int i = 0; i<M; i++){
                for(int j = 0; j<yHarmonic.length(); j++){
                    if(yHarmonic[j] >= kHist && yHarmonic[j] <= (kHist+shag)){
                        hist2[i] = hist2[i] + 1;
                    }
                }
                kHist += shag;
                *barSet2 << hist2[i];
            }

            series2->append(barSet2);
            QChart *chart2 = new QChart();
            QValueAxis *axisY2 = new QValueAxis();
            QValueAxis *axisX2 = new QValueAxis();
            axisY2->setRange(0, round(analysis.maxValue(hist2)));
            axisX2->setRange(0, analysis.maxValue(yHarmonic));
            axisX2->setTickCount(M);
            axisX2->setTickInterval(shag);
            chart2->addAxis(axisY2, Qt::AlignLeft);
            chart2->addAxis(axisX2, Qt::AlignBottom);
            chart2->legend()->setVisible(false);
            series2->attachAxis(axisY2);
            series2->attachAxis(axisX2);
            chart2->addSeries(series2);
            chart2->setAnimationOptions(QChart::SeriesAnimations);
            QChartView *chartView2 = new QChartView(chart2);
            chartView2->setRenderHint(QPainter::Antialiasing);
            ui->graphTaskFiveDensity->setContentsMargins(0,0,0,0);
            QHBoxLayout *lay2 = new QHBoxLayout(ui->graphTaskFiveDensity);
            lay2->addWidget(chartView2);
            //QThread::sleep(1);
        // ----------------------------------------
            break;
        }
        case 6:{
        // ---------------- Task 6 ----------------
            int offsetCoeff = 200;
            QVector<double> offsetted = processing.offset(yLinInc, offsetCoeff);

            ui->graphTaskSixOrig1->addGraph();
            ui->graphTaskSixOrig1->graph(0)->setData(x, yLinInc);
            ui->graphTaskSixOrig1->xAxis->setRange(0, LENGTH);
            ui->graphTaskSixOrig1->yAxis->setRange(analysis.minValue(yLinInc)-1, analysis.maxValue(yLinInc)+1);
            ui->graphTaskSixOrig1->replot();

            ui->graphTaskSixShifted1->addGraph();
            ui->graphTaskSixShifted1->graph(0)->setData(x, offsetted);
            ui->graphTaskSixShifted1->xAxis->setRange(0, LENGTH);
            ui->graphTaskSixShifted1->yAxis->setRange(analysis.minValue(offsetted)-offsetCoeff, analysis.maxValue(offsetted)+offsetCoeff);
            ui->graphTaskSixShifted1->replot();

            ui->graphTaskSixOrig2->addGraph();
            ui->graphTaskSixOrig2->graph(0)->setData(x, yHarmonic);
            ui->graphTaskSixOrig2->xAxis->setRange(0, LENGTH);
            ui->graphTaskSixOrig2->yAxis->setRange(analysis.minValue(yHarmonic)-1, analysis.maxValue(yHarmonic)+1);
            ui->graphTaskSixOrig2->replot();

            offsetCoeff = 5;
            offsetted = processing.offset(yHarmonic, 5);

            ui->graphTaskSixShifted2->addGraph();
            ui->graphTaskSixShifted2->graph(0)->setData(x, offsetted);
            ui->graphTaskSixShifted2->xAxis->setRange(0, LENGTH);
            ui->graphTaskSixShifted2->yAxis->setRange(analysis.minValue(offsetted)-offsetCoeff, analysis.maxValue(offsetted)+offsetCoeff);
            ui->graphTaskSixShifted2->replot();

            QVector<double> spiked = processing.spikes(yHarmonic);

            ui->graphTaskSixSpiked->addGraph();
            ui->graphTaskSixSpiked->graph(0)->setData(x, spiked);
            ui->graphTaskSixSpiked->xAxis->setRange(0, LENGTH);
            ui->graphTaskSixSpiked->yAxis->setRange(analysis.minValue(spiked)-1, analysis.maxValue(spiked)+1);
            ui->graphTaskSixSpiked->replot();
        // ----------------------------------------
            break;
        }
        case 7:{
        // ---------------- Task 7 ----------------
            int offsetCoeff = 5;
            QVector<double> offsetted = processing.offset(yHarmonic, offsetCoeff);
            QVector<double> spiked = processing.spikes(yHarmonic);

            ui->graphTaskSevenShifted->addGraph();
            ui->graphTaskSevenShifted->graph(0)->setData(x, offsetted);
            ui->graphTaskSevenShifted->xAxis->setRange(0, LENGTH);
            ui->graphTaskSevenShifted->yAxis->setRange(analysis.minValue(offsetted)-offsetCoeff, analysis.maxValue(offsetted)+offsetCoeff);
            ui->graphTaskSevenShifted->replot();

            QVector<double> antiShifted = processing.antiShift(offsetted);

            ui->graphTaskSevenAntiShifted->addGraph();
            ui->graphTaskSevenAntiShifted->graph(0)->setData(x, antiShifted);
            ui->graphTaskSevenAntiShifted->xAxis->setRange(0, LENGTH);
            ui->graphTaskSevenAntiShifted->yAxis->setRange(analysis.minValue(antiShifted)-offsetCoeff, analysis.maxValue(antiShifted)+offsetCoeff);
            ui->graphTaskSevenAntiShifted->replot();

            ui->graphTaskSevenSpiked->addGraph();
            ui->graphTaskSevenSpiked->graph(0)->setData(x, spiked);
            ui->graphTaskSevenSpiked->xAxis->setRange(0, LENGTH);
            //ui->graphTaskSevenSpiked->yAxis->setRange(analysis.minValue(spiked)-1, analysis.maxValue(spiked)+1);
            ui->graphTaskSevenSpiked->yAxis->setRange(-100, 100);
            ui->graphTaskSevenSpiked->replot();

            QVector<double> antispiked = processing.antiSpike(spiked);

            ui->graphTaskSevenAntiSpiked->addGraph();
            ui->graphTaskSevenAntiSpiked->graph(0)->setData(x, antispiked);
            ui->graphTaskSevenAntiSpiked->xAxis->setRange(0, LENGTH);
            ui->graphTaskSevenAntiSpiked->yAxis->setRange(-100, 100);
            ui->graphTaskSevenAntiSpiked->replot();
        // ----------------------------------------
            break;
        }
        case 8:{
        // ---------------- Task 8 ----------------
            QVector<double> trAdRa = processing.trendAddRandom(yLinInc);
            QVector<double> pot = processing.pickoutTrend(trAdRa);
            QVector<double> antr = processing.antiTrend(trAdRa);

            ui->graphTaskEightTrendAdd->addGraph();
            ui->graphTaskEightTrendAdd->graph(0)->setData(x, trAdRa);
            ui->graphTaskEightTrendAdd->xAxis->setRange(0, LENGTH);
            ui->graphTaskEightTrendAdd->yAxis->setRange(analysis.minValue(trAdRa)-1, analysis.maxValue(trAdRa)+1);
            ui->graphTaskEightTrendAdd->replot();

            ui->graphTaskEightPickTrend->addGraph();
            ui->graphTaskEightPickTrend->graph(0)->setData(x, pot);
            ui->graphTaskEightPickTrend->xAxis->setRange(0, LENGTH);
            ui->graphTaskEightPickTrend->yAxis->setRange(analysis.minValue(pot)-1, analysis.maxValue(pot)+1);
            ui->graphTaskEightPickTrend->replot();

            ui->graphTaskEightAntiTrend->addGraph();
            ui->graphTaskEightAntiTrend->graph(0)->setData(x, antr);
            ui->graphTaskEightAntiTrend->xAxis->setRange(0, LENGTH);
            ui->graphTaskEightAntiTrend->yAxis->setRange(analysis.minValue(antr)-1, analysis.maxValue(antr)+1);
            ui->graphTaskEightAntiTrend->replot();

            QVector<double> stdDev = processing.pdfTaskEight();

            ui->graphTaskEightStdDev->addGraph();
            ui->graphTaskEightStdDev->graph(0)->setData(x, stdDev);
            ui->graphTaskEightStdDev->xAxis->setRange(0, stdDev.length()-1);
            ui->graphTaskEightStdDev->yAxis->setRange(0, analysis.maxValue(stdDev)+0.05);
            ui->graphTaskEightStdDev->replot();

            QVector<double> stdDevDiv = processing.pdfTaskEight2(stdDev);

            ui->graphTaskEightStdDev_2->addGraph();
            ui->graphTaskEightStdDev_2->graph(0)->setData(x, stdDevDiv);
            ui->graphTaskEightStdDev_2->xAxis->setRange(0, stdDevDiv.length()-1);
            ui->graphTaskEightStdDev_2->yAxis->setRange(analysis.minValue(stdDevDiv)-0.5, analysis.maxValue(stdDevDiv)+0.5);
            ui->graphTaskEightStdDev_2->replot();
        // ----------------------------------------
            break;
        }
        case 9:{
        // ---------------- Task 9 ----------------
            QVector<double> taskNine1 = processing.pdfTaskNine(yHarmonic, 1);
            QVector<double> taskNine10 = processing.pdfTaskNine(yHarmonic, 10);
            QVector<double> taskNine100 = processing.pdfTaskNine(yHarmonic, 100);
            QVector<double> taskNine1000 = processing.pdfTaskNine(yHarmonic, 1000);


            ui->graphTaskNine1->addGraph();
            ui->graphTaskNine1->graph(0)->setData(x, taskNine1);
            ui->graphTaskNine1->xAxis->setRange(0, taskNine1.length()+1);
            ui->graphTaskNine1->yAxis->setRange(analysis.minValue(taskNine1), analysis.maxValue(taskNine1));
            ui->graphTaskNine1->replot();

            ui->graphTaskNine10->addGraph();
            ui->graphTaskNine10->graph(0)->setData(x, taskNine10);
            ui->graphTaskNine10->xAxis->setRange(0, taskNine10.length()+1);
            ui->graphTaskNine10->yAxis->setRange(analysis.minValue(taskNine10), analysis.maxValue(taskNine10));
            ui->graphTaskNine10->replot();

            ui->graphTaskNine100->addGraph();
            ui->graphTaskNine100->graph(0)->setData(x, taskNine100);
            ui->graphTaskNine100->xAxis->setRange(0, taskNine100.length()+1);
            ui->graphTaskNine100->yAxis->setRange(analysis.minValue(taskNine100), analysis.maxValue(taskNine100));
            ui->graphTaskNine100->replot();

            ui->graphTaskNine1000->addGraph();
            ui->graphTaskNine1000->graph(0)->setData(x, taskNine1000);
            ui->graphTaskNine1000->xAxis->setRange(0, taskNine1000.length()+1);
            ui->graphTaskNine1000->yAxis->setRange(analysis.minValue(taskNine1000), analysis.maxValue(taskNine1000));
            ui->graphTaskNine1000->replot();
        // ----------------------------------------
            break;
        }
        case 10:
        // ---------- Cardio ----------
            ui->graphCardio->addGraph();
            ui->graphCardio->graph(0)->setData(x, yCardio);
            ui->graphCardio->xAxis->setRange(0, 1200);
            ui->graphCardio->yAxis->setRange(analysis.minValue(yCardio)-1, analysis.maxValue(yCardio)+1);
            ui->graphCardio->replot();
        // -----------------------------
            break;

        default:
            return;
    }

// ---------- Очистка графиков ----------
    ui->graphLinInc->clearGraphs();
    ui->graphLinDec->clearGraphs();
    ui->graphExpInc->clearGraphs();
    ui->graphExpDec->clearGraphs();
    ui->graphCombinated->clearGraphs();
    ui->graphMyRandom->clearGraphs();
    ui->graphEmbedRandom->clearGraphs();
    ui->graphTaskThree->clearGraphs();
    ui->graphTaskFourRandom1->clearGraphs();
    ui->graphTaskFourRandom2->clearGraphs();
    ui->graphCovar->clearGraphs();
    ui->graphAutocorr->clearGraphs();
    ui->graphTaskFiveHarmo->clearGraphs();
    ui->graphTaskFivePolyharmo->clearGraphs();
    ui->graphTaskFiveHarmoAutocorr->clearGraphs();
    ui->graphTaskFivePolyAutocorr->clearGraphs();
    ui->graphTaskFiveCrosscorr->clearGraphs();
    ui->graphTaskSixOrig1->clearGraphs();
    ui->graphTaskSixShifted1->clearGraphs();
    ui->graphTaskSixOrig2->clearGraphs();
    ui->graphTaskSixShifted2->clearGraphs();
    ui->graphTaskSixSpiked->clearGraphs();
    ui->graphTaskSevenShifted->clearGraphs();
    ui->graphTaskSevenAntiShifted->clearGraphs();
    ui->graphTaskSevenSpiked->clearGraphs();
    ui->graphTaskSevenAntiSpiked->clearGraphs();
    ui->graphTaskSevenShifted->clearGraphs();
// --------------------------------------
}
