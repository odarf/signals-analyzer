#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "processing.h"
#include "analysis.h"
#include <QMainWindow>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    float randomCoeff = 0.0;
    float randomGenerator(float coeff);
    float embedRandom();
    float shift(float n);
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_pushButton_clicked();

private:
    Ui::MainWindow *ui;
};
#endif // MAINWINDOW_H